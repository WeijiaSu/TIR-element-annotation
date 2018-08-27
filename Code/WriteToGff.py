import pandas as pd
import os



def split(s):
    return s.split("_")[1]+"_"+s.split("_")[2]+"_"+s.split("_")[3]


def FlankingBlast(file):
    f=pd.read_table(file,header=None,sep="\t",comment="#")
    f["coor"]=f[0].apply(lambda x: split(x))
    coor=list(f["coor"])
    return list(set(coor))

def TArepeats(s):
    t=s.upper().count("T")
    a=s.upper().count("A")
    ta=t+a
    if (ta==len(s)):
        return True
    else:
        return False

def tir(s):
    return s.split("_")[0].split(":")[1]

def tsd(s):
    return s.split("_")[3].split(":")[1]

def RemovePutativeDTA(file,output):
    f=pd.read_table(file,header=None,sep="\t")
    f["tir"]=f[8].apply(lambda x:tir(x))
    print(f.shape)
    n=f.shape[0]
    f=f.loc[(f["tir"]!="CGGACAGTCCG") & (f["tir"]!="CGGACTGTCCG") & (f["tir"]!="GGACAGTCCGG") & (f["tir"]!="GGACTGTCCGG")]
    print(f.shape)
    print(f.shape[0]/n)
    f.to_csv(output,header=None,index=None,sep="\t")

############################################################################## Remove Entries with Homology in Flanking Sequences ##################################
print("############################################################################## Remove Putative DTA ##################################")
for genome in ["MaizeB73","MaizeMo17","MaizeCML247","Sorghum","Teosinte"]:
    for dataset in ["IRFHB","RefHB", "IRFML"]:
        file = RemovePutativeDTA("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/%s/%s_%s.gff3"%(genome,dataset,genome,dataset),
                                 "/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/%s/%s_%s.gff3" % (genome, dataset, genome, dataset) )
        print("Removing Putative DTA in %s"%(genome))
print("########################################################################## Finished : Remove Putative DTA ##############################")

######################################################################################################################################################################


def removeIRFhomo(file,removelist,outputname):
    f=pd.read_table(file,header=None,sep="\t")
    print("%s entires in total"%(str(f.shape[0])))
    f["coor"]=f[0].apply(str)+"_"+f[3].apply(str)+"_"+f[4].apply(str)
    keep=f.loc[(~(f["coor"].isin(removelist)))]
    keep.to_csv(outputname,header=None,index=None,sep="\t")
    return keep

############################################################################## Remove Entries with Homology in Flanking Sequences ##################################
print("############################################################################## Remove Entries with Homology in Flanking Sequences ##################################")
for genome in ["MaizeB73","MaizeMo17","MaizeCML247","Sorghum","Teosinte"]:
    for dataset in ["IRFHB","RefHB", "IRFML"]:
        blastname=genome+"_"+dataset+"_flankBlast"
        remove = FlankingBlast("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/%s/%s"%(genome,dataset,blastname))
        keep=removeIRFhomo("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/%s/%s_%s.gff3"%(genome,dataset,genome,dataset),remove,
                           "/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/%s/%s_%s_flankchecked.gff3"%(genome,dataset,genome,dataset))
        print("%s removed in %s" % (str(len(remove)), genome))
        print("%s retained in %s"%(str(keep.shape[0]),genome))
        print("copying files to the combine folder in %s"%(genome))
        cp ="cp %s %s"%("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/%s/%s_%s_flankchecked.gff3"%(genome,dataset,genome,dataset),
                        "/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/combine/"%(genome))
        os.system(cp)
print("########################################################################## Finished : Remove Entries with Homology in Flanking Sequences ##############################")

######################################################################################################################################################################


def renameSource(file):
    f = pd.read_table(file, header=None, sep="\t")
    f[1]="IRF_HC"
    f.to_csv(file,header=None,index=None,sep="\t")
    return f

def removeDupinSingle(file):
    f=pd.read_table(file,header=None,sep="\t")
    f=f.sort_values([0,3,4],ascending=[True,True,True])
    f=f.drop_duplicates([0,3,4],keep="last")
    return f

def RemoveTA(f):
    f["TIR"]=f[8].apply(lambda x:tir(x))
    f["TA"]=f["TIR"].apply(lambda x:TArepeats(x))
    sub=f.loc[f["TA"]==False]
    return sub[[0,1,2,3,4,5,6,7,8]]

def combineHBIRFhb(f1,f2):
    f=f1.append(f2,ignore_index=True)
    f=f.sort_values([0,3,4],ascending=[True,True,True])
    f = f.drop_duplicates([0, 3, 4], keep="first")
    return f

def combineAll(f1,f2,out):
    f = f1.append(f2, ignore_index=True)
    f.to_csv(out,header=None,index=None,sep="\t")
    return f

############################################################################## Process and Combining three gff files ################################################
print("############################################################################## Processing and Combining three gff files ##################################")
for genome in ["MaizeB73","MaizeMo17","MaizeCML247","Sorghum","Teosinte"]:
    os.chdir("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/combine"%(genome))
    print("In Genome %s"%(genome))
    print("Renaming IRFHB in %s"%(genome))
    renameSource("%s_IRFHB_flankchecked.gff3"%(genome))
    print("Removing duplicates in %s_HB"%(genome))
    f_hb=removeDupinSingle("%s_RefHB_flankchecked.gff3"%(genome))
    print("removing duplicates in %s_IRFHC"%(genome))
    f_irfref=removeDupinSingle("%s_IRFHB_flankchecked.gff3"%(genome))
    print("removing duplicates in %s_SB"%(genome))
    f_sb=removeDupinSingle("%s_IRFML_flankchecked.gff3"%(genome))
    print("Removing TArepeats in %s_IRFHC"%(genome))
    f_irfref=RemoveTA(f_irfref)
    print("Removing TArepeats in %s_SB"%(genome))
    f_sb = RemoveTA(f_sb)
    print("combine HB and IRFHC in %s"%(genome))
    comHBIRF=combineHBIRFhb(f_hb,f_irfref)
    print("combine three datasets in %s"%(genome))
    comAll=combineAll(comHBIRF, f_sb, "%s_combined.gff3"%(genome))
print("#######################################################################Finished: Processing and Combining three gff files ##################################")
######################################################################################################################################################################

def ProcessGff(file,output):
    f=pd.read_table(file,header=None,sep="\t")
    f["pri"]=0
    mask=f[2]=="DTM"
    f.loc[mask, "pri"] = 1
    mask=f[2]=="DTC"
    f.loc[mask, "pri"] = 2
    mask = f[2] == "DTA"
    f.loc[mask, "pri"] = 3
    mask = f[2] == "DTT"
    f.loc[mask, "pri"] = 4
    mask = f[2] == "DTH"
    f.loc[mask, "pri"] = 5
    f["copy3"] = f[3]
    f["copy4"] = f[4]
    f.copy3 = f.copy3.shift(1).fillna(value=0).astype("int64")
    f.copy4 = f.copy4.shift(1).fillna(value=0).astype("int64")
    mask = ((f[3] == f["copy3"]) & (f[4] == f["copy4"]))
    f.loc[mask, 1] = "Both"
    f = f.sort_values([0, 3, 4, "pri",1], ascending=[True,True, True, True, True])
    f = f.drop_duplicates([0, 3, 4,"pri"], keep="first")
    f = f.sort_values([0, 3, 4,1, "pri"], ascending=[True,True,True, True, True])
    f=f.drop_duplicates([0,3,4],keep="first")
    f["length"]=f[4]-f[3]+1
    f=f[[0,1,2,3,4,5,6,7,8,"pri","length"]]
    f.to_csv(output,header=None,index=None,sep="\t")

############################################################################## Preparing for Removing Overlaps ################################################
print("############################################################################## Preparing for Removing Overlaps ##################################")
for genome in ["MaizeB73","MaizeMo17","MaizeCML247","Sorghum","Teosinte"]:
    os.chdir("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/combine"%(genome))
    print("In Genome %s"%(genome))
    ProcessGff("%s_combined.gff3"%(genome), "%s_combine_all_process.gff3"%(genome))
    f = pd.read_table("%s_combine_all_process.gff3"%(genome), header=None, sep="\t")
    for i in range(1,11):
        chr=f.loc[f[0]==i]
        chr.to_csv("%s_combine_all_process_%s.gff3"%(genome,i),header=None,index=None,sep="\t")
print("############################################################################## Finished: Removing Overlaps ##################################")
######################################################################################################################################################################

def splitInfor(s,i):
    l=s.split("_")
    return float(l[i])


def ProcessSelect(file):
    f = pd.read_table(file, header=None, sep="\t")
    f = f.sort_values([0,3, 4, 9,10], ascending=[True, True, True, True,True])

    f["copy3"]=f[3]
    f["copy4"]=f[4]
    f["copy10"]=f[10]
    f["TIRp"] = f[8].apply(lambda x: splitInfor(x, 2))
    f["TSDp"] = f[8].apply(lambda x: splitInfor(x, 5))
    f["copy9"] = f[9]
    f["copyTIRp"] = f["TIRp"]
    f["copyTSDp"] = f["TSDp"]
    f.copyTIRp = f.TIRp.shift(1).fillna(value=0)
    f.copyTSDp = f.TSDp.shift(1).fillna(value=0)
    f.copy3=f.copy3.shift(1).fillna(value=0).astype("int64")
    f.copy4 = f.copy4.shift(1).fillna(value=0).astype("int64")
    f.copy10 = f.copy10.shift(1).fillna(value=0).astype("int64")
    f.copy9 = f.copy9.shift(1).fillna(value=0).astype("int64")
    f.copyTIRp = f.TIRp.shift(1).fillna(value=0)
    f.copyTSDp = f.TSDp.shift(1).fillna(value=0)
    f["3_copy4"]=f[3]-f["copy4"]
    f["4_copy4"]=f[4]-f["copy4"]
    f["len-len"]=f[10]-f["copy10"]
    f["pri-pri"]=f[9]-f["copy9"]
    f["tir-tir"]=f["TIRp"]-f["copyTIRp"]
    f["tsd-tsd"]=f["TSDp"]-f["copyTSDp"]
    f=f.sort_values([0,3,4,1],ascending=[True,True,True,True])
    return f



def getRemoveList(f):
    removeList=[]
    overlap=(((f[3]-f["copy3"]<=30) & (f[3]-f["copy3"]>=0)) | ((f["copy4"]-f[4]<=30) & (f["copy4"]-f[4]>=0) ))
    r1=f.loc[(f["3_copy4"]<=0)&(f["4_copy4"]<=0) & (overlap==True)]
    if r1.shape[0]==0:
        pass
    removeList=removeList+list(r1.index.values)
    r2=f.loc[(f["3_copy4"]<=0)&(f["4_copy4"]>=0)]
    if r2.shape[0]==0:
        pass
    for index, row in r2.iterrows():
        if(row["pri-pri"]>0):
            removeList.append(index)
        elif(row["pri-pri"]<0):
            removeList.append(index-1)
        else:
            if(row["tir-tir"]<0):
                removeList.append(index)
            elif(row["tir-tir"]>0):
                removeList.append(index-1)
            else:
                if (row["tsd-tsd"] < 0):
                    removeList.append(index)
                elif (row["tsd-tsd"] > 0):
                    removeList.append(index-1)
                else:
                    if(row["len-len"] <= 0):
                        removeList.append(index)
                    else:
                        removeList.append(index - 1)
    return removeList


def CheckOverlap(file):
    re_open = pd.read_table(file, header=None, sep="\t")
    re_open["copy3"] = re_open[3]
    re_open["copy4"] = re_open[4]
    re_open["copy10"] = re_open[10]
    re_open.copy3 = re_open.copy3.shift(1).fillna(value=0).astype("int64")
    re_open.copy4 = re_open.copy4.shift(1).fillna(value=0).astype("int64")
    re_open.copy10 = re_open.copy10.shift(1).fillna(value=0).astype("int64")
    re_open["3_copy4"] = re_open[3] - re_open["copy4"]
    re_open["4_copy4"] = re_open[4] - re_open["copy4"]
    re_open["len-len"] = re_open[10] - re_open["copy10"]
    overlap = (((re_open[3]-re_open["copy3"]<=30) & (re_open[3]-re_open["copy3"]>=0)) | ((re_open["copy4"]-re_open[4]<=30) & (re_open["copy4"]-re_open[4]>=0) ))
    r1 = re_open.loc[(re_open["3_copy4"] <= 0) & (re_open["4_copy4"] <= 0) & (overlap==True)]
    r2 = re_open.loc[(re_open["3_copy4"] <= 0) & (re_open["4_copy4"] >= 0)]
    if (r1.shape[0]!=0 or r2.shape[0]!=0):
        return False
    else:
        return True

def deleteOverlap(file,output):
    f=ProcessSelect(file)
    l=getRemoveList(f)
    if len(l)!=0:
        newf = f.drop(f.index[l])
        newf = newf[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
        newf.to_csv(output, header=None, index=None, sep="\t")
        deleteOverlap(output, output)
    else:
        newf = f[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
        newf.to_csv(output, header=None, index=None, sep="\t")
        return newf

############################################################################## Removing Overlaps ################################################
print("############################################################################## Removing Overlaps ##################################")

for genome in ["MaizeB73","MaizeMo17","MaizeCML247","Sorghum","Teosinte"]:
    os.chdir("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/combine"%(genome))
    print("In Genome %s"%(genome))
    for i in range(1, 11):
        newf = deleteOverlap("%s_combine_all_process_%s.gff3"%(genome,i), "%s_combine_all_process_%s_Oremoved.gff3.txt" % (genome,i))
    cat="cat *_combine_all_process_*_Oremoved.gff3.txt > %s_FinalAnn.gff3"%(genome)
    os.system(cat)
    final=pd.read_table("%s_FinalAnn.gff3"%(genome),header=0,sep=",")
    print("%s Final Size:"%(genome))
    print(final.shape)
print("############################################################################## Finished: Removing Overlaps ##################################")
######################################################################################################################################################################


def RemoveOerlap(file,outname):
    removeList=[]
    f=pd.read_table(file,header=None,sep="\t")
    f = f.sort_values([0, 3, 4, 9, 10], ascending=[True, True, True, True, True])

    f["copy3"] = f[3]
    f["copy4"] = f[4]
    f["copy10"] = f[10]
    f.copy3 = f.copy3.shift(1).fillna(value=0).astype("int64")
    f.copy4 = f.copy4.shift(1).fillna(value=0).astype("int64")
    f.copy10 = f.copy10.shift(1).fillna(value=0).astype("int64")
    f["3_copy4"] = f[3] - f["copy4"]
    f["4_copy4"] = f[4] - f["copy4"]
    f["len-len"] = f[10] - f["copy10"]
    f["copy_0"]=f[0]
    f["0_copy0"]=f[0]-f["copy_0"]

    f = f.sort_values([0, 3, 4, 1], ascending=[True, True, True, True])
    sub1=f.loc[(f["3_copy4"]<=0)&(f["len-len"]<=0)&(f["0_copy0"]==0)]
    removeList = removeList + list(sub1.index.values)
    sub2=f.loc[(f["3_copy4"]<=0)&(f["len-len"]>0)&(f["0_copy0"]==0)]
    removeList = removeList + [i-1 for i in list(sub2.index.values)]
    f=f.drop(f.index[removeList])
    f=f.drop(["copy3","copy4","copy10","3_copy4","4_copy4","len-len"],axis=1)
    f.to_csv(outname,header=None,index=None,sep="\t")
    return removeList

############################################################################## Deleting internal copies ################################################
print("############################################################################## Deleting internal copies ##################################")
for genome in ["MaizeB73","MaizeMo17","MaizeCML247","Sorghum","Teosinte"]:
    os.chdir("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/MakeTE/Result0710/%s/combine"%(genome))
    print("In Genome %s"%(genome))
    RemoveOerlap("%s_FinalAnn.gff3"%(genome), "%s_FinalAnn_forPlot.gff3"%(genome))
    f=pd.read_table("%s_FinalAnn_forPlot.gff3"%(genome),header=None,sep="\t")
    print("%s DistinctforPlot shape:"%(genome))
    print(f.shape)

print("############################################################################## Finished: Deleting internal copies ##################################")
######################################################################################################################################################################
