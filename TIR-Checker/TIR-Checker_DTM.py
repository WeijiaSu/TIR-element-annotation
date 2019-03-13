import sys
import os
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import os
from os.path import basename
import multiprocessing
from multiprocessing import Pool
import time

os.chdir("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/TE_ML/SorghumData/predition")

def compare(tir1,tir2):
    d=0
    for i in range(0,len(tir1)):
        if(tir1[i]!=tir2[i]):
            d+=1
    return d


TIR={}
TIR["DTA"]=11
TIR["DTT"]=13
TIR["DTM"]=40
TIR["DTH"]=14
TIR["DTC"]=13
#
TSD={}
TSD["DTA"]=8
TSD["DTT"]=2
TSD["DTM"]=9
TSD["DTH"]=3
TSD["DTC"]=3

def slidingWindow(seq1,seq2,tsdlength):
    set1=[]
    set2=[]
    for i in range(0,len(seq1)-tsdlength+1):
        set1.append(seq1[i:i+tsdlength])
        set2.append(seq2[i:i+tsdlength])
    return set1,set2


def GetDiff(set1,set2):
    tsd_diff={}
    for i in range(0,len(set1)):
        for j in range(0,len(set2)):
            name=str(i)+":"+str(j)
            diff=compare(set1[i],set2[j])
            tsd_diff[name]=diff
    return tsd_diff

def isTSD(tsd_dffset,fam):
    for i in tsd_dffset:
        if tsd_dffset[i]<TSD[fam]*0.2:
            return True
    return False


# Use CheckTIR1 when seq with TSD
def CheckTIR1(rec):
    dic=[]
    withTIR=[]
    FinalTIRlist=[]
    s = str(rec.seq)[9:-9]
    l=40
    s1 = s[0:l]
    s2_ = s[-l:]
    s2 = Seq(s2_).reverse_complement()
    d = compare(s1, s2)
    if d > l*0.2:
        dic.append(str(rec.id))
    else:
        withTIR.append(str(rec.id))
    FinalTIRlist.append(dic)
    FinalTIRlist.append(withTIR)
    return FinalTIRlist

# Use CheckTIR2 when seq with TSD

def CheckTIR2(rec):
    dic=[]
    withTIR=[]
    FinalTIRlist=[]
    s = str(rec.seq)
    l=40
    s1 = s[0:l]
    s2_ = s[-l:]
    s2 = Seq(s2_).reverse_complement()
    d = compare(s1, s2)
    if d > l*0.2:
        dic.append(str(rec.id))
    else:
        withTIR.append(str(rec.id))
    FinalTIRlist.append(dic)
    FinalTIRlist.append(withTIR)
    return FinalTIRlist


def CheckTSD(rec):
    noTSD=[]
    withTSD=[]
    FinalTSDlist=[]
    s = str(rec.seq)
    l=9
    s1 = s[0:9]
    last20=s[-9:]
    s2 = last20[0:9]
    set1,set2=slidingWindow(s1,s2,9)
    dff=GetDiff(set1,set2)
    TSDexist=isTSD(dff,"DTM")
    if(TSDexist==True):
        withTSD.append(rec.id)
    else:
        noTSD.append(rec.id)
    FinalTSDlist.append(noTSD)
    FinalTSDlist.append(withTSD)
    return FinalTSDlist


def TIRpercent(seq1,seq2):
    d=compare(seq1,seq2)
    p=(40-d)/40
    p=p*100
    p=round(p, 2)
    return p


def TSDpercent(seq1,seq2):
    d=compare(seq1,seq2)
    p=(9-d)/9
    p=p*100
    p=round(p, 2)
    return p

def getTSD(tsd_dffset,fam,set1,set2):
    for i in tsd_dffset:
        if tsd_dffset[i]<TSD[fam]*0.2:
            seq1=set1[int(i.split(":")[0])]
            seq2=set2[int(i.split(":")[1])]
            return seq1,seq2


# Only Write TIR-Positive sequences
def writeTofa(file, output, List_seqwithTIR):
    used = []
    w = open(output, "a+")
    record = list(SeqIO.parse(file, "fasta"))
    for rec in record:
        if str(rec.id) in List_seqwithTIR and str(rec.id) not in used:
            s = str(rec.seq)
            l = 40
            s1 = s[0:l]
            s2_ = s[-l:]
            s2 = Seq(s2_).reverse_complement()
            s2 = str(s2)
            p_tir = TIRpercent(s1, s2)
            w.write(">" + str(rec.id) + "_TIR:" + str(s1) + "_" + str(s2_) + "_" + str(p_tir) + "_" + "\n" + str(
                rec.seq) + "\n")
            used.append(str(rec.id))


# Taget directory
# Please change whis to your folder
os.chdir("/Users/Taget/directory")

if __name__ == '__main__':
    records = list(SeqIO.parse("input_fileName", "fasta"))
    pool1 = multiprocessing.Pool(16)
    L_tir = pool1.map(CheckTIR2, records)
    pool1.close()
    pool1.join()
    noTIR = [i[0] for i in L_tir]
    withTIR = [i[1] for i in L_tir]
    noTIR = [i[0] for i in noTIR if len(i) != 0]
    withTIR = [i[0] for i in withTIR if len(i) != 0]

    print("Number of sequence without TIR: " + str(len(noTIR)))
    print("Number of sequence with TIR: " + str(len(withTIR)))

# Please change these to your input and output file names
    writeTofa("input_fileName", "output_fileName", withTIR)

