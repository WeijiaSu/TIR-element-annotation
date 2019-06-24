#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
#
#SBATCH --time=30:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=5   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=weijia@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#SBATCH --job-name='TIR-Learner_Maize'

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

############## Load dependencies ##############
#module load python/3            ##############
#module load cd-hit              ##############
#module load gcc                 ##############
###############################################


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
############################################### Change to your parameters #############################################################
                                                                                                                     ################## 
genomeFile="/work/LAS/.../Genomes/MaizeB73/MaizeB73_chr1.fa" # Path to your genome file                              ##################
genomeName="MaizeB73_Chr1"  # Genome Name                                                                            ##################
path="/work/LAS/.../research/TIR-Learner" # Path to the TIR-Learner source code                                      ##################
dir="/work/LAS/.../research/13_Maize" # your target directory                                                        ##################
t=80 # Number of processor                                                                                           ##################      
grfp="/work/LAS/thomasp-lab/weijia/research/software/GenericRepeatFinder/bin" # Path to the GRF program              ##################
species="Maize" # Species one of the following : Maize , Rice or Others                                              ##################  
l=5000 # Max length of TIR                                                                                           ##################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################



echo "############################################################ Pre-Processing ###########################################################"


python3 $path/pre.py -g $genomeFile -name $genomeName
for i in Module1 Module2 Module3; do mkdir $i; done


######################################################################################################################################################
############################################################################Module   1 ###############################################################
######################################################################################################################################################
echo "############################################################ Module   1 Begin ###########################################################"

cd $dir"/Module1"
mkdir $genomeName
mkdir temp
echo "Module 1, Step 1: Blast Genome against Reference Library"

python3 $path/Module1/Blast_Ref.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1" -s $species
cp $genomeName/*blast* temp/

echo "Module 1, Step 2: Select 100% coverage entries from Blast results"
python3 $path/Module1/Fullcov.py  -name $genomeName -p $path -d $dir"/Module1" -s $species

echo "Module1, Step 3: Making blastDB and get candidate sequences"
#
python3 $path/Module1/GetSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module1, Step 4: Check TIR and TSD"

python3 $path/Module1/CheckTIRTSD.py -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module1, Step 5: Write to Gff3"

python3 $path/Module1/WriteToGff_M1.py -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "Module1, Step 6: Check Low Complexity"

python3 $path/Module1/Lowcomp_M1.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module1"

echo "############################################################ Module   1 Finished ###########################################################"
echo "############################################################ Module   2 Begin ###########################################################"


cd $dir"/Module2/"
mkdir $genomeName
mkdir temp

echo "Module 2, Step 1: Split Genome and  Run GRF program to find Inverted Repeats"
python3 $path/Module2/RunGRF.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2" -grfp $grfp -l $l
cp -r $genomeName/$genomeName* temp/

echo "Module 2, Step 2: Process GRF results"
python3 $path/Module2/ProcessGRFmite.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2 , Step 3 : GRF result blast reference sequences"
python3 $path/Module2/Blast_ref.py -name $genomeName -p $path -t $t -d $dir"/Module2" -s $species
cp $genomeName/*80 temp/



echo "Module 2 , Step 4: Get sequences from 80% similarity"
python3 $path/Module2/GetSeq_M2.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2 , Step 5: Check TIR and TSD"
python3 $path/Module2/CheckTIRTSD_M2.py -name $genomeName -p $path -t $t -d $dir"/Module2"

echo "Module 2 , Step 6: Write to gff3"
python3 $path/Module2/WriteToGff_M2.py -name $genomeName -p $path -t $t -d $dir"/Module2"
#
echo "Module 2 , Step 7: Remove Low complexity"
python3 $path/Module2/Lowcomp_M2.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module2"
#
echo "############################################################ Module   2 Finished ###########################################################"
echo "############################################################ Module   3 Begin ###########################################################"



cd $dir"/Module3/"

mkdir $genomeName
mv ../Module2/$genomeName/*_nonHomo.fa $genomeName/
mkdir temp


echo "Module 3, Step 1: Get dataset"
python3 $path/Module3/getDataset.py -name $genomeName -p $path -t $t -d $dir"/Module3"
cp $genomeName/*.csv temp/
cp $genomeName/*nonHomo.fa temp/
#
echo "Module 3, Step 2: ML prediction"
python3 $path/Module3/ML_Ensemble.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3" -s $species

echo "Module 3, Step 3: Check TIR/TSD"
python3 $path/Module3/CheckTIRTSD_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3"

echo "Module 3, Step 4: Write to Gff"
python3 $path/Module3/WriteToGff_M3.py -name $genomeName -p $path -t $t -d $dir"/Module3"

echo "Module 3, Step 5: Remove Low complexity"
python3 $path/Module3/Lowcomp_M3.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir"/Module3"

echo "############################################################ Module   3 Finished ###########################################################"
#
echo "############################################################ Post Processing  ###########################################################"
#
cd $dir
#
mkdir $genomeName
#
#
echo "Get Final GFF" 
python3 $path/CombineAll.py -name $genomeName -p $path -t $t -d $dir

mv *.gff3 $genomeName
rm *Low

echo "Get fasta file"
python3 $path/GetAllSeq.py -g $genomeFile -name $genomeName -p $path -t $t -d $dir

mkdir TIR-Learner-Result
mv $genomeName/*FinalAnn.gff3 TIR-Learner-Result/
mv $genomeName/*FinalAnn.fa TIR-Learner-Result/
rm -r $genomeName"_combine"
rm -r $genomeName


echo "############################################################ TIR-Learner is finished! ###########################################################"
#
