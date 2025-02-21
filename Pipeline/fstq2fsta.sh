#!/bin/bash

#SBATCH --time=1:00:00 
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4GB 
#BATCH -J fstq2fsta_testing -o /fslhome/orange77/NEON_processed/LogFiles/fstq2fsta/ --mail-user=orange7724@gmail.com --mail-type=BEGIN --mail-type=END --mail-type=FAIL

#################
## SCRIPT INFO ##
#################

# Converts fastq files to fasta formatted files
    #This same code (but updated) is in pipeline.sh. It isn't complicated, but it could be something that could be separated later on? I'll keep it for now.
# Called by ~/NEON_processed/Pipeline/Main.sh
# Calls: none
# Last Update: 3/3/2023 BY AR Thompson

#################

## PENDING UPDATES ##

# None

## GENERAL VARIABLES ##

SITE=BART_2016-07
NAME=BMI_Plate1WellH7_mms
PRE_PATH=~/compute/NEON
RAW_PATH=$PRE_PATH/Raw
PRCS_PATH=~/NEON_processed/Test_fasta

##  REFORMATTING VARIABLES ##

SFX_rawR1_1=_R1.fastq         # original fastq suffix
SFX_rawR2_1=_R2.fastq
SFX_rawR1_2=_R1_2.fastq        # backup copy suffix
SFX_rawR2_2=_R2_2.fastq
SFX_rawR1_fa=_R1.fasta         # fasta file suffix
SFX_rawR2_fa=_R2.fasta
 
## R1 & R2 raw fastq ##

mkdir -p $PRCS_PATH/$SITE/

# Make copies of raw fastq files so we can convert one copy to fasta
cp $RAW_PATH/$SITE/"$NAME""$SFX_rawR1_1" $PRCS_PATH/$SITE/"$NAME""$SFX_rawR1_2"
cp $RAW_PATH/$SITE/"$NAME""$SFX_rawR2_1" $PRCS_PATH/$SITE/"$NAME""$SFX_rawR2_2"

# Convert each read file from fastq to fasta format
paste - - - - < $PRCS_PATH/$SITE/"$NAME""$SFX_rawR1_2" | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $PRCS_PATH/$SITE/"$NAME""$SFX_rawR1_fa" 
paste - - - - < $PRCS_PATH/$SITE/"$NAME""$SFX_rawR2_2" | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $PRCS_PATH/$SITE/"$NAME""$SFX_rawR2_fa"

# Remove copies of raw fastq files to preserve space in original folder
rm -i $PRCS_PATH/$SITE/"$NAME""$SFX_rawR1_2" 
rm -i $PRCS_PATH/$SITE/"$NAME""$SFX_rawR2_2" 

