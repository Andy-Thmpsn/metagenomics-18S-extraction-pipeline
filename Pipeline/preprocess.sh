#!/bin/bash

#SBATCH --time=16:00:00 
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=24GB 
#BATCH -J taxaTarget-TOOL_17-07 -o /fslhome/orange77/NEON_processed/taxaTarget/TOOL_17-07/ --mail-user=andy_thompson@byu.edu --mail-type=BEGIN --mail-type=END --mail-type=FAIL
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE


#################
## SCRIPT INFO ##
#################

# Manages the preprocessing options (sort sequences, fastqc, merge, and trim) for the SSU_EXTRACT pipeline.
# Called by ~/NEON_processed/Pipeline/Main.sh
# Calls: none
# Last Update: 3/3/2023 BY AR Thompson

# Part 1 - preFQC, trm, pstFQC
 # mem > 12GB for some 
 # 16 > hours for some (or maybe its because they run out of memory, but don't exit?)
 # ntasks 1

#################


## VARIABLES ##

echo "Site: $SITE"
echo "Name: $NAME"
PRE_PATH=~/compute/NEON
RAW_PATH=$PRE_PATH/Raw
PRCS_PATH=~/NEON_processed
QC_PATH=$PRCS_PATH/QCtrl
TRMD_PATH=$PRE_PATH/Trimmed
MRGD_PATH=$PRE_PATH/Merged

SCRPT_PATH=$PRCS_PATH/Pipeline

THREADS=4


## FUNCTIONS ##


fix_order(){

 echo "Sorting sample files"

 ~/programs/fastq-tools-master/src/fastq-sort $RAW_PATH/$SITE/"$NAME"_R1.fastq > $RAW_PATH/$SITE/"$NAME"_R1_sorted.fastq
 ~/programs/fastq-tools-master/src/fastq-sort $RAW_PATH/$SITE/"$NAME"_R2.fastq > $RAW_PATH/$SITE/"$NAME"_R2_sorted.fastq
 
 
}

clean_up(){

 echo "Removing Sorted Fastq files from Raw Sequence File Directory"

 rm $RAW_PATH/$SITE/"$NAME"_R1_sorted.fastq
 rm $RAW_PATH/$SITE/"$NAME"_R2_sorted.fastq

}

fastqc(){
 
 if [ "$1" == "pre" ]; then
 
  echo "Running Fastqc on sequence files BEFORE trimming..."

   mkdir -p $QC_PATH/Pre/$SITE/
 
#   find $MRGD_PATH/$SITE/. -name '*.extendedFrags.fastq' -exec ~/programs/FastQC/FastQC/fastqc -o $QC_PATH/Pre/$SITE {} \;

   ~/programs/FastQC/FastQC/fastqc $MRGD_PATH/$SITE/"$NAME".extendedFrags.fastq -o $QC_PATH/Pre/$SITE

 echo "Fastqc complete!"
   
 

 elif [ "$1" == "post" ]; then

  echo "Running Fastqc on sequence files AFTER trimming..."
 
   mkdir -p $QC_PATH/Post/$SITE/
 
#   find $TRMD_PATH/$SITE/. -name '*.paired.fastq' -exec ~/programs/FastQC/FastQC/fastqc -o $QC_PATH/Post/$SITE {} \;

   ~/programs/FastQC/FastQC/fastqc $TRMD_PATH/$SITE/"$NAME".trimmed.fastq -o $QC_PATH/Post/$SITE

  echo "Fastqc complete!"  
 
 

 else echo "No matching input parameters found (Check capitalizations and spellings.)"; usage

 fi  

}

merge(){

 mkdir -p $MRGD_PATH/$SITE/

 printf $MRGD_PATH/$SITE/"$NAME"

 echo "\n"

 #cat $TRMD_PATH/$SITE/"$NAME"_R1.paired.fastq $TRMD_PATH/$SITE/"$NAME"_R2.paired.fastq > $TRMD_PATH/$SITE/"$NAME"_R12.paired.fastq

 if [ "$1" == "pre" ]; then

#  fix_order

  echo "Merging R1 and R2 files with FLASH"

#  flash -M 150 -t $THREADS -q -o "$NAME" -d $MRGD_PATH/$SITE/ $RAW_PATH/$SITE/"$NAME"_R1_sorted.fastq $RAW_PATH/$SITE/"$NAME"_R2_sorted.fastq 

  clean_up

 elif [ "$1" == "post" ]; then

  fix_order
    
  flash -M 150 -t $THREADS -q -o "$NAME" -d $MRGD_PATH/$SITE/ $TRMD_PATH/$SITE/"$NAME"_R1.paired.fastq $TRMD_PATH/$SITE/"$NAME"_R2.paired.fastq

  clean_up

 fi

}


trimming(){

 mkdir -p $TRMD_PATH/$SITE/

 java -jar ~/programs/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $THREADS -phred33 $MRGD_PATH/$SITE/"$NAME".extendedFrags.fastq $TRMD_PATH/$SITE/"$NAME".trimmed.fastq LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:50
 
}

usage (){

cat << EOF

Usage:

    -f FastQC (requires an argument)
        pre ~ sends results to the pre-trimming folders
        post ~ sends results to the post-trimming folders   

    -m Merge Paired-end Reads (requires an argument)

    -t Trim reads (parameters are currently hardcoded)
    
    -s fix_order (sort) sequences in fastq files

    -c Clean Up (remove sorted fastq files)

    -a Run all Pre-processing steps in sequence (NEEDS UPDATING; Not Functional)
        Sort Raw fastq files (fix_order)
        Merge Reads
            Defaults to "pre"
        FastQC (pre)
        Trim Reads
        FastQC (post)
        Clean up (remove sorted fastq files)
    
    -h Usage

EOF

}

## USER INTERACTION ##


while getopts 'f:m:tscah' OPTION; do
    case ${OPTION} in
        f) 
           echo "Running FastQC"
           fastqc "$OPTARG"
           ;;
        m) 
           echo "Merging Paired Reads..."
           merge "$OPTARG"
           ;;
        t)
           echo "Trimming Reads..."
           trimming
           ;;
 
        s) 
           echo "Sorting Read files"
           fix_order
           ;;
        c) 
           echo "Consolidating Raw Sequence Files" 
           clean_up
           ;;
        a)
           echo "Running all pre-processing steps in sequence"
           echo "Sorting Read files"
          # fix_order
           echo "Merging Paired Reads..."
          # merge "pre"
           echo "Running FastQC (pre trimming)"
          # fastqc "pre"
           echo "Trimming Reads ..." 
          # trimming 
           echo "Running FastQC (post trimming)"
          # fastqc "post"
           echo "Consolidating Raw Sequence Files" 
          # clean_up
           
           ;;
        h) 
           usage 
           ;;
        ?)
           usage
           #echo "Script usage: $(basename \$0) [-l] [-h] [-a somevalue]" >&2
           exit 1
           ;;
    esac
done
shift "$((OPTIND -1))"

