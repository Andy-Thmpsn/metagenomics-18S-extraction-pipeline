#!/bin/bash

#SBATCH --job-name=neon_pipeline
#SBATCH --output=logs/pipelineNEON_%j.out
#SBATCH --error=logs/pipelineNEON_%j.err
#SBATCH --time=16:00:00 
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=24GB 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=athompson@cmc.edu

#export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

#################
## SCRIPT INFO ##
#################

# Runs the core steps of the SSU_EXTRACT pipeline, namely  sequence extraction and identification (nhmmer and blastn)
    # Optionally runs Metaxa2, TaxaTarget, and EukDetect
# Called by ~/NEON_processed/Pipeline/Main.sh
# Calls: filter.py
# Last Update: 3/3/2023 BY AR Thompson

# Part 2 & 3 - rDNA extraction and Identification
 # mem > 24GB for a few, <20GB for most
 # <3 hours for most, <6 hours max
 # ntasks < 2

#################

## PENDING UPDATES ##

 ## Optimize hmmr and blastn search parameters


## VARIABLES ##
NAME="BMI_Plate42WellH10"
SITE="ABBY_2017-06"


SOFTWARE=/hopper/groups/yoosephlab/software
HOME_PATH=/hopper/groups/yoosephlab/athompson/NEON_18S_Xtrct
RAW_PATH=$HOME_PATH/NEON_raw
PRCS_PATH=$HOME_PATH/NEON_processed
QC_PATH=$PRCS_PATH/QCtrl
MRGD_PATH=$HOME_PATH/NEON_merged
TRMD_PATH=$HOME_PATH/NEON_trimmed
MTXA_PATH=$PRCS_PATH/Xtract/Metaxa2
HMMR_PATH=$PRCS_PATH/Xtract/SSU_BLAST

SCRPT_PATH=$HOME_PATH/Pipeline

THREADS=4

printf "\n\n*********************************************************\n"
printf "\nBEGIN rDNA EXTRACTION.\n"

printf "Processing: $NAME from $SITE\n\n" 



######################################################
##### PART 1  REFORMAT FILES FROM FASTQ TO FASTA #####
######################################################

printf "Part 1 - Reformat Files From FASTQ to FASTA.\n"

## 1. FLASH Merged fastq files

SFX_FQ=.trimmed.fastq
SFX_FQ_2=.trimmed.fastq2
SFX_FA=.trmd.fasta

## 2. Make copies of raw fastq files so we can convert one copy to fasta
 cp $TRMD_PATH/$SITE/"$NAME""$SFX_FQ" $TRMD_PATH/$SITE/"$NAME""$SFX_FQ_2" 
 #cp $TRMD_PATH/$SITE/"$NAME""$SFX_FQ_2" $TRMD_PATH/$SITE/"$NAME""$SFX_FQ" 


## 3. Convert each read file from fastq to fasta format
 paste - - - - < $TRMD_PATH/$SITE/"$NAME""$SFX_FQ_2" | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $TRMD_PATH/$SITE/"$NAME""$SFX_FA"

 printf "\tFastq Converted to Fasta Format!\n"

## 4. STORE FILES THAT WON'T BE IN USE ANYMORE IN TAR.GZ BALL

 printf "\tTidying up the data ...\n\t"

 tar -cvzf $TRMD_PATH/$SITE/"$NAME"_trm.tar.gz $TRMD_PATH/$SITE/"$NAME".trimmed.fastq2  --remove-files

 printf "\tFastq Files Zipped Away for Another Day!\n"


 printf "PART 1: COMPLETE!\n\n"



######################################################
##### PART 2 rDNA EXTRACTION FROM FASTA FILES    #####
######################################################

printf "\nPART 2: BEGIN rDNA EXTRACTION FROM FASTA FILES.\n"

## A. METAXA2.2.3 AGAINST DEFAULT SILVA DATABASE
# mkdir -p $MTXA_PATH/$SITE
# mkdir -p $MTXA_PATH/$SITE/$NAME
# module load perl/5.28
# metaxa2 -i $TRMD_PATH/$SITE/"$NAME"_R1.paired.fasta -o $MTXA_PATH/$SITE/$NAME/$NAME --cpu 4 --plus T --table T
#metaxa2 -i $TRMD_PATH/$SITE/"$NAME"_R1.paired.fasta -o $MTXA_PATH/$SITE/$NAME/$NAME -f fasta -E 1e-5 -N 1 --cpu 4 --plus T --table T --allow_single_domain 1e-5,0


## B1. HMMER WITH EUK_BARR.HMM PROFILE AGAINST PRE-PROCESSED SITE FILES
 
 SFX_hmmTBL=.hmm.tbl

 mkdir -p $HMMR_PATH/$SITE/

 echo $TRMD_PATH/$SITE/"$NAME""$SFX_FA"

 nhmmer --tblout $HMMR_PATH/$SITE/"$NAME""$SFX_hmmTBL" -o $HMMR_PATH/$SITE/"$NAME"_hmm.out --notextw --cpu $THREADS ~/programs/hmmer-3.2.1/euk_barr.hmm $TRMD_PATH/$SITE/"$NAME""$SFX_FA"

 printf "\tnhmmer rRNA extraction completed\n"
 

## C. EXTRACT CORRESPONDING BEST HIT 18S SEQUENCES FROM ORIGINAL FASTA FILE

 printf "\tExtracting Best Hits from Original Fasta File...\n"

SFX_hmmBH=_hmm.besthit.tbl
SFX_hmmQUERY=.hmm.queries
SFX_hmmFASTA=.hmm.fasta
MATCH=NB551 SeqID=16

# Grabs 18S E-value < 1e-5; highest bit score breaks ties
 awk '$3 ~/18S_rRNA/ && $13 ~/e/' $HMMR_PATH/$SITE/"$NAME""$SFX_hmmTBL" | sort -k1,1 -k14,14nr | sort -u -k1,1 > $HMMR_PATH/$SITE/"$NAME""$SFX_hmmBH"
 awk '{print ">"$1}' $HMMR_PATH/$SITE/"$NAME""$SFX_hmmBH" > $HMMR_PATH/$SITE/"$NAME""$SFX_hmmQUERY"

# Converts sequence ID in BestHit.tbl file ($1) to fasta format (>)
 grep -F -A1 --no-group-separator -f $HMMR_PATH/$SITE/"$NAME""$SFX_hmmQUERY" $TRMD_PATH/$SITE/"$NAME""$SFX_FA" > $HMMR_PATH/$SITE/"$NAME""$SFX_hmmFASTA"

 printf "\tBest Hits Extracted from Fasta File!\n"
 printf "PART 2: rDNA EXTRACTION COMPLETE!\n\n"



############################################
##### PART 3 - SEQUENCE IDENTIFICATION #####
############################################

 printf "\nPART 3: BEGIN SEQUENCE IDENTIFICATION.\n"

SFX_blstOUT=.hmmblst.tbl6
SFX_blstOUT_filtered=.hmmblst.tbl6.filt
SFX_blstBH=.hmmblst.BH.tbl
SFX_blstBHtax=.hmmblst.BH_txID.tbl
Seqs_lost_table=SeqsFailingCutoffs.tbl
PFX_ALL=All_
DATABASE=$HOME_PATH/database/pr2/pr2db_blstfmt
TAXID_File=$HOME_PATH/database/pr2/pr2_version_4.14.0_SSU_mothur.tax
SeqID_Cutoff=93.000
SeqLength_Cutoff=125 #in bp
QueryCoverage_Cutoff=0.90 #percent

## 1. BLAST EXTRACTED SEQUENCES AGAINST PR2 DATABASE
 blastn -num_threads $THREADS -db $DATABASE -query $HMMR_PATH/$SITE/"$NAME""$SFX_hmmFASTA" -outfmt 6 -out $HMMR_PATH/$SITE/"$NAME""$SFX_blstOUT"  # Sequences from euk_barr + nhmmr (custom)

####////!!!! THIS SECTION NEEDS TESTING/DOUBLE CHECKING !!!!!\\\\####

## 2. FILTER TBL6 BLAST OUTPUT FOR LOW %ID SEQUENCE HITS, LOW QUERY COVERAGE, and SHORT SEQUENCE LENGTH 2/10/2023; revised 8/27/2025
  #Argument input sequence matters: it should be "file_name", "SeqID", then "SeqLength"
  python $SCRPT_PATH/filter.py $HMMR_PATH/$SITE/"$NAME""$SFX_blstOUT" $SeqID_Cutoff $SeqLength_Cutoff $QueryCoverage_Cutoff > $HMMR_PATH/$SITE/"$NAME""$SFX_blstOUT_filtered" 
  echo "Min_ID: $SeqID_Cutoff"
  echo "Min_Length: $SeqLength_Cutoff "
  echo "Min_Query_Coverage: $QueryCoverage_Cutoff percent"


## 3. CALCULATE SEQUENCES LOST FROM FILTERING STEP AND ADD TO SEQS_LOST TABLE
  #NEEDS TESTING 8/22/2025
  #Seqs_lost = unique length of input file ($SFX_blstOUT) minus unique length of output file ($SFX_blstOUT_filtered)  
  SeqsHit=$(wc -l < "$HMMR_PATH"/"$SITE"/"$NAME""$SFX_blstOUT")
  SeqsAfterFiltering=$(wc -l < $HMMR_PATH/$SITE/"$NAME""$SFX_blstOUT_filtered")
  Seqs_lost=$(( $(wc -l < $HMMR_PATH/$SITE/"$NAME""$SFX_blstOUT") - $(wc -l < $HMMR_PATH/$SITE/"$NAME""$SFX_blstOUT_filtered") ))
  #Add $Seqs_lost to SFX_SeqsLost.tbl
  echo "$Site\t$NAME\t$SeqsHit\t$SeqsFiltered\t$SeqsLost"  >> $HMMR_PATH/$SITE/"$NAME""$Seqs_lost_table"
  echo "Seqs_Lost: $Seqs_lost"


## 4. CONVERT HITS FROM BLASTN TO BEST HIT (ONE HIT ONLY)
sort -k1,1 -k12,12nr -k11,11n $HMMR_PATH/$SITE/"$NAME""$SFX_blstOUT_filtered" | sort -u -k1,1 >  $HMMR_PATH/$SITE/"$NAME""$SFX_blstBH"                                
echo "Hits Converted..."

## 5. MERGE BEST HIT TABLE WITH ASSOCIATED TAXONOMIC LINEAGE
join --nocheck-order -1 1 -2 2 <(sort -k 1 $TAXID_File) <(sort -k 2 $HMMR_PATH/$SITE/"$NAME""$SFX_blstBH") | tr ' ' "\t" > $HMMR_PATH/$SITE/"$NAME""$SFX_blstBHtax"

echo "Taxonomy added ..."

## 6. ADD $NAME'
sed -i 's/;/\t/g' $HMMR_PATH/$SITE/"$NAME""$SFX_blstBHtax" ;  sed -i "s/$/\t$NAME/g" $HMMR_PATH/$SITE/"$NAME""$SFX_blstBHtax"

echo "Sample Names added" 
printf "\nPART 3: SEQUENCE IDENTIFICATION COMPLETE!\n\n"



#######################################################
##### PART 4 - EXTERNAL GENE DETECTION SOFTWARES  #####
######################################################

## EUKDETECT

# conda init bash
# conda activate eukdetect
# eukdetect --mode runall --configfile ~/programs/EukDetect/NEON_configfile.yml --cores 4

## TAXATARGET

# module load python/3.8
# python ~/programs/taxaTarget/run_pipeline_scripts/run_protist_pipeline_fda.py -r $RAW_PATH/$SITE/"$NAME"_R1.fastq -r $RAW_PATH/$SITE/"$NAME"_R2.fastq -e ~/programs/taxaTarget/run_pipeline_scripts/environment.txt -o $OUT_PATH/taxaTarget/$SITE/$NAME -t 12 





format_time() {
 ((h=${1}/3600))
 ((m=(${1}%3600)/60))
 ((s=${1}%60))
 printf "%02d:%02d:%02d\n" $h $m $s
}

printf "Script completed in $(format_time $SECONDS)\n"
printf "Processing for "$NAME" complete.\n\n"
printf "\n***********************************************\n"
