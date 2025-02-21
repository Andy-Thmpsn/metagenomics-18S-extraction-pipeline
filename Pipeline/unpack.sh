#!/bin/bash

#SBATCH --time=2:00:00 
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=20GB 
#BATCH -J taxaTarget-TOOL_17-07 -o /fslhome/orange77/NEON_processed/taxaTarget/TOOL_17-07/ --mail-user=andy_thompson@byu.edu --mail-type=BEGIN --mail-type=END --mail-type=FAIL
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

#################
## SCRIPT INFO ##
#################

# Manages options for unpacking compressed raw sequence files downloaded from external databses (e.g., NEON)
# Called by ~/NEON_processed/Pipeline/Main.
# Calls: Nothing.
# Last Update: 3/3/2023 BY AR Thompson

# max completion time: > 60 min
# nodes 1
# ntasks 2
# mem: >12 GB

#################

## VARIABLES ##

HOME_PATH=~/NEON_processed
RAW_PATH=~/compute/NEON/Raw
HMMR_PATH=$HOME_PATH/Xtract/SSU_BLAST                                              ## BEWARE: Path name and variable_name match what is in ~/NEON_processed/Pipeline/pipeline.sh. DO NOT CHANGE.
MTDA_PATH=$HMMR_PATH/mtgnmDnaXtrcts
QC_PATH=$HOME_PATH/QCtrl

#SiteList=AllSiteNames1.list
#SiteList=AllSiteNames.list

## ROUTINES ##


manage_subds () {

 find $RAW_PATH/$site_name/ -mindepth 2 -type f -print -exec mv {} $RAW_PATH/$site_name/. \; ## find all files in all subdirectories of a given directory and move them into parent directory (NOTE: does not remove the subdirectories)
 #find $RAW_PATH/$site_name/ -mindepth 1 -type d -print -delete                               ## this one should remove the directories (NOTE: the '-delete' option only removes empty directories.
 
}


unpack_tar () {

 ## 1. UNPACK MULTIPLE *.tar.gz FILES AT ONCE              
  # echo "Opening the Tar Ball..."
 # while read sample_name;
 #  do
 #   tar zxvf - -iC $RAW_PATH/$site_name/$sample_name/
 #  done<$MTDA_PATH/${site_name}.list
 echo "$site_name"   
  cat $RAW_PATH/$site_name/*.tar.gz | tar zxvf - -iC $RAW_PATH/$site_name/
  
  manage_subds $site_name

}


unpack_gz () {
 
 ## 1. UNPACK MULTIPLE *.gz FILES AT ONCE

  echo "Unpacking *.gz files: check 2"
 
  for f in $RAW_PATH/$site_name/*fastq.gz; do
   STEM=$(basename "${f}" .gz)
   gunzip -c "${f}" > $RAW_PATH/$site_name/"${STEM}"
  done
 
  manage_subds $site_name
} 


unpack_zip () {

 ## 1. UNPACK MULTIPLE *.zip FILES AT ONCE
  #echo "Site Name (inside fxn): $site_name"
  #cat $RAW_PATH/${site_name}/*.zip | unzip
  for f in $RAW_PATH/$site_name/*.zip; do
   unzip -o $f -d $RAW_PATH/$site_name/.;
  done
 
  manage_subds $site_name

}

usage (){

cat << EOF

Usage:
   -t  tar ~ UNPACK MULTIPLE *.tar.gz FILES AT ONCE

   -g  gz ~ UNPACK MULTIPLE *.gz FILES AT ONCE

   -z  zip ~ UNPACK MULTIPLE *.zip FILES AT ONCE
            
   -h Usage

EOF

}



while getopts 't:g:z:h' OPTION; do
    case "$OPTION" in
        t) 
           echo "unpacking tar files for ..."
           site_name="$OPTARG"
           unpack_tar "$site_name"
           ;;
        g) 
           site_name="$OPTARG"
           unpack_gz "$site_name"
           ;;
        z)
           site_name="$OPTARG"
           unpack_zip "$site_name"
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


