#!/bin/bash

#SBATCH --time=8:00:00 
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --mem=4GB 
#BATCH -J taxaTarget-TOOL_17-07 -o /fslhome/orange77/NEON_processed/taxaTarget/TOOL_17-07/ --mail-user=andy_thompson@byu.edu --mail-type=BEGIN --mail-type=END --mail-type=FAIL
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

#################
## SCRIPT INFO ##
#################

# Assesses some info about the OTU table run created by the SSU_EXTRACT pipeline, namely... 
    #I think I originally built this to produce information I didn't realize R could provide more easily. Still, I think there could be space for a stats.sh in the pipeline.
# Called by none
# Calls: none
# Last Update: 3/3/2023 BY AR Thompson

#################

## PENDING UPDATES ##

  ## Ready for expansion to other stat assessments (e.g.?)


## VARIABLES ##

HOME_PATH=~/NEON_processed
RAW_PATH=~/compute/NEON/Raw
LINKS_PATH=$HOME_PATH/DataLinks
HMMR_PATH=$HOME_PATH/Xtract/SSU_BLAST      ## BEWARE: Path name and variable_name match what is in ~/NEON_processed/Pipeline/pipeline.sh. DO NOT CHANGE.
MTDA_PATH=$HMMR_PATH/mtgnmDnaXtrcts
QC_PATH=$HOME_PATH/QCtrl
TRMD_PATH=~/compute/NEON/Trimmed

scripts=$HOME_PATH/Pipeline

SiteList=AllSiteNames.list                               
OTUtbl_long=AllSites_long.otu                            
AllTaxa_list=AllSites_list.taxa                          
AllTaxIDTable=AllSites.BH_txID.tbl                       

statsSITES=stats_reads.site
statsSAMPLES=stats_reads.sample

stats (){


echo "" >$HMMR_PATH/$statsSITES

## header
print "Note that sites with '0' reads were those that the pipeline failed to process.\n\n\n"


## total reads per site 

while read site; do

 printf "$site: " >> $HMMR_PATH/$statsSITES
 
 grep -c "$site" $HMMR_PATH/$AllTaxIDTable >> $HMMR_PATH/$statsSITES
 
done<$HMMR_PATH/$SiteList



## total reads per sample

while read sample; do

 printf "$sample: " >> $HMMR_PATH/$statsSAMPLES
 
 grep -c "$site" $HMMR_PATH/$AllTaxIDTable >> $HMMR_PATH/$statsSAMPLES
 
done<$HMMR_PATH/$SiteList



## taxa per site

while read site; do

 printf "$site: " >> $HMMR_PATH/$statsSITES
 
 grep -c "$site" $HMMR_PATH/$AllTaxIDTable >> $HMMR_PATH/$statsSITES
 
done<$HMMR_PATH/$SiteList


}



usage (){

cat << EOF

Usage:

    -f Get stats

    -h Usage

EOF

}


while getopts 'sh' OPTION; do
    case "$OPTION" in
        s) 
           echo "Getting Stats"
           stats
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


