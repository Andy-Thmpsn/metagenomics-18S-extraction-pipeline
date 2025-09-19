#!/bin/bash

#SBATCH --job-name=NEON_otutable
#SBATCH --output=logs/otutableNEON_%j.out
#SBATCH --error=logs/otutableNEON_%j.err
#SBATCH --time=24:00:00 
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=8GB 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=athompson@cmc.edu


export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE


#################
## SCRIPT INFO ##
#################

# Manages the postprocessing options (build super BH table, build OTU tables [tall and long], and add metadata to the OTU table ) for the SSU_EXTRACT pipeline.
# Called by ~/NEON_processed/Pipeline/Main.sh
# Calls: none
# Last Update: 3/3/2023 BY AR Thompson

# Resource Use: Building the OTU table 
 # mem < 10GB for if partition is 20 samples each 
 # 4-8 hours per partition when building the OTU table (long)
 # ntasks 1

#################

## PENDING UPDATES ##

 ## try different databases (SILVA/NCBI) and compare to PR2
 ## revisit my "Best Hit" algorithm. 

## VARIABLES ##

BASE_PATH=/hopper/groups/yoosephlab/athompson/NEON_18S_Xtrct
HMMR_PATH=${BASE_PATH}/NEON_processed/Xtract/SSU_BLAST
PARTITION_PATH=$HMMR_PATH/Partitions
SFX_blstBHtax=.hmmblst.BH_txID.tbl   ## BEWARE: Suffix and variable_name match what is in ~/NEON_processed/Main.sh. DO NOT CHANGE WITHOUT CHANGING THAT FILE TOO.
PFX_ALL=All_                         ## BEWARE: Prefix and variable_name match what is in ~/NEON_processed/Pipeline/pipeline.sh. DO NOT CHANGE WITHOUT CHANGING THAT FILE TOO.
BestHits=$PFX_All"*"$SFX_blstBHtax
AllTaxa_head=AllSites_header.taxa
OTUtbl_tall=AllSites_tall.otu
AllTaxonomyTable=AllSites.taxonomy
#METADATA_tbl=string
MTDA_PATH=$HMMR_PATH/mtgnmDnaXtrcts 
SiteList=AllSiteNames.list
AllTaxIDTable=AllSites.BH_txID.tbl
SampleList=AllSites.samples
AllTaxa_list=AllSites_list.taxa
OTUtbl_long_pre=AllSites_long_pre.otu
OTUtbl_long_final=AllSites_long_final.otu

echo "Sample List: $SampleList"
echo "OTU Table: $OTUtbl_long_pre"
echo "Taxa List: $AllTaxa_list"
echo "Tax ID Table: $AllTaxIDTable"


## FUNCTIONS ##


gather_tables () {

    printf "\tGathering Files to Make OTU table...\n"

 ## Gather BestHit_files from all sites
    find $HMMR_PATH/. -mindepth 2 -type f -name All_"*"$SFX_blstBHtax -exec cat {} \; > $HMMR_PATH/$AllTaxIDTable                                            
 
 ## Output site list as a file (when making the long OTU table)
    awk -v OFS='\t' '{print $(NF-1), $NF}' $HMMR_PATH/$AllTaxIDTable | sort -u > $HMMR_PATH/$SampleList                                        

    awk '{print $8}' $HMMR_PATH/$AllTaxIDTable | sort -u > $HMMR_PATH/$AllTaxa_list # | head -n 10

 ## Get a list of all 1) Families, 2) Genera, or 3) Species in BestHit
    #awk -v OFS="\t" '{print $2, $3, $4, $5, $6, $7}' $HMMR_PATH/$AllTaxIDTable | sort -u | # tee $HMMR_PATH/$AllTaxa_list | #TO FAMILY
    awk -v OFS="\t" '{print $2, $3, $4, $5, $6, $7, $8}' $HMMR_PATH/$AllTaxIDTable | sort -u -k7,7 | #head -n 10 | tee check_list.taxa | #TO GENUS
    #awk -v OFS="\t" '{print $2, $3, $4, $5, $6, $7, $8, $9}' $HMMR_PATH/$AllTaxIDTable | sort -u | # tee $HMMR_PATH/$AllTaxa_list | #TO SPECIES                            

 
    ## The following code is for transposing the taxonomic table. It works, but I found it on the internets and don't exactly know what it's doing. I would like to figure that out at some point...
    awk '
    {       
        for (i=1; i<=NF; i++) {
            a[NR,i] = $i
       }   
    }
    NF>p { p = NF }
    END {
        for (j=1; j<=p; j++) {
            str=a[1,j]
            for(i=2; i<=NR; i++){
                str=str" "a[i,j];
            }
            print str
        }
    }' > $HMMR_PATH/$AllTaxa_head

    printf "\tAll Files Gathered!\n"
    printf "\tFormatting OTU Table...\n"
}

OTU_tall () {
    ## MAKE OTU TABLE
    # OTU as Y/tall
    #awk '{print $(NF-1)}' $HMMR_PATH/$AllTaxIDTable | sort -u | tee $HMMR_PATH/$SampleList |  paste -sd "\t" | sed '1s/^/\t/'> $HMMR_PATH/$OTUtbl_tall           ## Print header (list of sites) for OTU table - output site list as a file as well
    ##while read Taxon;
    #do
    #   printf '%s\t' $Taxon >> $HMMR_PATH/$OTUtbl_tall
    #   while read sample;
    #   do
    #    grep -c -e "$Taxon.*$sample" $HMMR_PATH/$AllTaxIDTable | tr '\n' '\t' >> $HMMR_PATH/$OTUtbl_tall                                     ## 
    #   done<$HMMR_PATH/$SampleList
    #   printf "\n" >> $HMMR_PATH/$OTUtbl_tall
    #done< $HMMR_PATH/$AllTaxa_list
echo "Making a Tall OTU table!"
}


OTU_long () {   ## OTU as X/long

    # when sample list is split into partitions of 20 each, this approach takes between 1 and 9 hours per partition to run (!!!!)
    # there must be a faster way...

awk -v OFS="\t" '
{
    taxon = $8      # column for taxon
    sample = $21     # column for sample
    site = $22      # column for site
    count[sample, site, taxon]++
}
END {
    for (key in count) {
        split(key, arr, SUBSEP)
        print arr[1], arr[2], arr[3], count[key]
    }
}' $HMMR_PATH/$AllTaxIDTable > $HMMR_PATH/$OTUtbl_long_pre

printf "\tOTU Table Complete!\n"

exit

# OLD, SUPER SLOW WAY (UPDATED9/8/2025)
#clear previous OTU tbl (long) file
truncate -s 0 $PARTITION_PATH/$OTUtbl_long_pre  
 #printf "" > $HMMR_PATH/$OTUtbl_long_pre

readarray -t samples < $PARTITION_PATH/$SampleList
readarray -t taxa < $PARTITION_PATH/$AllTaxa_list

for sample in "${samples[@]}"; do
#  echo $sample 
  printf '\n%s\t' "$sample" >> $PARTITION_PATH/$OTUtbl_long_pre
  
  for taxon in ${taxa[@]}; do
#    echo $taxon
    grep -c -e "$taxon.*$sample" $PARTITION_PATH/$AllTaxIDTable | tr '\n' '\t' >> $PARTITION_PATH/$OTUtbl_long_pre   
  done
done


printf "\tOTU Table Complete!\n"

}

metadata (){
    ## Add Metadata
    ## add Sampling/collection date ($7), Site Name ($??), metadata barcode ($9), Sample ID ($NF) from metadata to OTU table
        ## Note this won't handle all files
        ## Use Find to gather all into one table?
   #rm $MTDA_PATH/AllSites.mtgnmDnaXtrct_tbl.temp
   #awk '{print $2}' $HMMR_PATH/$SampleList | sort -u > ~/NEON_processed/Pipeline/site.temp
   # while read i; do
   #  echo "$i"
   #  awk '{print $10, $6, $i}' $MTDA_PATH/"$i".mtgnmDnaXtrct.tbl >> AllSites.mtgnmDnaXtrct_tbl.temp
   # done<~/NEON_processed/Pipeline/site.temp
   ## NEEDS REDOING!!!##
   # awk '{print $NF, $7}' $MTDA_PATH/*.metadata2.tbl | sort -k 1 > $MTDA_PATH/AllSites.metadata_tbl.temp
   # tail -n+2 $HMMR_PATH/$OTUtbl_long | sort -k1 > $HMMR_PATH/AllSites_longOTU.temp
   # join --nocheck-order -1 1 -2 1 $MTDA_PATH/AllSites.metadata_tbl.temp $HMMR_PATH/AllSites_longOTU.temp > AllSites_long_otu_tbl.test 
 
## Build a version of the metadata.tbl for each site that mitigates differences in metadata file formats across sampling years   
 while read site_name; do
     
  paste -d ' ' $MTDA_PATH/"$site_name".metadata.tbl $MTDA_PATH/"$site_name".list > $MTDA_PATH/"$site_name".metadata2.tbl
   
 done<$HMMR_PATH/$SiteList

## Merge Partitioned OTU tables
#cat $PARTITION_PATH/AllSites_long.otu.* > $HMMR_PATH/$OTUtbl_long_pre

## Merge OTU table and metadata into new OTU table (without headers) 
 join -t $'\t' -a 2 -e "NA" -o auto -1 1 -2 1 \
     <(awk '{print $NF "\t" $7 "\t" $9}' $MTDA_PATH/*.metadata2.tbl | sort -k 1,1) \
     <(sort -k 1,1 $HMMR_PATH/$OTUtbl_long_pre) \
 > $HMMR_PATH/AllSites_longOTU.temp


## Add Metadata Column Names to U+Taxonomy Column headers
 awk 'BEGIN{a=0}{b=++a; printf "OTU_%d\t",b}END{printf "\n"}' $HMMR_PATH/$AllTaxa_list >> $HMMR_PATH/$AllTaxa_head                 
 
# awk 'BEGIN{a=0; printf "Raw_File_Name\tCollection_Date\tMetadata_Tag\tSite_Name\tPhylum\tClass\Family\tGenus\t"}END{printf "\n"}' $HMMR_PATH/$AllTaxa_list >> $HMMR_PATH/$AllTaxa_head                 
 
## Merge Metadata Column Name + OTU/Taxonomy Header with new OTU table to make FINAL OTU table
#cat $HMMR_PATH/$AllTaxa_head $HMMR_PATH/AllSites_longOTU.temp > $HMMR_PATH/$OTUtbl_long_final

## Remove metadata2.tbls (for clarity and storage space)
# rm $MTDA_PATH/*.metadata2.tbl

}

taxonomy () {
echo "Subdividing OTU table!"


# Pull out all X



    ## Merge Taxonomy with Raw OTU Table
    #awk '{print $2,$3,$4,$5,$6,$7,$8,$9}' $HMMR_PATH/AllSites.BH_txID.tbl | sort -u -k 7 > $HMMR_PATH/temp.taxonomy
    #join --nocheck-order -1 1 -2 7 <(sort -k 1 $HMMR_PATH/$OTUtbl_long) <(sort -k 7 $HMMR_PATH/temp.taxonomy | tr ' ' "\t") > $HMMR_PATH/$AllTaxonomyTable
    #join --nocheck-order -1 7 -2 1 <(awk '{print $2,$3,$4,$5,$6,$7,$8}' Xtract/SSU_BLAST/AllSites.BH_txID.tbl | sort -u -k 7 | tr ' ' "\t") <(sort -k 1 Xtract/SSU_BLAST/AllSites.otu) > $HMMR_PATH/$AllTaxonomyTable
    ## grep through taxonomy table for taxa from list
    ## sort new list
    ## cat new list to top of OTU table (double check that OTU table is sorted the same as the first one)
    ## couldn't you also just add the taxonomy when you first make the OTU table header?
}



usage (){
cat << EOF
Usage:
    -g gather_tables
    -t OTU_tall
    -l OTU_long
    -m Add Metadata
    -x Add Taxonomy
    -a The whole Enchilada
    -u Usage
EOF
}

while getopts 'gtlmxau' OPTION; do
    case "$OPTION" in
        g) 
           echo "Gathering tables..."
           gather_tables
           ;;
        t) 
           echo "Making a tall OTU table..."
           OTU_tall
           ;;
        l)
           echo "Making a long OTU table..." 
           OTU_long
           ;;
        m)    
           echo "Adding metadata..."
           metadata
           ;;
        x) 
           echo "Adding taxonomy..."
           taxonomy
           ;;
        a) 
           echo "The whole enchilada..."
           #default long OTU table
           gather_tables
           OTU_tall
           metadata
           taxonomy
           ;;
        u) 
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

#usage 

## OUTPUT PROGRAM RUN TIME
format_time() {
 ((h=${1}/3600))
 ((m=(${1}%3600)/60))
 ((s=${1}%60))
 printf "%02d:%02d:%02d\n" $h $m $s
}

printf "Script completed in $(format_time $SECONDS)\n"

exit


