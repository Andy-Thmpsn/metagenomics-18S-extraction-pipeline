#!/bin/bash

#SBATCH --job-name=NEON_main
#SBATCH --output=logs/mainNEON_%j.out
#SBATCH --error=logs/mainNEON_%j.err
#SBATCH --time=8:00:00 
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --mem=4GB 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=athompson@cmc.edu

export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

#################
## SCRIPT INFO ##
#################

# Manages the SSU_EXTRACT pipeline and all its subscripts (unpack.sh, preprocess.sh, pipeline.sh, postprcs.sh, filter.py, stats.sh, concatenate.py, and fstq2fsta.sh)
# Called by USER
# Calls: unpack.sh, preprocess.sh, pipeline.sh, postprcs.sh, and fstq2fsta.sh
# Last Update: 3/3/2023 BY AR Thompson

#################


## PENDING UPDATES ##

   # See Readme file (Guide2Pipeline.txt) for full list

   ## Try different databases (SILVA/NCBI) and compare to PR2
   ## Revisit my "Best Hit" algorithm. 
   ## Optimize hmmr and blastn search parameters
   ## Double check that all the steps of the pipeline are working as intended
   ## Make a separate script file for the metadata reformatting section (mostly for clarity)


## VARIABLES ##

HOME_PATH=/hopper/groups/yoosephlab/athompson/NEON_18S_Xtrct
PROCESSED_PATH=$HOME_PATH/NEON_processed
RAW_PATH=~$HOME_PATH/NEON_raw
LINKS_PATH=$HOME_PATH/DataLinks
HMMR_PATH=$PROCESSED_PATH/Xtract/SSU_BLAST  ## WARNING: Path name and variable_name match what is in ~/NEON_processed/Pipeline/pipeline.sh. DO NOT CHANGE.
MTDA_PATH=$HMMR_PATH/mtgnmDnaXtrcts
QC_PATH=$PROCESSED_PATH/QCtrl
TRMD_PATH=$HOME_PATH/NEON_trimmed
LOG_PATH=$HOME_PATH/LogFiles

scripts=$HOME_PATH/Pipeline


## INPUT SITE FILE ##

#SiteList=AllSiteNames_practice.list
#SiteList=AllSiteNames_legacy.list
#SiteList=AllSiteNames3.list
SiteList=AllSiteNames.list
#SiteList=legacy_files
#SiteList=Error_sites.list


## FUNCTIONS ##

get_data () {                ## DOWNLOAD DATA FROM NEON WEBSITE ## ALL OPERATIONS TESTED AND WORKING 7/27/2022
 
 if [ "$1" == "split" ]; then                   ## WORKS!!! 7/27/2022

  echo "Splitting data links text file into individual sites..."
  
  filename=FirstLook.txt
  
  dos2unix ${LINKS_PATH}/${filename}                                                                             ## Converts text files made in Notepad++ to linux formatting

  awk '{ if ( $1 ~ /.+_20[12].{1}-/){head=$1} else {print > head".file"}}' ${LINKS_PATH}/$filename              ## Splits notepad++ combined NEON link file into separate files that can be run through wget
  
  mv $HOME_PATH/*_*-*.file $LINKS_PATH/.


 elif [ "$1" == "extract" ]; then

    echo "Extracting data links from metadata.tbl files..."

    while read name;       ### WARNING: be careful with this. You need to manually select the line you want below (for legacy or non files) and make sure it doesn't overwrite the other type of file already there. In reality, you shouldn't need to use this again...
    do
        awk '{print $(NF-1)}' $MTDA_PATH/og_files/${name}.rawDataFiles.tbl | sort -u > ${name}.file         #for "legacy" files
        #awk '{print $NF}' $i | sort -u > ${i/.rawDataFiles.tbl/.file}         #forregular files
    mv ${name}.file $LINKS_PATH/.
  done<$HMMR_PATH/$SiteList


 elif [ "$1" == "get" ]; then          ## NEEDS TO CHECK FOR LIST FIRST; IF NONE, THEN IT SHOULD CALL pre-process("list") ## WORKS!!! 7/20/2022
  
  echo "Downloading raw data files en masse..."
  
  while read i; do
  
   echo "Site Name: $i"
   mkdir -p $RAW_PATH/$i/   
   wget -i $LINKS_PATH/"$i".file -P $RAW_PATH/$i/
  
  #done< <(ls $LINKS_PATH/*.file | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1)
  done<$HMMR_PATH/$SiteList
 

 elif [ "$1" == "one" ]; then          ## NEEDS TO CHECK FOR LIST FIRST; IF NONE, THEN IT SHOULD CALL pre-process("list") ## WORKS!!! 7/20/2022
  
  read -p "Site Name: " site_name
  
  echo "Downloading raw data files for $site_name..."
  mkdir -p $RAW_PATH/${site_name}/   
  wget -i $LINKS_PATH/${site_name}.file -P $RAW_PATH/${site_name}/
 

 fi

}

unpack () {

 ## 1. UNPACK ALL COMPRESSED FILES (*.TAR.GZ, *.GZ, *.ZIP)
 if [ "$1" == "all" ]; then
  echo "Unpacking All"
  
  while read site_name;
   do
  echo "Site Name: ${site_name}"
  #TAR
   if compgen -G "${RAW_PATH}/${site_name}/*fastq.tar.gz"; then         # Compgen returns true if specified extension is found in path
   echo "Unpacking *tar.gz files"
   sbatch $scripts/unpack.sh -t $site_name
  #GZ
   elif compgen -G "${RAW_PATH}/${site_name}/*.fastq.gz"; then
     echo "Unpacking *.gz files"
     sbatch $scripts/unpack.sh -g $site_name
  
  #ZIP
   elif compgen -G "${RAW_PATH}/${site_name}/*.zip"; then
  #  unpack_zip $site_name
    sbatch $scripts/unpack.sh -z ${site_name}
   echo "not today georgie!"
   fi
  done<$HMMR_PATH/$SiteList
 

 elif [ "$1" == "manage" ]; then

  while read site_name; 
  do
   ls $RAW_PATH/${site_name}/*_R*.fastq | sed -r 's/(.*_)R[12]*.fastq/\1/g' | cut -d '/' -f 8 | sed 's/.$//' | sort -u > $MTDA_PATH/${site_name}_temp1.list

   > $MTDA_PATH/${site_name}_to.remove              # clears output file

   while read list;
   do
    sed 's/${list}//g' $MTDA_PATH/${site_name}_temp1.list >> $MTDA_PATH/${site_name}_to.remove  # this should output only the files that don't match the target fastqs, but its not working with the variable for some reason
#    grep -v "$list" $MTDA_PATH/${site_name}_temp1.list >> $MTDA_PATH/${site_name}_to.remove  # this should output only the files that don't match the target fastqs, but its not working with the variable for some reason
   done<$MTDA_PATH/${site_name}.list
  
   #rm $MTDA_PATH/${site_name}_temp1.list

 #  while read remove_name;
 #   do
 #    rm $RAW_PATH/${site_name}/${remove_name}*fastq

 #   done<$MTDA_PATH/${site_name}_to.remove

  done<$HMMR_PATH/$SiteList

 else

   read -p "Site Name: " site_name

   if [ "$1" == "tar" ]; then sbatch $scripts/unpack.sh -t "$site_name"
   elif [ "$1" == "gz" ]; then sbatch $scripts/unpack.sh -g  "$site_name"
   elif [ "$1" == "zip" ]; then sbatch $scripts/unpack.sh -z "$site_name"

   else echo "No matching input parameters found"; usage

   fi

 fi
}




pre_process () {
 
 ## 1. CREATE MASTER LIST OF SITE NAMES            ## WORKS!!! 7/20/2022
  ## Make list of names from the names of the original metadata files downloaded from NEON, then send it to a text file
 if [ "$1" == "sites" ]; then
  echo "Making a List, Checking it Twice ..."
  cd $MTDA_PATH/og_files/
  ls *.csv | cut -d '.' -f 1 | sort -u | cat > $HMMR_PATH/$SiteList 
    

 ## 2. RENAME METADATA FILE NAMES                            ##WORKS!! 7/20/2022
 elif [ "$1" == "rename" ]; then
  echo "Changing NEON File Names..."
  cd $MTDA_PATH/
  for k in NEON*.csv
   do mv "$k" "$(echo "$k" | sed -E 's/NEON\..+\.([A-Za-z]+)\..+\.mms_([A-Za-z]+)\.([0-9]{4}-[0-9]{2})\..+\.csv/\1_\3.\2.csv/g')"
  done

 ## 3. CONVERT NEON METADATA FILES TO HUMAN READABLE FORMAT
    #metagenomeDnaExtraction
    #RawDataFiles

 elif [ "$1" == "reformat" ]; then
  echo "Formatting NEON Metadata..."

 #metagenomeDnaExtraction
  for i in $MTDA_PATH/*.metagenomeDnaExtraction.csv; do          
    echo "MetagenomeDnaExtraction File: "$i

 ## Add Header to New metagenomeDnaExtraction File (.tbl)
    head -1 $i | tr ' ' '_' | tr -d '"' | tr ',' '\t' > ${i/.csv/.tbl}                        

 ## Filter for "metagenome" samples and add to new file (.tbl)
    cat $i | tr ' ' '_' | tr -d '"' | tr ',' '\t' | grep -P "\tmetagenomics\t" >> ${i/.csv/.tbl} #   
  
  done

  #RawDataFiles
  for q in $MTDA_PATH/*.rawDataFiles.csv  ## TESTED/WORKS ###         # Clean this up when you get a chance (e.g., split the join command so it's easier to parse, put comments etc.)
   do
   echo "RawDataFile: "$q

   ## Modern Data (post 2015)   
   #Convert to human readable, make list file
   cat $q | tr ' ' '_' | tr -d '"' | tr ',' '\t' | tee ${q/.csv/.tbl} | awk '{print $11, $12}' | tail -n+2 | sed 's/_R[0-9].*[fastq].*//g' >  ${q/.rawDataFiles.csv/.list}

   # Merge List and metagenomeDnaExtraction.tbl to create metadata.tbl
   join -1 10 -2 1 <(sort -k10,10 ${q/.rawDataFiles.csv/.metagenomeDnaExtraction.tbl}) <(sort -u -k1,1 ${q/.csv/.list}) > ${q/.rawDataFiles.csv/.metadata.tbl}
 
   ## Legacy Data (pre 2016)  ## Working 9/20/2022 ##
#   cat $q | tr ' ' '_' | tr -d '"' | tr ',' '\t' | tee ${q/.csv/.tbl} | awk '{print $11, $10}' | tail -n+2 | sed 's/-Sequences.zip//g' >  ${q/.rawDataFiles.csv/.list}
#   join -1 10 -2 2  <(sort -k10,10 ${q/.rawDataFiles.csv/.metagenomeDnaExtraction.tbl}) <(sort -u -k2,2 ${q/.rawDataFiles.csv/.list}) > ${q/.rawDataFiles.csv/.metadata.tbl}
  
#    awk '/BMI_*Plate/ {print $0}' ${q/.rawDataFiles.csv/.metadata.tbl} > ${q/.rawDataFiles.csv/.list}
  done

 
## 4. ??????
 elif [ "$1" == "samples" ]; then 
 
  while read site_name;
  do
   #echo "$site_name"
 
   # Form sample list from trimmed R1 files  8/18/2022   
   ls $TRMD_PATH/${site_name}/*_R1.paired.fasta | sed -r 's/(.*_)R1.paired.fasta/\1/g' | cut -d '/' -f 8 | sort -u > $MTDA_PATH/${site_name}_temp1.list
    ls -p $RAW_PATH/${site_name}/*.fastq | sed -r 's/(.*)_R[12]*.fastq/\1/g' | cut -d '/' -f 8 | sort -u > $MTDA_PATH/${site_name}_temp1.list              # form list from raw unzipped files 
    > $MTDA_PATH/${site_name}.list


## Legacy Data ONLY 
   while read metadata_name;
   do
    # Compare fastq file names in the $RAW_PATH/$SITE/ folder to the 'metagenome' ones from ... (?)   
    grep -e "${metadata_name}" $MTDA_PATH/"$site_name"_temp1.list >> $MTDA_PATH/"$site_name".list  
    #  the metadata files to make a final list of names that can be used to process the fastq files. The list 

   done< <(awk '{print $NF}' $MTDA_PATH/"$site_name".metadata.tbl)
   done<$HMMR_PATH/$SiteList
   #rm $MTDA_PATH/"$site_name"_temp1.list 

exit
## Modern Data
 # This compares the list of fastq file names in the $RAW_PATH/$SITE/. folder to the 'metagenome' ones from  the metadata files to make a final list of names that can be used to process the fastq files. The list is made twice, first by comparing column $1, and then column $NF. Depending on the site, the information is in either of those columns (it appears the formatting was not consistent over time or across sites).   
   while read metadata_name; do
    
     # echo "$metadata_name"
    
     grep -e "${metadata_name}_" $MTDA_PATH/"$site_name"_temp1.list >> $MTDA_PATH/"$site_name"_temp2.list
   
   done< <(awk '{print $1}' $MTDA_PATH/"$site_name".metadata.tbl)

   while read metadata_name;   do

     # echo "$metadata_name"
    
     grep -e "${metadata_name}_" $MTDA_PATH/${site_name}_temp1.list >> $MTDA_PATH/"$site_name"_temp2.list

   done< <(awk '{print $NF}' $MTDA_PATH/"$site_name".metadata.tbl)
   
     
   #rm $MTDA_PATH/"$site_name"_temp1.list 

#  done<$HMMR_PATH/$SiteList

exit

for i in $MTDA_PATH/*_temp2.list
   do
     # echo $i 
       sed 's/.$//' $i | sort -u > ${i/_temp2.list/.list}                        # This removes any duplicate namess in the $MTDA_PATH/*.list files and trims the "_" from their ends.

       rm $i
   done

  #mv $MTDA_PATH/*.csv $MTDA_PATH/og_files/.
  #mv $MTDA_PATH/*.rawDataFiles.tbl $MTDA_PATH/og_files/.
  #mv $MTDA_PATH/*.metagenomeDnaExtraction.tbl $MTDA_PATH/og_files/.

else echo "No matching input parameters found"; usage
fi

}

fastq2fasta () {

 ## 1. RUN FSTQ2FSTA.SH FOR ALL SAMPLES IN ALL SITES
 if [ "$1" == "all" ]; then
 
  echo "Converting files from fastq to fasta format in all sites ..."
 
  cp $scripts/fstq2fsta.sh ~/NEON_processed/Backups/fstq2fsta.sh
 
  while read site_name; do
   echo "Site Name: $site_name"  
 
   mkdir -p $LOG_PATH/fstq2fsta/$site_name/
  
   echo "Submitting Site: " $site_name >> file.log
   
   while read sample_name; do
   
    echo "Submitting: " $sample_name >> file.log
  
    sbatch --export=SITE=$site_name,NAME=$site_name $scripts/fstq2fsta.sh
  
   done<$MTDA_PATH/"$site_name".list
 
  done<$HMMR_PATH/$SiteList                        ## !/!/!/!/!/ I'll have to check for this input

 
 ## 2. RUN FSTQ2FSTA.SH FOR ALL SAMPLES IN A SINGLE, USER SPECIFIED SITE 
 elif [ "$1" == "one" ]; then
  
  read -p "Site Name: " site_name
 
  echo "Converting fastq to fasta format for $site_name ..."
  cp $scripts/fstq2fsta.sh ~/NEON_processed/Backups/fstq2fsta_backup.sh

  mkdir -p $LOG_PATH/fstq2fsta/$site_name/
 

  while read new_name; do

   echo "Submitting: " $new_name
  
    sbatch --export=SITE=$site_name,NAME=$new_name $scripts/fstq2fsta.sh

  done<$MTDA_PATH/"${site_name}".list

 else echo "No matching input parameters found"; usage
 
 fi

}


main_pipeline () {

## 1. CHANGE SITE_NAME IN PIPELINE.SH (Performs QC, Trim, Merge, hmmer, blastn, etc.)      ## /!/!/!/!/!/ Not sure how to make this one work yet

## 1. RUN PIPELINE.SH FOR ALL SAMPLES IN ALL SITES
 if [ "$1" == "all" ]; then
 
  echo "Running pipeline for all samples in all sites ..."
  
  while read site_name; do
  
     echo "Site Name: $site_name"  
  
     mkdir -p $LOG_PATH/Pipeline/$site_name/
  

    while read sample_name; do
   

      sbatch --export=SITE=$site_name,NAME=$sample_name $scripts/pipeline.sh
  

    done<$MTDA_PATH/"$site_name".list
 

  done<$HMMR_PATH/$SiteList                        ## !/!/!/!/!/ I'll have to check for this input


 ## 2. RUN PIPELINE.SH FOR ALL SAMPLES IN A SINGLE SPECIFIED SITE (pipeline.sh performs QC, Trim, Merge, hmmer, blastn, etc.)
 elif [ "$1" == "one" ]; then
  
  read -p "Site Name: " site_name
 
  echo "Running pipeline for $site_name ..."

  mkdir -p $LOG_PATH/Pipeline/$site_name/

  while read new_name; do

   echo "Submitting: " $new_name

   sbatch --export=SITE=$site_name,NAME=$new_name $scripts/pipeline.sh

  done<$MTDA_PATH/"${site_name}".list


 ## 3. RUN FASTQC WITH PREPROCESS.SH FOR ALL SAMPLES (preprocess.sh performs QC, Trim, Merge, and raw FASTQ file sorting)

 elif [ "$1" == "fastqc" ]; then

  read -p "pre or post: " option

  echo ${option}

  while read site_name; do

      while read sample_name; do

#       sed -i -E "s/SITE=[[:alnum:]]+_[[:digit:]]{4}-[[:digit:]]{2}/SITE=$site_name/" $scripts/preprocess.sh # Change $SITE variable in preprocess.sh
 
        sbatch --export=SITE=$site_name,NAME=$sample_name $scripts/preprocess.sh "-f$option"

      done<$MTDA_PATH/"$site_name".list

  done<$HMMR_PATH/$SiteList


 ## 4. RUN FLASH WITH PREPROCESS.SH TO MERGE SAMPLES IN ALL SITES (preprocess.sh performs QC, Trim, Merge, and raw FASTQ file sorting)

 elif [ "$1" == "merge" ]; then

  read -p "Pre or post (trimming): " option

  echo ${option}

  while read site_name; do

     while read sample_name; do

#       sed -i -E "s/SITE=[[:alnum:]]+_[[:digit:]]{4}-[[:digit:]]{2}/SITE=$site_name/" $scripts/preprocess.sh # Change $SITE variable in preprocess.sh
        
        sbatch --export=SITE=$site_name,NAME=$sample_name $scripts/preprocess.sh "-m$option"

     done<$MTDA_PATH/"$site_name".list

  done<$HMMR_PATH/$SiteList


 ## 5. RUN TRIMMOMATIC WITH PREPROCESS.SH FOR ALL SAMPLES (preprocess.sh performs QC, Trim, Merge, and raw FASTQ file sorting)

 elif [ "$1" == "trim" ]; then

  while read site_name; do

    while read sample_name; do

        sbatch --export=SITE=$site_name,NAME=$sample_name $scripts/preprocess.sh "-t"

    done<$MTDA_PATH/"$site_name".list

  done<$HMMR_PATH/$SiteList

 else echo "No matching input parameters found"; usage
 
 fi



}

super_tables ()  {
 
 echo "Creating Site-Specific Best Hit Summary Tables!"

 ## BEWARE: Suffix and variable_name match what is in ~/NEON_processed/Pipeline/pipeline.sh. DO NOT CHANGE.
 SFX_blstBHtax=.hmmblst.BH_txID.tbl
 PFX_ALL=All_          ##THREE STRINGS IN A ROW NOT WORKING; CHANGE MANUALLY  ## BEWARE: Prefix and variable_name match what is in ~/NEON_processed/Pipeline/pipeline.sh. DO NOT CHANGE.

## 7. MAKE COMBINED BH_txID.tbl FOR $SITE AND ADD $SITE NAME TO END OF LINE

  #printf "Checking to see if a combined BH_txID.tbl currently exists.\n"

  while read site_name; do

    if [ -f "$HMMR_PATH/${site_name}/All_${site_name}${SFX_blstBHtax}" ]; then
     rm "$HMMR_PATH/${site_name}/All_${site_name}${SFX_blstBHtax}"
    fi
  
    cat $HMMR_PATH/${site_name}/*"$SFX_blstBHtax" | awk -v s=${site_name} -v OFS="\t" '{print $0, s}' > $HMMR_PATH/${site_name}/All_${site_name}${SFX_blstBHtax}

  done<$HMMR_PATH/$SiteList 

## GATHER ALL SITE-SPECIFIC SUMMARY TABLES INTO ONE SUPER TABLE

}


OTU_tbl () {

 FILE=postprcs.sh
 LOGFILE=buildOTU.log

 > $LOGFILE

 ## Build super table, create reference docs (taxa list and sample list) and taxonomic header for the final OTU table
 
  $scripts/postprcs.sh "-a"
 
 ## Create copies of reference docs for each partition, then run the OTU table building subroutine in postprcs.sh on each partition seperately
 for i in $HMMR_PATH/Partitions/AllSites.samples.*; do

  printf "********************************************************\n" >> $LOGFILE
  printf "New File: \n${i}\n\n" >> $LOGFILE

  A="$(cut -d'.' -f3 <<< ${i})" 
  A=".$A"
  B="$(cut -d'/' -f8 <<< ${i})"

  printf "Suffix: ${A}\n" >> $LOGFILE
  printf "File Name: ${B}\n\n" >> $LOGFILE
  
  # sed -i -E "s/SampleList=AllSites.samples\.*[[:alnum:]]{0,2}/SampleList=${B}/g" $scripts/${FILE}

  cp $HMMR_PATH/AllSites.BH_txID.tbl $HMMR_PATH/Partitions/AllSites.BH_txID.tbl${A}
  cp $HMMR_PATH/AllSites_list.taxa $HMMR_PATH/Partitions/AllSites_list.taxa${A}

  sbatch --export=SampleList_Part=${B},OTUtbl_long_Part=AllSites_long.otu${A},AllTaxa_list_Part=AllSites_list.taxa${A},AllTaxIDTable_Part=AllSites.BH_txID.tbl${A} $scripts/${FILE} "-l"

   # sbatch $scripts/${FILE} -l 

 done

}


multi_qc () {

 if [ "$1" == "pre" ]; then

    echo "Running MultiQC on project fastqc files PRE trimming..."
    cd $QC_PATH/Pre/
    multiqc .

 elif [ "$1" == "post" ]; then
    echo "Running MultiQC on project fastqc files POST trimming..."
    cd $QC_PATH/Post/
    multiqc .

 else echo "No matching input parameters found"; usage
 fi

} 



usage (){

cat << EOF

Usage:

    -d get_Data
       split ~ prepare wget links from master file
       extract ~ pull links from *.rawDataFile.tbl files
       get ~ download data from links
       one ~ download data from a user specified file

    -u unpack (requires an argument)
       all ~ UNPACKS ALL ZIPPED files at once, whether they are .tar, .gz, or .zip
       tar ~ UNPACK MULTIPLE *.tar.gz FILES AT ONCE
       gz ~ UNPACK MULTIPLE *.gz FILES AT ONCE
       zip ~ UNPACK MULTIPLE *.zip FILES AT ONCE
       manage ~ remove fastq files that are not metagenome files (to save space)

    -s subdir  Manage fastq files in subdirectories to parent folder for downstream processing

    -p pre-Process (requires an argument)
       sites ~ CREATE LIST SITES TO PROCESS  [!!//!!NOTE: THIS IS DEPENDENT ON 'samples' BEING RUN FIRST]
       rename ~ RENAME .csv METADATA FILES 
       reformat ~ PROCESS METADATA FILES
       samples ~ CREATE LIST OF usamples LES TO PROCESS PER SITE [!!//!!NOTE: THIS IS DEPENDENT ON 'reformat' BEING RUN FIRST]

    -f fastq to fasta 
       one ~ convert fastq files to fasta format in one site
       all ~ convert fastq files to fast format in ALL sites in Site names list in HMMR_PATH/
    
    -m run Main pipeline (requires an argument)
       one ~ run the full pipeline on one site (program will request site_name input after it is run)
       all ~ run the full pipeline on ALL sites in Site names list in HMMR_PATH/
       fastqc ~ run only the fastqc part of the pipeline (runs on all samples/sites provided)
       merge ~ merge samples with FLASH (note that this sorts samples first; NEON sequence order in R1 and R2 files do not match)
       trim ~ trim samples with trimmomatic

    -t build Super Tables
        First gathers all best hit tables from all samples in each site into one, site-specific table, then ...
        Second collates these summary tables from all sites into one mega table comprising all sites and all samples.

    -o build OTU table (long) through partitioning
    
    -q run multiqc on fastqc files
        pre ~ run multiqc on "pre" trimmed fastqc files
        post ~ run multiqc on "post" trimmed fastqc files

    -h Usage

EOF

}


while getopts 'd:u:sp:f:m:toq:h' OPTION; do
    case "$OPTION" in
        d) 
           echo "Getting Data"
           get_data "$OPTARG"
           ;;
        u) 
           echo "Unpacking Compressed Files..."
           unpack "$OPTARG"
           ;;
        s)
           echo "Managing fastq files in subdirectories..."
           echo "Site Name: "
           read "$site_name"
           manage_subds "$site_name"
          ;;
        p) 
           echo "Preparing Data for Processing"
           #arg1="$OPTARG"
           #echo $arg1
           pre_process "$OPTARG"           
           ;;
        f)
           echo "converting Fastq to Fasta"
           fastq2fasta "$OPTARG"
           ;; 
        m) 
           #SITE=TOOL_2018-08
           #echo "The Site provided is $OPTARG"
           main_pipeline "$OPTARG"
           ;;
        t)    
           echo "Building a Super Table now..."
           super_tables
           #sbatch $scripts/postprcs.sh
           ;;
        
        o) 
           OTU_tbl
           ;;
        q)
           multi_qc "$OPTARG"
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
