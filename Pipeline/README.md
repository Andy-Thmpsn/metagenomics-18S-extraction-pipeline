# metagenomics-18S-extraction-pipeline
LAST UPDATED: 2/20/2025
LAST UPDATED: 7/6/2023

### OVERVIEW OF SCRIPTS

main.sh
preprocess.sh
    unpack.sh
    fstq2fsta.sh
pipeline.sh
postprcs.sh
    concatenate.sh
    stats.sh

### SCRIPT USAGE

 ##Main.sh##
    #Runs all other scripts
    #Called by None
    #Calls: All

    Usage:
    -d get_Data         ## DOWNLOAD DATA FROM NEON WEBSITE
        split ~ prepare wget links from master file
        extract ~ pull links from *.rawDataFile.tbl files
        get ~ download data from links
        one ~ download data from a user specified file

    -u unpack           ## DECOMPRESSES DOWNLOADED NEON FILES
        all ~ Decompresses ALL ZIPPED files at once, whether they are .tar, .gz, or .zip
        tar ~ Decompresses MULTIPLE *.tar.gz files AT ONCE
        gz ~ Decompresses MULTIPLE *.gz files AT ONCE
        zip ~ Decompresses MULTIPLE *.zip files AT ONCE
        manage ~ remove fastq files that are not metagenome files (to save space)

    -s subdir

    -p pre-Process      ## Prepares sample files for downstream processing
        reformat ~ PROCESS METADATA FILES
        samples ~ CREATE LIST OF usamples LES TO PROCESS PER SITE [!!IMPORTANT: THIS IS DEPENDENT ON 'reformat' BEING RUN FIRST]
        sites ~ CREATE LIST OF SITES TO PROCESS  [!!IMPORTANT: THIS IS DEPENDENT ON 'samples' BEING RUN FIRST]
        rename ~ RENAME .csv METADATA FILES 

    -f fastq to fasta   ## CONVERT FASTQ FORMAT TO FASTA
        one ~ convert fastq files to fasta format in ONE SITE
        all ~ convert fastq files to fast format in ALL SITES in Site names list in HMMR_PATH/

    -m run processing pipeline  ## 
        one ~ run the full pipeline on ONE SITE (program will request site_name input after it starts)
        all ~ run the full pipeline on ALL SITES in Site names list in HMMR_PATH/
        fastqc ~ run only the fastqc part of the pipeline (runs on all samples/sites provided)
        merge ~ merge samples with FLASH (note that this sorts samples first; NEON sequence order in R1 and R2 files do not match)
        trim ~ trim samples with trimmomatic

    -t build Super Tables
        First, gathers ALL BEST HIT TABLES from ALL SAMPLES in each site into one, site-specific table
        Second, collates summary tables from ALL SITES into one MEGA TABLE comprising ALL SITES and ALL SAMPLES.
    
    -o build OTU table (long) through partitioning

    -q run multiqc on fastqc files
        pre ~ run multiqc on "pre" trimmed fastqc files
        post ~ run multiqc on "post" trimmed fastqc file

    -u Usage

        unpack.sh
        fstq2fsta.sh

 ##preprocess.sh##

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


 ##pipeline.sh##

    # Runs the core steps of the SSU_EXTRACT pipeline, namely 
        sequence extraction and identification (nhmmer and blastn)
    # Optionally runs Metaxa2, TaxaTarget, and EukDetect
    # Called by ~/NEON_processed/Pipeline/Main.sh
    # Calls: filter.py
    # Last Update: 3/3/2023 BY AR Thompson

    # Part 2 & 3 - rDNA extraction and Identification
    # mem > 24GB for a few, <20GB for most
    # <3 hours for most, <6 hours max
    # ntasks < 2



 ##postprcs.sh##

    # Manages the postprocessing options (build super BH table, build OTU tables [tall and long], and add metadata to the OTU table ) for the SSU_EXTRACT pipeline.
    # Called by ~/NEON_processed/Pipeline/Main.sh
    # Calls: none
    # Last Update: 3/3/2023 BY AR Thompson


    Usage:
        -g gather_tables
        -t OTU_tall
        -l OTU_long
        -m Add Metadata
        -x Add Taxonomy
        -a The whole Enchilada
        -u Usage

    concatenate.py

    stats.sh


