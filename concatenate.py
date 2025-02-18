#!/bin/python3

#################
## SCRIPT INFO ##
#################

# This script concatenates the R1 and R2 fasta sequences in a list of fasta files passed at the cmd line 
# NOT CURRENTLY IN USE. Concatenation of raw NEON sequence files is not needed as I discovered how to merge them (they come unsorted and need to be sorted first).
    # This might be a useful script in the future, although I think there is a faster/more efficient way to accomplish it.
# NOT CALLED BY ANYTHING. Would be called by ~/NEON_processed/Pipeline/Main.sh
# Calls: none
# Last Update: 12/2/2022 by AR Thompson

#################

## PENDING UPDATES ##

## Part B is working as intended; however, header lines in paired R1 and R2 files from NEON do not match. I need to investigate that before proceeding with Part A.


## PROGRAM OUTLINE ##

# 1. Collate all fasta file names (paths?) into an R1 and R2 document, respectively (happens prior to program running). ORDER OF FILENAMES MUST MATCH EXACTLY BETWEEN R1 AND R2 FILES!
        ## I don't understand why this step matters ... 2/10/2023
# 2. Open the R1 and R2 files
# 2. Read in the names of each file (R1 and R2) and then add each name to a vector
# 3. Take the nth element of R1_file_vector and R2_file_vector and open the files those elements refer to.
# 4. Go through the open R1 and R2 files:
#    i. check that the seqIDs at the nth line in each file matches
#    ii. copy the seqID to the new file
#    iii. grab the next line from each file (should be the associated DNA Sequence)
#    iv. concatenate the two sequences, then output the sequence to the next line in the new file
#    v. repeat for the rest of the lines in the open R1 and R2 files
#    vi. close R1 and R2 files and open the next ones in the respective vectors
# 5. End Program


## SCRIPT ##

import sys
import os
import argparse

### A. LOAD FILENAMES FROM COMMAND LINE ###

parser=argparse.ArgumentParser()

parser.add_argument("input_file_R1", type=str, help="R1 Input file")
parser.add_argument("input_file_R2", type=str, help="R2 Input file")

args = parser.parse_args()

in_R1 = args.input_file_R1
in_R2 = args.input_file_R2

### B. CONCATENATE SEQUENCES ###

#i = 0
#while i < 1:
#while i < length(fastaList_Fwd):

#   in_R1 = fastaList_Fwd[i]                          # load nth element of FWD fasta list vector into R1 filename variable
#   in_R2 = fastaList_Rev[i]                          # load nth element of REV fasta list vector into R2 filename variable

# in_R1 = SITE + sample

 in_R1 = "ABBY_2017-06/BMI_Plate42WellG11_test_R1.fasta"
 in_R2 = "ABBY_2017-06/BMI_Plate42WellG11_test_fail_R1.fasta"

 with open(in_R1) as R1_file:               # automatically closes input file when block code completes
  R1_lines = R1_file.read().splitlines()    # reads all lines from files into list variable; removes '\n'
 with open(in_R2) as R2_file:
  R2_lines = R2_file.read().splitlines()
   
 if (len(R1_lines) != len(R2_lines)):         # WORKS!
  print ("Error: R1 and R2 files have differing number of lines. Program terminated.")
  sys.exit()

 R12_new_lines=list()
 j = 0
 while j < len(R1_lines):
   cnct_line = ""                                                              #WORKS! 
   line_R1 = R1_lines[j]                           # read line_r1 from R1_file #WORKS!
   line_R2 = R2_lines[j]                           # read line_r2 from R2_file #WORKS!

   if ((line_R1.find(">") == 0) and (line_R2.find(">") == 0)):   # if line is a header line... #WORKS!
   #  print ("Found the Header Lines!")
    if line_R1 != line_R2:                                                  #WORKS!
        print ("Error: Header lines do not match! Program terminated.")                           # header lines from R1 and R2 should match at this point; exit if not
        sys.exit()
    else:
        R12_new_lines.append(line_R1)           # add R1 header line to vector (header could come from either R1 or R2)
        print ("Moved to next line")                
        j += 1                                  # move to next line

   elif line_R1.find(">") != 1 and line_R2.find(">") != 1:  # this could be more robust, by maybe checking for ATCGNs? #WORKS!
   #  print ("Not Header lines!")
     cnct_line += line_R1         #WORKS!
     cnct_line += line_R2         #WORKS
     R12_new_lines.append(cnct_line)                     # !!QUESTION!! does vector.append() add an EOL character automatically? #No it does not # WORKS!
     j += 1                                      # move to next line
    
   else:
     print ("Error: Unknown Problem with lines from files caused program to terminate.")                              # there should only be cases where both lines have ">" or neither does. Any other case indicates a proble.
     sys.exit()

## Open R12 File and write output to it
 print(R12_new_lines)
 #out_R12 = in_R1 + "_R12.fasta"                            # I will probably have to do some more intricate parsing here
 #  out_R12 = "ABBY_2017-06/BMI_Plate42WellG11_test_R12.fasta"
 #  out_file = open(out_R12, 'w')                   # the new file name has to match the input file name
 #  out_file.write('\n'.join(R12_new_lines))                      # print header to out_file
 #  out_file.close()
   

 i += 1

sys.exit()



