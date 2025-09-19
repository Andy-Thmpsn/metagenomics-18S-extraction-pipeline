import sys
import os
import argparse

#################
## SCRIPT INFO ##
#################

# Filters nhmmr output table for sequences with % IDs and ?lengths and query coverage below supplied cutoff values (see pipeline.sh).
# Called by ~/NEON_processed/Pipeline/pipeline.sh
# Calls: Nothing.
# Last Update: 3/3/2023 BY AR Thompson

#################



parser=argparse.ArgumentParser()

parser.add_argument("input_file", type=str, help="Input file to be processed")
parser.add_argument("Min_Identity", type=float, help="Sequence ID minimum cutoff")
parser.add_argument("Min_Length", type=int, help="Sequence minimum length cutoff")
parser.add_argument("Query_Coverage", type=float, help="Minimum query coverage cutoff")

args = parser.parse_args()

File = args.input_file
IDCutoff = args.Min_Identity
Len_Cutoff = args.Min_Length
Query_Cutoff = args.Query_Coverage

exit

with open(File) as f:
  for line in f:
      line_values = line.strip().split()
      if len(line_values) >= 4:                    #Checks for 4 or greater columns, to avoid throwing an error
         #Filters sequence match identities below 93%
        if ( 
           float(line_values[2]) >= IDCutoff and 
           int(line_values[3]) >= Len_Cutoff and
           int(line_values[8]) - int(line_values[7]) >= (int(line_values[3]) * Query_Cutoff)
        ):
           print(line, end="")      
