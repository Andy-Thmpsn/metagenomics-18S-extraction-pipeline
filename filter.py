import sys
import os
import argparse

#################
## SCRIPT INFO ##
#################

# Filters nhmmr output table for sequences with % IDs and lengths below supplied cutoff values (see pipeline.sh).
# Called by ~/NEON_processed/Pipeline/pipeline.sh
# Calls: Nothing.
# Last Update: 3/3/2023 BY AR Thompson

#################



parser=argparse.ArgumentParser()

parser.add_argument("input_file", type=str, help="Input file to be processed")
parser.add_argument("Min_Identity", type=float, help="Sequence ID minimum cutoff")
parser.add_argument("Min_Length", type=int, help="Sequence minimum length cutoff")

args = parser.parse_args()

File = args.input_file
IDCutoff = args.Min_Identity
Len_Cutoff = args.Min_Length

with open(File) as f:
  for line in f:
      line_values = line.strip().split()
      if len(line_values) >= 4:                    #Checks for 3 or greater columns, to avoid throwing an error
        if float(line_values[2]) >= IDCutoff and int(line_values[3]) >= Len_Cutoff:    #Filters sequence match identities below 93%
           print(line, end="")
#        if int(line_values[3]) >= Len_Cutoff:
#            print(line, end="")
