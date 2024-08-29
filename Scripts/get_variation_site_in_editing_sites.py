#!/usr/bin/env python3
# Caculate the number of variation site in editing sites
import sys
input_file = open(sys.argv[1], 'r')
output_file = open(sys.argv[2], 'w+')
for line in input_file:
      if line.strip().split()[8] == "+":
            site = str(int(line.strip().split()[1]) - int(line.strip().split()[4]) + 1)
            output_file.write(site + "\n")
      else:
            site = str(int(line.strip().split()[5]) - int(line.strip().split()[1]) + 1)
            output_file.write(site + "\n")
input_file.close() 
output_file.close()
