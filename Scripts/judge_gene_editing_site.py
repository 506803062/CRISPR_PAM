#!/usr/bin/env python3
# Whether gene has editing sites
import sys
input_file1 = open(sys.argv[1], 'r')
input_file2 = open(sys.argv[2], 'r')
output_file = open(sys.argv[3], 'w+')
gene_list1 = []
gene_list2 = []
for line1 in input_file1:
      gene_list1.append(line1.strip().split()[3])
gene_list1 = list(set(gene_list1))
for line2 in input_file2:
      gene_list3 = line2.strip().split()[9].split(',')
      gene_list2.extend(gene_list3)
gene_list2 = list(set(gene_list2))
for a in gene_list1:
      if a in gene_list2:
            output_file.write(a + "\t" + "Yes" + "\n")
      else:
            output_file.write(a + "\t" + "No" + "\n")
input_file1.close()
input_file2.close() 
output_file.close()
