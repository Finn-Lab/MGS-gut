#!/usr/bin/env python2

import os
import sys

if len(sys.argv) < 2:
    print "usage: python script.py taxonomy.tab (up to genus rank)"
    sys.exit()

ranks = ["Phylum", "Class", "Order", "Family", "Genus"]
tax_counts = {}

for r in ranks:
    tax_counts[r] = {}

with open(sys.argv[1], "rU") as f:
    linen = 0
    for line in f:
        linen += 1
        if linen > 1:
            line = line.strip("\n")
            cols = line.split("\t")
            # check if taxonomy file relates to UMGS or HGR
            if "UMGS" in cols[0]:
                tax = line.split("\t")[2:]
            else:
                tax = line.split("\t")[1:]
            # remove NAs
            for n,i in enumerate(tax):
                if i == "NA":
                    tax = tax[:n]
            # go through each lineage and count taxa
            for n,i in enumerate(tax):
                if n == 0:
                    i = "p__"+i
                    if i not in tax_counts["Phylum"].keys():
                        tax_counts["Phylum"][i] = 1
                    else:
                        tax_counts["Phylum"][i] += 1
                elif n == 1:
                    i = "c__"+i
                    if i not in tax_counts["Class"].keys():
                        tax_counts["Class"][i] = 1
                    else:
                        tax_counts["Class"][i] += 1
                elif n == 2:
                    i = "o__"+i
                    if i not in tax_counts["Order"].keys():
                        tax_counts["Order"][i] = 1
                    else:
                        tax_counts["Order"][i] += 1
                elif n == 3:
                    i = "f__"+i
                    if i not in tax_counts["Family"].keys():
                        tax_counts["Family"][i] = 1
                    else:
                        tax_counts["Family"][i] += 1
                elif n == 4:
                    i = "g__"+i
                    if i not in tax_counts["Genus"].keys():
                        tax_counts["Genus"][i] = 1
                    else:
                        tax_counts["Genus"][i] += 1

# print final output
for d in tax_counts.keys():
    for t in tax_counts[d].keys():
        print "%s\t%s\t%i" % (d, t, tax_counts[d][t])
print "Kingdom\tk__Bacteria\t%i" % (linen-1)
