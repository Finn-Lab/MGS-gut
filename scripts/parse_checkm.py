#!/usr/bin/env python2

import sys
import os

if len(sys.argv) < 3:
    print "ERROR! usage: python script.py bins_qa.tab bins_taxonomy.tab use_prefix"
    sys.exit()

try:
    if sys.argv[3]:
        if sys.argv[3] == "use_prefix":
            fields = os.path.abspath(sys.argv[1]).split("/")
            for n,i in enumerate(fields):
                if i == "metaspades":
                    store = n
            pref = os.path.abspath(sys.argv[1]).split("/")[store-1]
        else:
            pref = sys.argv[3]
except:
    pref = ""

bins = {}

with open(sys.argv[1], "rU") as f:
    for line in f:
        line = line.strip("\n")
        cols = line.split("\t")
        if cols[0] != "Bin Id":
            name = "%s%s" % (pref, cols[0])
            complet = float(cols[11])
            cont = float(cols[12])
            heter = float(cols[13])
            bins[name] = [complet, cont, heter]

with open(sys.argv[2], "rU") as f:
    for line in f:
        line = line.strip("\n")
        cols = line.split("\t")
        if cols[0] != "Bin Id":
            name = "%s%s" % (pref, cols[0])
            taxa = cols[3]
            try:
                bins[name].append(taxa)
            except:
                continue

for ele in bins.keys():
    print "%s_%s\t%.2f\t%.2f\t%.2f\t%s" % (pref, ele, bins[ele][0], bins[ele][1], bins[ele][2], bins[ele][3])
