#!/usr/bin/env python2
import sys

if len(sys.argv) < 2:
    print "ERROR! usage: python script.py dist.tab"
    sys.exit()

hits = {}

with open(sys.argv[1], "r") as f:
    for line in f:
        cols = line.strip("\n").split("\t")
        if len(cols) == 5: 
            ref = cols[0]
            mag = cols[1]
            dist = float(cols[2])
            if mag not in hits.keys():
                hits[mag] = [ref,dist]
            else:
                if hits[mag][1] > dist:
                    hits[mag] = [ref,dist]

for ele in hits.keys():
    print "%s\t%s\t%.2f" % (ele, hits[ele][0], hits[ele][1])
    
