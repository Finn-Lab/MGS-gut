#!/usr/bin/env python2

import os
import sys
import math
import numpy

if len(sys.argv) < 3:
    print "usage: script.py bwa_depth.tab bwa_depth-pos.tab"
    sys.exit()

name = os.path.basename(sys.argv[1]).split("_")[0]
spec_stats = {}

# recover genome length and total counts per species (for average depth)
with open(sys.argv[1], "r") as f:
    for line in f:
        if line[0] != "*": # exclude unmapped line
            cols = line.strip("\n").split("\t")
            species = "_".join(cols[0].split("_")[:-1]) # species name
            clength = int(cols[1]) # contig length
            counts = int(cols[2]) # read counts
            if species not in spec_stats.keys(): # first time species is found
                spec_stats[species] = [clength, counts, []] # dict of length, counts and covPos
            else:
                spec_stats[species][0] += clength
                spec_stats[species][1] += counts

# recover covered positions (for calc mean depth, coverage and evenness)
with open(sys.argv[2], "r") as f:
    for line in f:
        if line[0] != "*": # exclude unmapped line
            cols = line.strip("\n").split("\t")
            species = "_".join(cols[0].split("_")[:-1]) # species name
            pos = int(cols[1]) # covered position
            depth = int(cols[2]) # depth of covered position
            spec_stats[species][2].append(depth)

print "Genome\t%s_Length\t%s_Counts\t%s_MeanDepth\t%s_Coverage\t%s_CoeffVar" % (name, name, name, name, name)

# combine stats and print per species
for species in spec_stats.keys():
    length = spec_stats[species][0]
    counts = spec_stats[species][1]
    covPos = spec_stats[species][2]
    covBases = len(covPos)
    coverage = float(covBases)/length*100
    #if covBases < length: # include zeroes
    #    missing = length - covBases
    #    covPos.extend([0]*missing)
    if covBases < 1:
        covSd = 0
    else:
        covSd = numpy.std(covPos)
    meanDepth = float(sum(covPos))/max(len(covPos),1)
    if float(meanDepth) > 0:
        cV = covSd/meanDepth
    else:
        cV = 0
    print "%s\t%i\t%i\t%.2f\t%.2f\t%.2f" % (species, length, counts, meanDepth, coverage, cV)
            
