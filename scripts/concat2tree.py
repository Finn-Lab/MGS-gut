#!/usr/bin/env python2

import sys
import os
import argparse
import glob
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def build_tree(tree):
        cmd_rax = ["raxmlHPC-PTHREADS-SSE3",
                   "-T", "32",
                   "-s", tree,
                   "-m", "PROTGAMMAAUTO",
                   "-p", "12345",
                   "-n", args.output_name+"_phylogeny"]
        cmd_fast = ["FastTree", tree]
        if args.tree_method == "raxml":
            subprocess.check_call(cmd_rax)
        elif args.tree_method == "fasttree":
            with open("FastTree_"+args.output_name+"_phylogeny.nwk", "w") as out:
                subprocess.check_call(cmd_fast, stdout=out)
        else:
            print "ERROR! -m needs to be 'raxml' or 'fasttree'"
            sys.exit(1)

def concat(args):
	seqs = {}
	for concat in glob.glob(os.path.join(args.directory, "*.fa*")):
		fasta_in = open(concat, "r")
		for record in SeqIO.parse(fasta_in, "fasta"):
			if record.id not in seqs.keys():
				seqs[record.id] = record.seq
			else:
				seqs[record.id] = seqs[record.id]+record.seq

        new_align = args.output_name+"_align.faa"
        if not os.path.isfile(new_align):	
	    fasta_out = open(new_align, "w")
	    for element in seqs.keys():
                newRecord = SeqRecord(seqs[element], id=element, description="")
		SeqIO.write(newRecord, fasta_out, "fasta") 
	tree_input = new_align
	return tree_input
	

if __name__ == '__main__':
        parser = argparse.ArgumentParser(usage='Concate alignments and build tree')
        parser.add_argument('-d', dest='directory', \
                                  help='Directory with protein alignments')
	parser.add_argument('-o', dest='output_name', \
                                  help='Name for the output files')
        parser.add_argument('-m', dest='tree_method', \
                                  help='Tree building method ("raxml" or "fasttree"); if raxml use 32 threads')
        if len(sys.argv) == 1:
                parser.print_help()
                sys.exit(1)
        else:
                args = parser.parse_args()
                tree_input = concat(args)
		build_tree(tree_input)
