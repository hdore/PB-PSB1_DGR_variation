#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import datetime
import logging
import sys
import os
import argparse
#from sys import stdout
#from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
#from Bio.SeqFeature import SeqFeature, FeatureLocation

#import cPickle
#import random
import numpy as np
#import matplotlib.pyplot as plt
#from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord


__progName__ = "calc_pi_from_pileup"
__author__ = "Hugo Doré"
__copyright__ = "Copyright 2021, Hugo Doré, Wilbanks Lab, UC Santa Barbara"
__credits__ = ["Hugo Doré"]
__license__ = "CC-BY-NC-SA"
__status__ = "Dev"
__version__ = "1.0"
__last_modification__ = "2021-04-22"

def parse_regions(input_file) :
    logging.info("Getting regions coordinates...")
    regions = {}
    with open(input_file, 'r') as  input_fh :
        for line in input_fh:
            if line[0] == '#':
                continue
            DGR, VRnum, VRid, target_start, target_end, VR_start, VR_end  = line.strip().split('\t')
            regions[VRid]={"DGR":DGR, "VRnum":VRnum, "target_start":target_start, "target_end":target_end, "VR_start":VR_start, "VR_end":VR_end}
    
    return regions

def parse_mpileup(input_file, regions, cov_cutoff) : #I don't take into account insertions and deletions here...
    logging.info("Calculating diversity metrics...")
    results={}
    
    with open(input_file, 'r') as  input_fh :
        lines=input_fh.readlines() 
        
        for region in regions :
            start = int(regions[region]["VR_start"])
            end = int(regions[region]["VR_end"])
            results[region] = {}
            results[region]["pi"]=[]
            results[region]["ratio"]=[]
            results[region]["n_or_del"]=[]
            
            for i in range(start-1,end) : # will start at start-1 and stop at end-1 (in range(x,y), y is excluded). There is a 1 difference between the list index and the genome position
                scaff, pos, ref, cov, reads, qual = lines[i].strip().split('\t')
                cov = int(cov)
                logging.debug("Coverage from the mpileup file at this position (position %s) is %s"%(pos, cov))
                if cov >= cov_cutoff: #only consider positions for which coverage is >= cutoff (in the future, maybe put a condition on the mean coverage of the region of interest?
                    ref=ref.lower()
                    reads=reads.lower()
                        #remove the insertion and deletion 
                    logging.debug("mapped reads are %s"%(reads))
                    while '^' in reads :
                        logging.debug("Reads with read start is %s"%(reads))
                        pos = reads.index("^")
                        new_reads = reads[:pos]+reads[pos+1+1:]
                        reads = new_reads 
                        logging.debug("Read without read start is %s"%(reads))
                    while "-" in reads: 
                        pos = reads.index("-")
                        logging.debug("Number of deleted bases is %s"%(reads[pos+1]))
                        length = int(reads[pos+1])
                        new_reads = reads[:pos]+reads[pos+length+1+1:]
                        reads = new_reads 
                        logging.debug("mapped reads are %s"%(reads))
                    while "+" in reads: 
                        pos = reads.index("+")
                        length = int(reads[pos+1])
                        new_reads = reads[:pos]+reads[pos+length+1+1:]
                        reads = new_reads 
                    reads=reads.replace(".", ref).replace(",", ref)
                    ref_count = reads.count(ref)
                    a_count = reads.count("a")
                    t_count = reads.count("t")
                    g_count = reads.count("g")
                    c_count = reads.count("c")
                    n_or_del_count = reads.count("n") + reads.count("*")
                    bases_sum = a_count + t_count+g_count+c_count
                    if bases_sum != cov :
                        logging.info("Warning: the base count does not match coverage value!")
                        #print("modified reads is", reads)
                        logging.info("a_count is %s and t_count is %s and g_count is %s and c_count is %s"%(a_count, t_count, g_count, c_count))
                        logging.info("n + deletions count is %s"%(n_or_del_count))
                        logging.info("the modified reads are %s and cov is %s and sum of bases is %s \n"%(reads, str(cov), str(bases_sum)))
                    results[region]["pi"].append(1-((a_count/cov)**2+(t_count/cov)**2+(g_count/cov)**2+(c_count/cov)**2))
                    results[region]["ratio"].append((cov-ref_count)/cov)
                    results[region]["n_or_del"].append(n_or_del_count)
                    
                else :
                    results[region]["pi"].append(-1)
                    results[region]["ratio"].append(-1)
                    results[region]["n_or_del"].append(-1)
                
    return results
    

def write_output(results, regions, output_file) :
    logging.info("Writing the results...\n")
    out_fh = open(output_file, 'w')
    out_fh.write("VRid\ttarget_id\tVR_start\tVR_end\tmean_pi\tmean_nonref_ratio\tpi_list\tnonref_ratio_list\tn_and_del_count_list\n")
    for region in results : #mean for a position is calculated for positions with enough coverage only
        if len([i for i in results[region]["pi"] if i >= 0 ]) > 0 : #if there is data for pi, there is some for ratio
            pi_mean = str(np.mean([i for i in results[region]["pi"] if i >= 0 ]))
            ratio_mean = str(np.mean([i for i in results[region]["ratio"] if i >=0]))
        else :
            pi_mean = "NA"
            ratio_mean = "NA"
        out_fh.write(region+"\t"+regions[region]["DGR"]+"\t"+regions[region]["VR_start"]+"\t"+regions[region]["VR_end"]+"\t"+pi_mean+"\t"+ratio_mean+"\t"+",".join([str(i) for i in results[region]["pi"]])+"\t"+",".join([str(i) for i in results[region]["ratio"]])+"\t"+",".join([str(i) for i in results[region]["n_or_del"]])+"\n")
    out_fh.close()

if __name__ == '__main__':
    
    ##OPTION PARSER##
    parser = argparse.ArgumentParser(description='This script [link to python script] allows to calculate the nucleotide diversity (pi) and the proportion of non-reference alleles from a bam file pileup, at positions defined by the user.',
        epilog='Last modification on '+__last_modification__+' by '+__author__)
    parser.add_argument("-r","--regions_file",dest = "regions_file", help = "Name or path to input file with the VR regions coordinates.")
    parser.add_argument("-o","--output_file",dest = "output_file",help = "Name or path for output tab-delimited file.")
    parser.add_argument("-i --input_file", dest = "input_file", help = "Name or path of input_file from samtools mpileup with -a and --fasta-ref options.")
    parser.add_argument("-c --cov_cutoff", dest = "cov_cutoff", default=5, help = "Int. The coverage cutoff at a position to calculate statistics.")
    #parser.add_argument("-v","--verbose",dest="verbose",help="Enable verbose output.",action="store_true",default=False)
    parser.add_argument("-d", "--debug_mode", dest="debug_mode", action="store_true", default=False, help="Enable debug mode (very verbose).")
    args = parser.parse_args()
    
    
    
    
    if args.debug_mode is True :
        logging.basicConfig(level=logging.DEBUG)
    else :
        logging.basicConfig(level=logging.INFO)
    
    logging.info(__progName__ + " " + __version__ + " by " + __author__)
    logging.info("Last modification on " + __last_modification__ + "\n")
    logging.info(time.strftime("%Y-%m-%d %X %Z") +"\n")
    
    
    
    if args.regions_file is None :
        logging.error("Missing argument -r --regions_file (Name or path to input file with the VR regions coordinates). Exiting.")
        sys.exit(-1)
    if args.output_file is None :
        logging.error("Missing argument -o --output_file (Name or path for output tab-delimited file). Exiting.")
        sys.exit(-1)
    if args.input_file is None :
        logging.error("Missing argument -i --input_file (Name or path of input_file from samtools mpileup with -a and --fasta-ref options). Exiting.")
        sys.exit(-1)
    
    
    # Launch!
    regions = parse_regions(args.regions_file)
    results = parse_mpileup(args.input_file, regions, int(args.cov_cutoff))
    write_output(results, regions, args.output_file)
    
    

