#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import datetime
import logging
import sys
import os
import argparse
#from sys import stdout
from Bio.Seq import Seq
#from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
#from Bio.SeqFeature import SeqFeature, FeatureLocation

#import cPickle
#import random
import numpy as np
#import matplotlib.pyplot as plt
#from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord


__progName__ = "nt_prop_from_pileup"
__author__ = "Hugo Doré"
__copyright__ = "Copyright 2021, Hugo Doré, Wilbanks Lab, UC Santa Barbara"
__credits__ = ["Hugo Doré"]
__license__ = "CC-BY-NC-SA"
__status__ = "Dev"
__version__ = "2.0"
__last_modification__ = "2023-02-17"

def parse_regions(input_file) :
    logging.info("Getting regions coordinates...")
    regions = {}
    with open(input_file, 'r') as  input_fh :
        for line in input_fh:
            if line[0] == '#': #The header line must start with # if there is one.
                continue
            DGR, VRnum, VRid, target_start, target_end, VR_start, VR_end, TR_start, TR_end  = line.strip().split('\t')
            regions[VRid]={"DGR":DGR, "VRnum":VRnum, "target_start":target_start, "target_end":target_end, "VR_start":VR_start, "VR_end":VR_end, "TR_start":TR_start, "TR_end":TR_end}
    
    return regions

def parse_mpileup(input_file, regions, cov_cutoff) : #I don't take into account insertions and deletions here...
    logging.info("Harnessing nucleotide counts and calculating diversity metrics...\n")
    results={}
    
    with open(input_file, 'r') as  input_fh :
        lines=input_fh.readlines() 
        
        for region in regions :
            logging.info("Working with region %s \n"%(region))
            ## Get TR sequence 
            TR_start = int(regions[region]["TR_start"])
            TR_end = int(regions[region]["TR_end"])
            if TR_start > TR_end : #if the TR is on the minus strand
                logging.debug("TR is on the minus strand")
                TR_new_start = TR_end - 1 #this is because we change to a 0-based index in python (There is a 1 difference between the python list index and the genome position
                TR_new_end = TR_start - 1 + 1 # Change of index AND (in range(x,y), y is excluded)
                TR_strand = - 1 
            else : #if the VR is on the + strand 
                TR_new_start = TR_start - 1
                TR_new_end = TR_end
                TR_strand = 1
            TR_seq = ""
            for i in range(TR_new_start, TR_new_end, 1) : # This will extract the sequence that is on the + strand
                TR_seq += lines[i].strip().split('\t')[2]
            if TR_strand == -1 : # Then we reverse complement if TR was on the - strand
                TR_seq = str(Seq(TR_seq).reverse_complement())
            
            logging.info("Extracted TR sequence is %s \n"%(TR_seq))
            
            
            start = int(regions[region]["VR_start"])
            end = int(regions[region]["VR_end"])
            ## Determine VR start and end depending on strand 
            if start > end : #if the VR is on the minus strand
                logging.debug("VR is on the minus strand")
                new_start = end - 1 #this is because we change to a 0-based index in python (There is a 1 difference between the python list index and the genome position
                new_end = start - 1 + 1 # Change of index AND (in range(x,y), y is excluded)
                strand = - 1 
            else : #if the VR is on the + strand 
                new_start = start - 1
                new_end = end
                strand = 1
            
            TR_index = 0
            results[region] = {}
            #results[region]["pi"]=[]
            #results[region]["ratio"]=[]
            #results[region]["n_or_del"]=[]
            for i in range(new_start,new_end, 1) : # /!\ Going on + strand whatever the VR strand is. Coordinates are always + strand.
                scaff, pos, ref, cov, reads, qual = lines[i].strip().split('\t')
                cov = int(cov)
                logging.debug("Coverage from the mpileup file at this position (position %s) is %s"%(pos, cov))
                results[region][pos] = {}
                if cov >= cov_cutoff: #only consider positions for which coverage is >= cutoff (in the future, maybe put a condition on the mean coverage of the region of interest?
                    ref=ref.lower()
                    logging.debug("Reference base is %s"%(ref))
                    reads=reads.lower()
                        #remove the insertion and deletion 
                    logging.debug("mapped reads are %s"%(reads))
                    while '^' in reads :
                        logging.debug("Reads with read start is %s"%(reads))
                        position = reads.index("^")
                        new_reads = reads[:position]+reads[position+1+1:]
                        reads = new_reads 
                        logging.debug("Read without read start is %s"%(reads))
                    while "-" in reads: 
                        position = reads.index("-")
                        logging.debug("Number of deleted bases is %s"%(reads[position+1]))
                        length = int(reads[position+1])
                        new_reads = reads[:position]+reads[position+length+1+1:]
                        reads = new_reads 
                        logging.debug("mapped reads are %s"%(reads))
                    while "+" in reads: 
                        position = reads.index("+")
                        length = int(reads[position+1])
                        new_reads = reads[:position]+reads[position+length+1+1:]
                        reads = new_reads 
                    reads=reads.replace(".", ref).replace(",", ref)
                    logging.debug("mapped reads after replacing to ref base are %s"%(reads))
                    ref_count = reads.count(ref)
                    a_count = reads.count("a")
                    t_count = reads.count("t")
                    g_count = reads.count("g")
                    c_count = reads.count("c")
                    n_or_del_count = reads.count("n") + reads.count("*")
                    logging.debug("a_count is %s and t_count is %s and g_count is %s and c_count is %s"%(a_count, t_count, g_count, c_count))
                    
                    bases_sum = a_count + t_count+g_count+c_count
                    if bases_sum != cov :
                        logging.info("Warning: the base count does not match coverage value!")
                        #print("modified reads is", reads)
                        logging.info("a_count is %s and t_count is %s and g_count is %s and c_count is %s"%(a_count, t_count, g_count, c_count))
                        logging.info("n + deletions count is %s"%(n_or_del_count))
                        logging.info("the modified reads are %s and cov is %s and sum of bases is %s \n"%(reads, str(cov), str(bases_sum)))
                    #results[region][pos]["pi"] = 1-((a_count/cov)**2+(t_count/cov)**2+(g_count/cov)**2+(c_count/cov)**2) #This cov value counts n, which is not clean when calculating the pi... I guess this could be discussed.
                    #results[region][pos]["ratio"] = (cov-ref_count)/cov #This cov value counts n, which is not clean when calculating the ratio
                    results[region][pos]["pi"] = 1-((a_count/bases_sum)**2+(t_count/bases_sum)**2+(g_count/bases_sum)**2+(c_count/bases_sum)**2)
                    results[region][pos]["ratio"] = (bases_sum-ref_count)/bases_sum
                    results[region][pos]["n_or_del"] = n_or_del_count
                    
                    results[region][pos]["ref"] = ref
                    if TR_index < len(TR_seq) :
                        if strand == 1 :
                            results[region][pos]["TR_ref"] = TR_seq[TR_index]
                        elif strand == -1 :
                            results[region][pos]["TR_ref"] = TR_seq[-TR_index-1] #because we are walking through the VR on the + strand. -1 because index is 0-based but reverse index is -1-based!
                    else : 
                        results[region][pos]["TR_ref"] = "X"
                    results[region][pos]["a_ct"] = a_count 
                    results[region][pos]["t_ct"] = t_count 
                    results[region][pos]["c_ct"] = c_count 
                    results[region][pos]["g_ct"] = g_count 
                    
                    
                else :
                    results[region][pos]["pi"] = "NA" #np.nan
                    results[region][pos]["ratio"] = "NA"
                    results[region][pos]["n_or_del"] = "NA"
                    results[region][pos]["ref"] = ref
                    if TR_index < len(TR_seq) :
                        if strand == 1 :
                            results[region][pos]["TR_ref"] = TR_seq[TR_index]
                        elif strand == -1 :
                            results[region][pos]["TR_ref"] = TR_seq[-TR_index] #because we are walking through the VR on the + strand
                    else : 
                        results[region][pos]["TR_ref"] = "X"
                    results[region][pos]["a_ct"] = "NA"
                    results[region][pos]["t_ct"] = "NA"
                    results[region][pos]["c_ct"] = "NA"
                    results[region][pos]["g_ct"] = "NA"
                
                results[region][pos]["VR_strand"] = str(strand)
                TR_index += 1
                
        
    return results
    

def write_output(results, regions, output_file) :
    logging.info("Writing the results...\n")
    out_fh = open(output_file, 'w')
    out_fh.write("VRid\ttarget_id\tVR_start\tVR_end\tVR_strand\tposition\tref_allele\tTR_nt\tnt_value\tnt_count\tpi\tnonref_ratio\n")
    for region in results : #mean for a position is calculated for positions with enough coverage only
        for pos in results[region] :
            for nt in ["n_or_del", "a_ct", "t_ct", "g_ct", "c_ct"] : 
                out_fh.write(region+"\t"+regions[region]["DGR"]+"\t"+regions[region]["VR_start"]+"\t"+regions[region]["VR_end"]+"\t"+results[region][pos]["VR_strand"]+"\t"+str(pos)+"\t"+results[region][pos]["ref"].upper()+"\t"+results[region][pos]["TR_ref"]+"\t"+nt+"\t"+str(results[region][pos][nt])+"\t"+str(results[region][pos]["pi"])+"\t"+str(results[region][pos]["ratio"])+"\n")
    out_fh.close()

if __name__ == '__main__':
    
    ##OPTION PARSER##
    parser = argparse.ArgumentParser(description='Takes a table file with the coordinates of the regions of interest, and a pileup file obtained by running samtools mpileup on a bam file (with -a and --fasta-ref options). Outputs a table of each nt proportion at each position of each region. It also calculates nucleotide diversity (pi) and non-reference allele proportion as a bonus. v1 modifies pi and non-reference allele proportion to use the sum of A + T + C + G instead of the coverage calculated by mpileup, which includes Ns. V2 brings a correction for when VR is on the minus strand.',
        epilog='Last modification on '+__last_modification__+' by '+__author__)
    parser.add_argument("-r","--regions_file",dest = "regions_file", help = "Name or path to input file with the VR and TR regions coordinates.")
    parser.add_argument("-o","--output_file",dest = "output_file",help = "Name or path for output tab-delimited file.")
    parser.add_argument("-i", "--input_file", dest = "input_file", help = "Name or path of input_file from samtools mpileup with -a and --fasta-ref options.")
    parser.add_argument("-c", "--cov_cutoff", dest = "cov_cutoff", default=5, help = "Int. The coverage cutoff at a position to calculate statistics. Default: 5.")
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
    
    logging.info("Congrats, it's all done!")
    

