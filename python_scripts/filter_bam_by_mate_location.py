#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import datetime
import logging
import sys
import os
import argparse
#from sys import stdout
from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
#from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

#import cPickle
#import random
#import numpy as np
#import matplotlib.pyplot as plt
#from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord
import pysam

__progName__ = "filter_bam_by_mate_location"
__author__ = "Hugo Doré"
__copyright__ = "Copyright 2021, Hugo Doré, Wilbanks Lab, UC Santa Barbara"
__credits__ = ["Hugo Doré"]
__license__ = "CC-BY-NC-SA"
__status__ = "Dev"
__version__ = "1.0"
__last_modification__ = "2022-12-13"

def bam_list_parser(bam_list_file) :
    logging.info("Parsing bam list file...\n\n")
    bam_list_fh = open(bam_list_file, "r") 
    bam_list = [x.strip() for x in bam_list_fh.readlines()] 
    bam_list_fh.close()
    
    return bam_list

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

# def get_region(genome_file, contig, my_start, my_end) :
    # genome_seqrecords = SeqIO.parse(genome_file, "fasta")
    # for rec in genome_seqrecords : 
        # if rec.id == contig :
            # genome_seq = rec.seq
            # region_seq = genome_seq[my_start:my_end]
    
    # return region_se

def bam_list_wrapper(bam_list, regions, mq_threshold, window, suffix) : #for loop on all bam in the list ##contigs (add that option later if we want to accept bam files with multiple contigs)
    for bam in bam_list :
        logging.info("Working on bam file %s"%(bam))
        get_subsetted_bam(regions, bam, mq_threshold, window, suffix) ##contigs

def get_subsetted_bam(regions, bam, mq_threshold, window, suffix) : #get subsetted bam for a given bam and all given regions 
    samfile = pysam.AlignmentFile(bam, "rb")
    # first prepare your output bam 
    output_bam = bam.replace(".bam", suffix+".bam").replace(".sorted", "")
    logging.info("\nCreating new header...\n")
    original_header = samfile.header #.to_dict()
    logging.debug("Original header:\n\n\n%s"%(original_header))
    new_header = {}
    for field in original_header :
        new_header[field] = original_header[field]
    
    logging.info("\Writing reads with sufficient mapQ whose mate is in the chosen window for each region of interest in output file:\n%s\n"%(output_bam))
    subsetted_bam = pysam.AlignmentFile(output_bam, "wb", header=new_header)
    for region in regions :
        VR_start = int(regions[region]["VR_start"])
        VR_end = int(regions[region]["VR_end"])
        my_start=min(VR_start, VR_end) #to ensure that "start" is in coordinates on the genome, not relative to gene orientation. This won't cause issues because we want to keep reads whose mate within X kb upstream OR downstram
        my_end=max(VR_start, VR_end)
        read_in_region_ct = 0
        read_kept_ct = 0
        for read in samfile.fetch('psb-scaff03', my_start, my_end):
            read_in_region_ct += 1
            if read.is_proper_pair : #is_paired:
                mate = samfile.mate(read) # A warning from pysam help that might need to be considered: Calling this method will change the file position. This might interfere with any iterators that have not re-opened the file
                logging.debug("Read is paired. Read is \n%s\nand mate id is %s\n"%(read.query_name, mate.query_name))
                if (not mate.is_unmapped) and (mate.reference_end + 1 >= my_start-window) and (mate.reference_start + 1 <= my_end+window) and (read.mapping_quality > mq_threshold) : #reference_start and reference_end are 0-based. So I added 1
                    logging.debug("Found a read in region of interest with sufficient mapQ AND a mate within the window!")
                    read_kept_ct += 1
                    subsetted_bam.write(read) # I just write the read that maps to the region of interest. But if you also want the mate, uncomment the following
                    # subsetted_bam.write(mate)
            else :
                logging.debug("Read %s has was registered as not paired in the bam file"%(read.query_name))
        logging.info("\n\nIn bam file %s, %s reads were kept out of %s mapping in region %s"%(bam, read_kept_ct, read_in_region_ct, region))
    


if __name__ == '__main__':
    
    ##OPTION PARSER##
    parser = argparse.ArgumentParser(description='Takes as input a list of bam files (1 per line), start and end coordinates of regions of interest and a window size, and creates new bam files with only reads mapping within the region of interest whose mate is within the defined window.',
        epilog='Last modification on '+__last_modification__+' by '+__author__)
    parser.add_argument("-b","--bam_list",dest = "bam_list", help = "Name or path to input list of bam files. Must be 1 file/line, with full path.")
    parser.add_argument("-r","--regions_file",dest = "regions_file", help = "Name or path to input file with the VR regions coordinates.")
    parser.add_argument("-s","--suffix", dest = "suffix", default="subsetted", help = "Suffix for output subsetted bam file. The '.bam' of the input bam will be replaced by 'suffix.bam'. Default is 'subsetted'")
    parser.add_argument("--mq_threshold",dest = "mq_threshold", default=0, help = "Int. Minimum mapping quality for a read to be considered. Default is 0.")
    parser.add_argument("--window",dest = "window", default=1000, help = "Int. Window in which the mate is looked for. The mate will be looked for in the sequence defined by [region_start - window, region_end + window]. Default is 1000.")
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
    
    
    
    if args.bam_list is None :
        logging.error("Missing argument -b --bam_list (Name or path to input list of bam files. Must be 1 file/line, with full path)... Exiting.")
        sys.exit(-1)
    if args.regions_file is None :
        logging.error("Missing argument -r --regions_file (Name or path to input file with the VR regions coordinates). Exiting.")
        sys.exit(-1)
        
    
    # Launch!
    
    regions = parse_regions(args.regions_file)
    bam_list = bam_list_parser(args.bam_list)
    bam_list_wrapper(bam_list, regions, int(args.mq_threshold), int(args.window), args.suffix) #for loop on all bam in the list ##contigs (add that option later if we want to accept bam files with multiple contigs)
    
    
    
    
