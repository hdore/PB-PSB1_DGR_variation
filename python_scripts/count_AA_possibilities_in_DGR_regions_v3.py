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

__progName__ = "count_AA_possibilities_in_DGR_regions"
__author__ = "Hugo Doré"
__copyright__ = "Copyright 2021, Hugo Doré, Wilbanks Lab, UC Santa Barbara"
__credits__ = ["Hugo Doré"]
__license__ = "CC-BY-NC-SA"
__status__ = "Dev"
__version__ = "1.0"
__last_modification__ = "2023-04-12"

#####
## Let's say I have the genome sequence, position of target, position of VR and position of TR

def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

def get_genome(genome_file) :
    genome_seqrecord = SeqIO.read(genome_file, "fasta")
    genome_seq = genome_seqrecord.seq
    
    return genome_seq


# Parse the input table 
def parse_regions(input_file, genome_seq, TR_from_file) :
    logging.info("Getting regions coordinates...")
    regions = {}
    target_seq={}
    A_pos_in_target={}
    TR_info={}
    with open(input_file, 'r') as  input_fh :
        for line in input_fh:
            if line[0] == '#':
                continue
            DGR, target, VRnum, VRid, target_start, target_end, VR_start, VR_end, VRsize, TR_start, TR_end, TR_seq_from_file = line.strip().split('\t')
            regions[VRid]={"DGR":DGR, "target":target, "VRnum":VRnum, "target_start":target_start, "target_end":target_end, "VR_start":VR_start, "VR_end":VR_end}
            # Get sequence of gene 
            
            if (target_start == "NA") or (target_end == "NA") or (VR_start == "NA") or (VR_end == "NA") or (TR_start == "NA") or (TR_end == "NA"):
                continue 
            
            if (int(target_end) - int(target_start)) > 0 : # This is if the gene is on + strand. Otherwise, take reverse complement
                target_seq[VRid] = genome_seq[int(target_start)-1:int(target_end)] #check if that should be end-1? 
            else :
                target_seq[VRid] = genome_seq[int(target_end)-1:int(target_start)].reverse_complement()
            
            if (int(TR_end) > int(TR_start)) : 
                TR_seq_from_coord = genome_seq[int(TR_start)-1:int(TR_end)] #check if that should be end-1? OK
            else : 
                TR_seq_from_coord = genome_seq[int(TR_end)-1:int(TR_start)].reverse_complement() #check if that should be end-1? OK
            
            if not TR_seq_from_coord == TR_seq_from_file :
                logging.warning("TR sequence you gave in input file does not match TR sequence extracted for coordinates. This might be on purpose if their is a mismatch between TR and VR you want to accomodate for, but we wanted to let you know. We are looking at VR %s, the TR sequence from input file is %s and the sequence from coordinates is %s.\n"%(VRid, TR_seq_from_file, TR_seq_from_coord))
            
            if TR_from_file :
                logging.info("Extracting info for VR: %s. You chose to continue with the TR sequence you gave in the file, and this is what we will do.\n"%(VRid))
                TR_seq = TR_seq_from_file
            else :
                logging.info("Extracting info for VR: %s. You chose to continue with the TR sequence extracted from TR coordinates, and this is what we will do.\n"%(VRid))
                TR_seq = TR_seq_from_coord
            A_pos_in_TR = findOccurrences(str(TR_seq), "A") #0-based index
            
            TR_info[VRid] = {}
            TR_info[VRid]["seq"] = TR_seq
            A_pos_in_target[VRid] = []
            for pos in A_pos_in_TR :
                A_pos_target = abs(int(VR_start) - int(target_start)) + pos 
                A_pos_in_target[VRid].append(A_pos_target) # I think that works with 0-based index here ? Should work on both strands, because we took reverse complement for -strand, so all that counts is distance from gene start # if on + strand 
                TR_info[VRid][A_pos_target] = pos
            logging.debug("A_pos_in_target[VRid] for VRid = %s is %s"%(VRid, A_pos_in_target[VRid]))
            
    
    return regions, target_seq, A_pos_in_target, TR_info
# Get sequence of gene 
# genes_seq={}
# genes_seq["DGR1_target1_VR1"] = genome[start-1:end] #check if that should be end-1? 
                                                # # This is if the gene is on + strand. Otherwise, take reverse complement




# Get T coordinates in TR seq
# TR_seq = genome[start-1:end] #check if that should be end-1? 
# A_pos_in_TR = {}
# A_pos_in_TR["DGR1_target1_VR1"]=findOccurrences(T_seq, "T") #0-based index

# # Get T coordinates in VR seq, then in gene seq 
# A_pos_in_VR = {}
# A_pos_in_VR["DGR1_target1_VR1"] = A_pos_in_TR["DGR1_target1"] + VR_start # If on + strand 
# A_pos_in_target = {}
# A_pos_in_target["DGR1_target1_VR1"] = A_pos_in_TR["DGR1_target1"] + abs(VR_start - target_start) # Should work on both strands, because we took reverse complement for -strand, so all that counts is distance from gene start # if on + strand 


# Then do the translation work 
#translation_table = dico{} 
def counA_possible_AA(target_seq, A_pos_in_target, TR_info) :
    logging.info("Analyzing data... \n")
    result = {}
    for VR in target_seq: #actually for each VR in a target
        logging.debug("Current VR is %s"%(VR))
        processed_pos = []
        seq = target_seq[VR]
        codon_start_list = []
        codonID = 0
        nb_possible_combinations = 1
        TR_seq = TR_info[VR]["seq"]
        for A_pos in A_pos_in_target[VR] :
            logging.debug("A_pos is %s (type: %s)"%(A_pos, type(A_pos)))
            if A_pos in processed_pos : 
                logging.debug("Already processed. Skipping.")
                continue 
            else : 
                codon_name = VR+"_codon"+str(codonID)
                codon_list=[]
                if A_pos%3 == 0 : 
                    place_in_codon = 1
                    #ref_codon = seq[A_pos:A_pos+3] 
                    TR_pos = TR_info[VR][A_pos]
                    if TR_pos+2 <= len(TR_seq)-1 : # if the codon does not end outside of TR (len(TR_seq)-1 is the index of the last letter in string)
                        ref_codon = TR_seq[TR_pos:TR_pos+3]
                    elif TR_pos+2 == len(TR_seq)-1+1 : #if there is one letter missing for the codon 
                        ref_codon = TR_seq[TR_pos:TR_pos+2] + seq[A_pos+2]
                    elif TR_pos+2 == len(TR_seq)-1+2 : #if there is 2 letters missing for the codon 
                        ref_codon = TR_seq[TR_pos] + seq[A_pos+1:A_pos+3]
                    
                    logging.debug("ref_codon is %s"%(ref_codon))
                    if (A_pos + 1 not in A_pos_in_target[VR]) and (A_pos + 2 not in A_pos_in_target[VR]) :
                        for nt in ["A", "C", "T", "G"] :
                            codon_list.append(Seq(nt + ref_codon[1] + ref_codon[2]))
                        processed_pos.append(A_pos)
                        var_pos = str(A_pos + 1) #+ 1 to get 1-based position (because that's what I will write)
                    elif (A_pos + 1 in A_pos_in_target[VR]) and (A_pos + 2 not in A_pos_in_target[VR]) :
                        place_in_codon = str(1)+";"+str(2)
                        for nt in ["A", "C", "T", "G"] :
                            for nt2 in ["A", "C", "T", "G"] :
                                codon_list.append(Seq(nt + nt2 + ref_codon[2]))
                        processed_pos.extend([A_pos, A_pos + 1])
                        var_pos = str(A_pos + 1)+";"+str(A_pos+1 + 1) #+ 1 to get 1-based position (because that's what I will write)
                    elif (A_pos + 1 not in A_pos_in_target[VR]) and (A_pos + 2 in A_pos_in_target[VR]) :
                        place_in_codon = str(1)+";"+str(3)
                        for nt in ["A", "C", "T", "G"] :
                            for nt2 in ["A", "C", "T", "G"] :
                                codon_list.append(Seq(nt + ref_codon[1] + nt2 ))
                        processed_pos.extend([A_pos, A_pos + 2])
                        var_pos = str(A_pos + 1)+";"+str(A_pos+2 + 1) #+ 1 to get 1-based position (because that's what I will write)
                    elif (A_pos + 1 in A_pos_in_target[VR]) and (A_pos + 2 in A_pos_in_target[VR]) :
                        logging.warning("The codon is made of 3 As in the TR. This is highly unlikely. VR is %s and first A pos is %s"%(VR, str(A_pos)))
                        place_in_codon = str(1)+";"+str(2)+";"+str(3)
                        for nt in ["A", "C", "T", "G"] :
                            for nt2 in ["A", "C", "T", "G"] :
                                for nt3  in ["A", "C", "T", "G"] :
                                    codon_list.append(Seq(nt + nt2 + nt3 ))
                        processed_pos.extend([A_pos, A_pos + 1, A_pos + 2])
                        var_pos = str(A_pos + 1)+";"+str(A_pos+1 + 1)+";"+str(A_pos+2 + 1) #+ 1 to get 1-based position (because that's what I will write)
                    AA_list = [str(codon.translate(table=11)) for codon in codon_list]
                    codon_start_pos = A_pos 
                elif A_pos%3 == 1 :
                    place_in_codon = 2
                    #ref_codon = seq[A_pos-1:A_pos+2]
                    TR_pos = TR_info[VR][A_pos]
                    if TR_pos - 1 < 0 : # if the first letter of codon is before the start of TR 
                        ref_codon = seq[A_pos-1] + TR_seq[TR_pos:TR_pos+2]
                    elif TR_pos+1 <= len(TR_seq)-1 : # if the codon does not end outside of TR (len(TR_seq)-1 is the index of the last letter in string)
                        ref_codon = TR_seq[TR_pos-1:TR_pos+2]
                    elif TR_pos+1 == len(TR_seq)-1+1 : #if there is one letter missing for the codon 
                        ref_codon = TR_seq[TR_pos-1:TR_pos+1] + seq[A_pos+2]
                    logging.debug("ref_codon is %s"%(ref_codon))
                    
                    codon_list = []
                    if A_pos - 1 in A_pos_in_target[VR] :
                        logging.error("This should not be possible. This position is second in codon and the first position in codon is an A. This codon should already have been processed. VR is %s and current position is %s"%(VR, str(A_pos)))
                        sys.exit(-1)
                    elif (A_pos + 1 not in A_pos_in_target[VR]) :
                        for nt in ["A", "C", "T", "G"] :
                            codon_list.append(Seq(ref_codon[0] + nt + ref_codon[2]))
                        processed_pos.append(A_pos)
                        var_pos = str(A_pos + 1) #+ 1 to get 1-based position (because that's what I will write)
                    elif (A_pos + 1 in A_pos_in_target[VR]) :
                        place_in_codon = str(2)+";"+str(3)
                        for nt in ["A", "C", "T", "G"] :
                            for nt2 in ["A", "C", "T", "G"] :
                                codon_list.append(Seq(ref_codon[0] + nt + nt2))
                        processed_pos.extend([A_pos, A_pos + 1])
                        var_pos = str(A_pos + 1)+";"+str(A_pos+1 + 1) #+ 1 to get 1-based position (because that's what I will write)
                    AA_list = [str(codon.translate(table=11)) for codon in codon_list]
                    codon_start_pos = A_pos -1 
                    
                elif A_pos%3 == 2 :
                    place_in_codon = 3
                    #ref_codon = seq[A_pos-2:A_pos+1]
                    TR_pos = TR_info[VR][A_pos]
                    if TR_pos - 2 >= 0 : # if the codon starts within the TR (it has to end within TR since the 3rd position is an A is the TR)
                        ref_codon = TR_seq[TR_pos-2:TR_pos+1]
                    elif TR_pos - 2 == -1 : # if the first letter of codon is before the start of TR 
                        ref_codon = seq[A_pos-2] + TR_seq[TR_pos-1:TR_pos+1]
                    elif TR_pos - 2 == -2 : #if the 2 first letters of codon is before the start of TR 
                        ref_codon = seq[A_pos-2:A_pos] + TR_seq[TR_pos]
                    logging.debug("ref_codon is %s"%(ref_codon))
                    
                    codon_list = []
                    if (A_pos - 2 in A_pos_in_target[VR]) or (A_pos - 1 in A_pos_in_target[VR]) :
                        logging.error("This should not be possible. This position is third in codon and the first or second position in codon is an A. This codon should already have been processed. VR is %s and current position is %s"%(VR, str(A_pos)))
                        sys.exit(-1)
                    else :
                        for nt in ["A", "C", "T", "G"] :
                            codon_list.append(Seq(ref_codon[0] + ref_codon[1] + nt)) #using Biopython to translate
                        processed_pos.append(A_pos)
                        var_pos = str(A_pos + 1) #+ 1 to get 1-based position (because that's what I will write)
                    AA_list = [str(codon.translate(table=11)) for codon in codon_list]
                    codon_start_pos = A_pos - 2
                logging.debug("place in codon was is %s"%(str(place_in_codon)))
                logging.debug("processed_pos is %s"%(processed_pos))
                nb_possible_stop = AA_list.count("*")
                nb_possible_AA = len(set(AA_list))
                codon_start_list.append(codon_start_pos)
                codonID+=1
                
                result[codon_name] = [VR, var_pos, codon_start_pos+1, ref_codon, place_in_codon, nb_possible_AA, nb_possible_stop] # I add 1 to positions to transform them to 1-based
                nb_possible_combinations *= nb_possible_AA
        if len(codon_start_list) > len(set(codon_start_list)) :
            #logging.warning("/!\/!\ 2 variable positions are in the same codon! This is not yet implemented and each are treated separately (no combination made)\nVR is: %s\npositions number is %s, and codon number is %s"%(VR, len(codon_start_list), len(set(codon_start_list))))
            logging.error("/!\/!\ 2 variable positions are in the same codon were not considered together! This should not happen any more.\nVR is: %s\npositions number is %s, and codon number is %s"%(VR, len(codon_start_list), len(set(codon_start_list))))
            sys.exit(-1)
        
        for codon in result :
            if result[codon][0] == VR :
                result[codon].append(nb_possible_combinations) #this will give the same value for all codons in a VR (this is the VR number of possible combinations)
        
    return result
    
def write_output(results, regions, output_file) :
    logging.info("Writing the results...\n")
    out_fh = open(output_file, 'w')
    out_fh.write("#The codon sequence is based on TR sequence (but based on where codon starts in the target). This is what makes sense since DGR replace the whole VR with the TR sequence, including non-variable positions.\n")
    out_fh.write("CodonID\tDGR\tTargetID\tVRnumber\tVRid\tVR_start\tVR_end\tA_pos_in_target_1based\tcodon_start_1based\tref_codon\tplace_in_codon\tnb_possible_AA\tnb_possible_stop\tVR_nb_possible_combinations\n")
    for codon in results : #mean for a position is calculated for positions with enough coverage only
        VRid = results[codon][0]
        out_fh.write(codon+"\t"+regions[VRid]["DGR"]+"\t"+regions[VRid]["target"]+"\t"+regions[VRid]["VRnum"]+"\t"+VRid+"\t"+regions[VRid]["VR_start"]+"\t"+regions[VRid]["VR_end"]+"\t"+"\t".join([str(i) for i in results[codon][1:]])+"\n")
    out_fh.close()
    


if __name__ == '__main__':
    
    ##OPTION PARSER##
    parser = argparse.ArgumentParser(description='Takes a fasta file with a single sequence and a table of DGR coordinates with target, VR and TR coordinates and TR sequence, and outputs the number of possible amino-acids for each putatively variable position in each VR of each target. V1 implements a way to take into account multiple variable positions within a codon. v2 considers codon sequence from the TR, not the target (but still using target sequence to know where codons start and end). v3 adds the possibility to use TR sequence from the input file to accomodate potential mismatches (in our case a 3bp deletion in the VR), and adds codon sequence to output file.',
        epilog='Last modification on '+__last_modification__+' by '+__author__)
    parser.add_argument("-r","--regions_file",dest = "regions_file", help = "Name or path to input file with the VR regions coordinates. Must be 'DGR	Target_ID	VR_number	VR_ID	Target start	Target end	VR start	VR end	VRsize	TR start	TR end	TR seq'.")
    parser.add_argument("-o","--output_file",dest = "output_file",help = "Name or path for output tab-delimited file.")
    parser.add_argument("-i --input_file", dest = "input_file", help = "Name or path of the fasta file of the genome (for now only single-sequence fasta files are accepted).")
    #parser.add_argument("-v","--verbose",dest="verbose",help="Enable verbose output.",action="store_true",default=False)
    parser.add_argument("--tr_from_file", dest="tr_from_file", action="store_true", default=False, help="Flag. If set, the program will use the TR sequence provided in the input regions file. Otherwise, it extracts the TR sequence from the fasta file based on the TR coordinates. By default the flag is not set and the sequence is extracted based on TR coordinates.")
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
        logging.error("Missing argument -r --regions_file (Name or path to input file with the VR regions coordinates. Must be 'DGR	Target_ID	VR_number	VR_ID	Target start	Target end	VR start	VR end	VRsize	TR start	TR end'.. Exiting.")
        sys.exit(-1)
    if args.output_file is None :
        logging.error("Missing argument -o --output_file (Name or path for output tab-delimited file). Exiting.")
        sys.exit(-1)
    if args.input_file is None :
        logging.error("Missing argument -i --input_file (Name or path of the fasta file of the genome (for now only single-sequence fasta files are accepted). Exiting.")
        sys.exit(-1)
    
    
    # Launch!
    genome_seq = get_genome(args.input_file)
    regions, target_seq, A_pos_in_target, TR_info = parse_regions(args.regions_file, genome_seq, args.tr_from_file)
    result = counA_possible_AA(target_seq, A_pos_in_target, TR_info)
    write_output(result, regions, args.output_file)
    
    



