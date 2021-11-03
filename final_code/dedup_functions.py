#!/usr/bin/env python

import argparse
import re

def get_args():                                         #function to parse command line arguments
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument("-s", "--sam", help="input SAM file to deduplex", required=True)
    parser.add_argument("-p", "--paired_end", help="designates if file is paired end", required=False)
    parser.add_argument("-u", "--umi", help="designates file containing the list of UMIs", required=False)
    return parser.parse_args() 


def store_umi_list(umi_file_path):                      # function to parse input file of known UMIs (Unique Molecular Identifiers)
    umi_list = []                                       # (if provided) and store them to a list
    for umi in umi_file_path:
        umi_list.append(umi.strip('\n'))
    return(umi_list)

def get_umi(sam_entry):                                 # function to retrieve the UMI from the QNAME a SAM file entry
    umi = sam_entry.split('\t')[0].split(':')[-1]
    return umi 

def get_position(sam_entry):                            # function to retrieve the start position from a SAM file entry
    pos = sam_entry.split('\t')[3]
    return pos

def get_rname(sam_entry):                               # function to retrieve the RNAME (chromosome) from a SAM file entry
    rname = sam_entry.split('\t')[2]
    return rname

def get_strand(sam_entry):                              # function to retireve the strand (forward (+) or reverse (-)) from a 
    flag = int(sam_entry.split('\t')[1])                # SAM file entry
    if ((flag & 16) == 16):
        return ('-')
    else:
        return ('+')  
    
def is_mapped(sam_entry):                               # function to determine if SAM file read is mapped using the bitwise flag 
    flag = int(sam_entry.split('\t')[1])                # returns 'True' is mapped, 'False' if unmapped
    if ((flag & 4) != 4):
        return True
    else:
        return False

def softclip_adjust(sam_entry):                         # function to calculate the true start position of a SAM file read based on
    cigar = sam_entry.split('\t')[5]                    # CIGAR string sequence. If the strand is a forward strand, starting position              
    pos = int(get_position(sam_entry))                  # is calculated by subtracting the number of softclipped (S) bases on the left
    adjusted_pos = pos                                  # of the CIGAR string. If the strand is a reverse strand, starting position is 
    right_softclipped = 0                               # calculated by summing all numeric values in the CIGAR string and adding this
    inserted = 0                                        # this to the starting position reported in the SAM file, and then subtracting 
                                                        # any bases that were softclipped on the right side (S) or inserted (I)
    if get_strand(sam_entry) == '+':                    # if the strand is forward
        if 'S' in cigar:                                # if it is softclipped 
            clip_num = cigar.split('S')[0]              # find the num bases softclipped     
            if clip_num.isdigit():                      
                adjusted_pos = pos - int(clip_num)      # subtract softclipped bases from position
                return adjusted_pos                     # return adjusted position
            else:
                return pos                              # else return original position
        else:
            return pos
        
    elif get_strand(sam_entry) == '-':                  # else if the strand is reverse
        temp = re.split('([A-Z])', cigar.strip("\n"))   # split cigar string into list of elements
        for i in temp:                                  # for every element in the list of CIGAR                   
            if i.isdigit():                             # is element is a number                            
                adjusted_pos += int(i)                  # add that value to the adjusted position           
        if 'I' in temp:                                 # if I is in the temp list
            indices = [index for index, element in enumerate(temp) if element == "I"] # enumerate list indicies
            for index in indices:                       # find value of insterted in CIGAR
                inserted += int(temp[index-1])
            
        adjusted_pos = adjusted_pos - inserted          # subtract that value from adjusted pos
            
        return adjusted_pos
        
        
def check_duplicate(sam_entry, dict, umi_list = None):  # this function checks if the SAM entry is a duplicate, by comparing it to 
    umi = get_umi(sam_entry)                            # a nested dictionary of visited UMIs/Chromosome/Strand/Start combinations 
    chr = get_rname(sam_entry)                  
    strand = get_strand(sam_entry)                      # nested dictionary anatomy: 
    start = softclip_adjust(sam_entry)                  #{UMI1:{'chrom + strand': [start1, start2], 'chrom + strand': [start]}, UMI2: {'chrom + strand': [start, start]}}

    key = chr + strand                                  # dictionary keys are a string of chromosome number plus strand
   
    if umi_list != None and umi in umi_list or umi_list == None: # if known UMI file was provided vs randomers    
        if umi in dict:                                          # if the UMI is in the dictionary
            if key in dict[umi]:                                 # if the key is in the respective UMI dictionary
                if start in dict[umi][key]:                      # if the start position is in the list of the key values
                    return "duplicate"                           # we have found a duplicate read  
                else:                                               
                    dict[umi][key].append(start)                 # else add the unvistied start position to the list   
                    return "good umi"
            else:
                dict[umi][key] = [start]                         # else create a new inner dictionary key value pair 
                return "good umi"
        else: 
            dict[umi] = {}                                       # else create a new outer dictionary  
            dict[umi][key] = [start]                             # create a new inner dictionary key value pair 
            return "good umi"
    else:
        return "bad umi"                                         # UMI was not found in the known UMI list
   

def new_dict_entry(line, dict):                                  # function to create a new dictionary entry
    key = get_umi(line)
    key2 = get_rname(line) + get_strand(line)
    dict[key] = {}
    dict[key][key2] = [softclip_adjust(line)]
    return dict

