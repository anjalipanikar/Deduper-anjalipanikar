#!/usr/bin/env python

########################################################################################################
# Program: De-duplexing Algorithm for Mapped and Aligned Illumina Sequencing Data 
#  
# Creator: Geethanjali Panikar 
#
# Description: An algorithm to filter potential PCR duplicated sequence reads from a trimmed and aligned SAM file.
# Algorithm additionally filters data for unwwanted 'unmapped' reads and if the option is chosen, reads 
# that do not contain a Unique Molecular Identifier (UMI) present in the provided UMI file.
# 
# Input: 1 required files: The SAM file to be deduplexed. 
#        1 optional file: A file of known UMIs
#
# Output: to output directory (deduped/), 4 files: (1) deduplexed data (.deduplexed.sam)
#                                                  (2) extracted duplicate reads (.duplicates.sam)
#                                                  (3) extracted unmapped reads (.unmapped.sam)
#                                                  (4) extracted bad reads that did not match known UMIs
#  
# Arguments to Run Program: $ python3 panikar_deduper.py -s <SAM file> -u <optional UMI file> 
#########################################################################################################

import os
import sys
import dedup_functions 

# -------------------------------------------- COLLECT ARGUMENTS AND PRE-PROCESS INPUT SAM FILE WITH UNIX  --------------------------------------------- # 

# collect arguments 

args = dedup_functions.get_args()                                   
sam_file = args.sam
umi_file = args.umi

if args.paired:
    print('Cannot process paired-end data. Run this program with single-end reads only.')
    sys.exit(1)

# declare file names 

sorted_sam = "sorted_temp.sam"
deduped_sam = "deduped.sam"
final_sam = sam_file.split("/")[-1].split(".sam")[0] + '.deduplexed.sam'
duplicates_sam = sam_file.split("/")[-1].split(".sam")[0] + '.duplicates.sam'
bad_sam = sam_file.split("/")[-1].split(".sam")[0] + '.bad_umi.sam'
unmapped_sam = sam_file.split("/")[-1].split(".sam")[0] + '.unmapped.sam'
bam_file = sam_file.split("/")[-1].split(".sam")[0] + ".bam"
sorted_bam = sam_file.split("/")[-1].split(".sam")[0] + ".sorted.bam"
sorted_no_headers = "sorted_no_headers.sam"

# execute uniux commands to use samtools to sort the same file, translocate headers to final output, 
# and remove headers from file to parse 

os.system("samtools view -S -b " + sam_file + " > " + bam_file)
os.system("samtools sort " + bam_file + " -o " + sorted_bam)
os.system("samtools view -h " + sorted_bam + " > " + sorted_sam)
os.system("grep '^[[:space:]]*@' " + sorted_sam + ' > ' + final_sam)
os.system("grep -v '^[[:space:]]*@' " + sorted_sam + " > " + sorted_no_headers)


# -------------------------------------------------------- OPEN FILES AND INITIALIZE VARIABLES ---------------------------------------------------------- # 


# set umi filehandle and umi list addresses to 'None' as default 

uf = None 
umi_list = None

# if an known umi file has been inputted, create a file handle for the file and 
# set a adress value for the umi file handle. In addition, create a list of known
# umis from the file and store in the umi list.

if umi_file != None:
    uf = open(umi_file, "r")
    umi_list = dedup_functions.store_umi_list(uf)

# open file handle connections 

sf = open(sorted_no_headers, "r")
sf2 = open(deduped_sam, "w")
dp = open(duplicates_sam, "w")
bad = open(bad_sam, "w")
unmapped = open(unmapped_sam, "w")

# create a dictionary to hold umis and corresponding chr #, strand, and start position

umi_dict = {}
file_lines = 1 
duplicates = 0

# ------------------------------------------------------------ ALGORITHM FOR DEDUPLEXING DATA ----------------------------------------------------------- # 

#read in first SAM entry of the file 

first_line = sf.readline()

# create a dictionary entry for the first SAM entry

dedup_functions.new_dict_entry(first_line, umi_dict)

# set the value of prev_chrom equal to the chrom of this SAM entry

prev_chrom = dedup_functions.get_rname(first_line) 

# write the first SAM entry to the output deduplexed file

sf2.write(first_line)

# for each line in the SAM file

for line in sf:
    
    current_chrom = dedup_functions.get_rname(line)                         # set current chromosome variable to current chrom
    if prev_chrom != current_chrom:                                         # if the prev chrom variable is not equal to the current chrom variable                                         
        umi_dict.clear()                                                    # clear the dictionary and start new
        prev_chrom = current_chrom                                          # set the previous chrom to current chrom                                       
    else: 
        prev_chrom = current_chrom                                          # if the prev chrom equals current chrom, set 
                                                                            # prev chrom to current chrom
    
    if dedup_functions.is_mapped(line) == True:                             # if the SAM entry is mapped, continue 
                                                                            
        flag = dedup_functions.check_duplicate(line, umi_dict, umi_list)    # set flag variable to output of check_duplicate()
        
        if flag == 'duplicate':                                             # if the SAM entry is a duplicate, write to duplicates file
            dp.write(line)
            duplicates += 1
        elif flag == 'good umi':                                            # if the SAM entry is a not a duplicate, write to deduplexed file 
            sf2.write(line)
        elif flag == 'bad umi':                                             # if the SAM entry is an unknown UMI, write to bad file
            bad.write(line)
    else:
        unmapped.write(line)                                                # if unmapped, write to unmapped output
        continue

    file_lines += 1

# close files 

sf.close()
sf2.close()
dp.close()
bad.close()
unmapped.close()

# concatenate deduped sam file to final sam file with headers
cmd = 'cat ' + deduped_sam + ' >> ' + final_sam
os.system(cmd)

# remove unwanted files 

os.system("rm -rf " + sorted_sam)
os.system("rm -rf " + deduped_sam)
os.system("rm -rf " + bam_file)
os.system("rm -rf " + sorted_bam)
os.system("rm -rf " + sorted_no_headers)
    
