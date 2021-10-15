# Deduplexing Pseudocode  

### Anjali Panikar 

<br/>



### ***Define the problem*** 
After quality filtering and trimming reads produced in an RNA-seq experiment, we align these reads to a genome and alignment information is outputted in the form of SAM file. The SAM file in it's raw form however, will not be completely ready to use in downstream analyses due to the potential presence of **PCR Duplicates** within the file. PCR duplicates of RNA reads are an unwanted by-product of PCR library amplification, and if they are not removed from the SAM file they may be falsley over-represented during reference genome mapping. PCR duplicates can arise during library prep for various reasons but usually by chance, and for certain types of libraries and fragments this chance can increase. A common reason for fragment PCR duplication is that the fragment is an already highly expressed segment of RNA, so the chances of it being duplicated increased substantially compared to regularly expressed fragments. In addition, a low starting library concentration will inevitably require more PCR cycles, further increasing the chance that a mistake will be made during sequence amplification. Multiple things can be done to cleam these duplicates from the SAM files, and there are specific SAM file indicators that we can look at and compare to determine true duplicates versus highly expressed genes. 


### ***Pseudocode***

#### > Pre-Run File Sorting: 
Using unix commands, the SAM file will be sorted by chromosome number and start position in ascending order, to make file parsing in the latter steps more efficient. The SAM file headers can also be removed since they will not need to be referenced for deduping.

#### > Deduper Algorithm:

#### #function definitions 
```
store_umi_list(umi_file_path): 
  opens the umi filepath and returns a list of all the umis present in the file.
  
  return umi_list
  
  ex.) input: file contains AAA TTT CCC GGG seperated by newline charachters.
       output: return statement returns [AAA, TTT, CCC, GGG]. 
       
get_umi(sam_entry):
  given an input SAM entry, the function will parse the entry for the umi header and return it. 
  
  return umi
  
  ex.) input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC 0 2 76814284 36 71M {...} XO:Z:UU
       output: CTGTTCAC 

get_position(sam_entry):
  given a an input SAM entry, the function will parse the entry for the read position and return it.
  
  return position
  
  ex.) input: input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC 0 2 76814284 36 71M {...} XO:Z:UU
       output: 76814284 
  
get_seq(sam_entry):
   given an input SAM entry, the function will parse the entry for the sequence and return it. 
   
   return rname
   
   ex.) input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC 2 2 76814284 36 71M AATCT {...} XO:Z:UU
        output: 2AATCT

get_rname(sam_entry):
   given an input SAM entry, the function will parse the entry for the chromosome number and return it. 
   
   return rname
   
   ex.) input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC 2 2 76814284 36 71M {...} XO:Z:UU
        output: 2
 
get_strand(sam_entry):
   given an input SAM entry, the function will parse the entry for the bitwiseflag and return ('+')  
   if the strand is the reference strand and ('-') if the strand is the revcomp strand.
   
   return strand_designation
   
   ex.) input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC 0 2 76814284 16 71M {...} XO:Z:UU
        output: '-'
        
is_mapped(sam_entry):
    given an input SAM entry, the function will parse the entry for the bitwiseflag and return 'TRUE'  
    if the strand is mapped and 'FALSE' if the strand is unmapped.
    
    return boolean
    
    ex.) input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC 0 2 76814284 4 71M {...} XO:Z:UU
         output: FALSE
 
 softclip_adjust(sam_entry):
    given an input SAM entry, the function will parse the entry for the CIGAR string and determine whether   
    the alignment was soft-clipped. If it was, it will return the number of nucleotides soft-clipped. If it 
    was not, it will return 0. 
    
    return num
 
     ex.) input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC 0 2	76814284 36 3S71M {...} XO:Z:UU
          output: 3
  
 check_duplicate(sam_entry, dictionary):
    given an input SAM entry, the function will use the other checker functions to determine whether the  
    unique umi/chrom/pos/strand combination already exists in the dictionary. if it does, the function will  
    return 'TRUE' to designate this alignment as a duplicate. If it does not, the function will return 'FALSE'
    
    return boolean
   
    ex.) input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC 0 2 76814284	36 71M {...} XO:Z:UU
          output: 'TRUE' 
 
```

#### #main algorithm
```

read in SAM file from arguments
read in UMI file from arguments
open output file to write duplicates to
open output file to write bad alignments to 

umi_list = store_umi_list(umi file)

initialize dictionary to store information as key:[UMI, CHROM, POS, STRAND] and value: 'alignment sequence'

Read the current line and store in 'line' variable
write get_umi(line), get_rname(line), get_position(line), get_position(line) and get_seq(line) variables to dictionary 

While it is not the end of the file:
      
      if get_umi(line) in umi_list:
            if is_mapped(line) == 'TRUE':
                if check_duplicate(line, dict) == 'TRUE':
                      write current line to output file
                elif check_duplicate(line, dict) == 'FALSE':
                    write get_umi(line), get_rname(line), get_position(line), get_position(line) and get_seq(line) variables to dictionary 
            else:
                 write current line to bad alignments file
      else:
          write current line to bad alignments file 
                 
```
