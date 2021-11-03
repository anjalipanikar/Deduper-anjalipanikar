# Deduplexing Answers

### Anjali Panikar 


### ***Pre-Processing*** 

Reads were sorted and headers were translocated from the SAM file to the output file. Samtools sort was used to achieve this. Commands were executed within the main python script (panikar_deduper.py) using the os.system() function.

```
$ samtools view -S -b <sam file> > <bam file>  
$ samtools sort <bam file> -o <sorted bam file>  
$ samtools view -h <sorted bam file> > <sorted sam file>   
$ "grep '^[[:space:]]*@' sorted_sam > final_sam  
$ "grep -v '^[[:space:]]*@' sorted_sam > sorted_no_headers  
```  
### ***Main Algorithm Outputs***
  
4 files are outputted by the main algorithm:     
(1) final SAM file with deduplexed reads and all headers (.deduplexed.sam)  
(2) file with duplicate reads (.duplicates.sam)  
(3) file with unmapped reads (.unmapped.sam)  
(4) file with bad UMIs that don't match the UMI file (.bad_umis.sam)  

Counts from outputs:    

Unique reads: 13,719,048   
command: ```cat C1_SE_uniqAlign.deduplexed.sam | grep -v "^@" | wc -l```    

Header lines: 65  
command: ```cat C1_SE_uniqAlign_deduped.sam | grep -E "^@" | wc -l```    

Duplicate reads: 4,467,362  
command: ```cat C1_SE_uniqAlign.duplicates.sam | wc -l```    

Wrong UMIs: 0
command: ```cat C1_SE_uniqAlign.bad_umis.sam | wc -l```

Unmapped reads: 0
command: ```cat C1_SE_uniqAlign..unmapped.sam | wc -l```

Uniq reads per chromosome: 
```
Number of Reads   Chromosome
697508 1
564903 10
1220389 11
359951 12
467659 13
387239 14
437465 15
360923 16
517566 17
290506 18
571665 19
2787018 2
547615 3
589839 4
562160 5
510818 6
1113183 7
576463 8
627488 9
202002 MT
317853 X
2247 Y
3 JH584299.1
656 GL456233.2
6 GL456211.1
4 GL456221.1
1 GL456354.1
5 GL456210.1
4 GL456212.1
294 JH584304.1
2 GL456379.1
3 GL456367.1
1 GL456239.1
1 GL456383.1
5450 MU069435.1
1 GL456389.1
21 GL456370.1
1 GL456390.1
1 GL456382.1
17 GL456396.1
3 GL456368.1
3 MU069434.1
111 JH584295.1
```
command: ```cat C1_SE_uniqAlign.deduplexed.sam | grep -v "^@" | cut -f 3 | uniq -c```



