******************************************************************
* Author: Lucie Hutchins, Scientific Software Engineer
* Institution : The Jackson Lab (Graber Lab)
*
* Date Sept 2009
*
* Note: This directory contains all the text files generated from
*    musdiv SNP_Annotation probes alignment analysis
*******************************************************************

We have 623,124 SNPs sampled into four pairs of 25 basepair-long probes each(total 8 probes for each SNP):
Probes with probe_id ending in:
a. ) *_1  and  *_0  are the same
b. ) *_2  and  *_3  are the same
c. ) *_4  and  *_5  are the same (antisense)
d. ) *_6  and  *_7  are the same (antisense)

We aligned the SNP annotation probes against C56BL/6J build 37 using 'PASS' alignment tool (setting : -fid 90 -t 2 -gff -info_gff -no_qual )
 and generated two alignment files. One containing probes that aligned uniquely on the genome (single hit probes -> 605,532) and another for probes that have
multiple alignments either on the same chromosome or on different chromosomes.

PASS Alignment:
1. Probes with genome alignment:
   a. *_1 (*_0) : 623,079 (99.99% of total)
   b. *_2 (*_3) : 622,556  (99.9% of total)
   c. *_4 (*_5) : 623,079  (99.99% of total)
   d. *_6 (*_7) : 622,422  (99.88% of total)
   
2. Probes that aligned to the genome only once (single hit)
   Probe count:
  a. )  *_1 (*_0) : 574,184 (92.2% of aligned probes) 
  b. )  *_2 (*_3) : 575,066(92.3% of aligned probes)     
  c. )  *_4 (*_5) : 573,766 (92.1% of aligned probes)
  d. )  *_6 (*_7) : 574,655 (92.3% of aligned probes)


We then used the master single alignment file to compute the mismatches files using probes *_0,*_2,*_4,*_6 of each SNP.
with the SNP expexted pattern for these probes: 100 010 100 010
A total of 4 files generated:
1. snp_annotation_singlehit_0-2-4-6_mismatch.txt
2. snp_annotation_singlehit_0-2-4-6_alignment.txt
3. snp_annotation_singlehit_0-2-4-6_100010100010_patern.txt
4. snp_annotation_singlehit_0-2-4-6_100010100010_patern_mis.txt  


---------------------------------------------------------------------------------------------
File : snp_annotation_singlehit_0-2-4-6_mismatch.txt (605,532 SNPs)
Description :
     This file (commas separated file) contains the mismatches info of all
       the probes for each SNP
     The file was generated running the program generate_table.pl using the main
     single hit file of the snpAnnotation probes (snp_annotation_singlehit.psl)

Note: the _0-2-4-6_ in the file name specify the probes# (probes *_0,*_2,*_4, and *_6 were used)

Fields:
   1. SNPID,
   2. P_0-0mm,
   3. P_0-1mm,
   4. P_0-2mm,
   5. P_2-0mm,
   6. P_2-1mm,
   7. P_2-2mm,
   8. P_4-0mm,
   9. P_4-1mm,
   10.P_4-2mm,
   11.P_6-0mm,
   12.P_6-1mm,
   13.P_6-2mm
   14.has_expected_pattern [0/1]
---------------------------------------------------------------------------------------
File : snp_annotation_singlehit_0-2-4-6_alignment.txt (605,532 SNPs)
Description : 
     This file (commas separated file) contains the alignment info of all
       the probes for each SNP
     The file was generated running the program generate_table.pl using the main
     single hit file of the snpAnnotation probes (snp_annotation_singlehit.psl)

Note: the _0-2-4-6_ in the file name specify the probes# (probes *_0,*_2,*_4, and *_6 were used)

Fields:
   1. SNPID,
   2. snpChr,
   3. snpPos,
   4. P_0-chr,
   5. P_0-tSart,
   6. P_0-tEnd,
   7. P_0_alignR,
   8. P_2-chr,
   9. P_2-tSart,
   10.P_2-tEnd,
   11.P_2_alignR,
   12.P_4-chr,
   13.P_4-tSart,
   14.P_4-tEnd,
   15.P_4_alignR,
   16.P_6-chr,
   17.P_6-tSart,
   18.P_6-tEnd,
   19.P_6_alignR 

-----------------------------------------------------------------------------
File : snp_annotation_singlehit_0-2-4-6_100010100010_patern.txt (521,928 SNPs)
Description :
     This file (commas separated file) contains the mismatches info of all
     the probes for each SNP and includes cases where all the probes of a given
     SNP have the expexted pattern : 100 010 100 010
       
     The file was generated running the program generate_table.pl using the main
     single hit file of the snpAnnotation probes (snp_annotation_singlehit.psl)

Note: the _0-2-4-6_ in the file name specify the probes# (probes *_0,*_2,*_4, and *_6 were used)

Fields:
   1. SNPID,
   2. P_0-0mm,
   3. P_0-1mm,
   4. P_0-2mm,
   5. P_2-0mm,
   6. P_2-1mm,
   7. P_2-2mm,
   8. P_4-0mm,
   9. P_4-1mm,
   10.P_4-2mm,
   11.P_6-0mm,
   12.P_6-1mm,
   13.P_6-2mm
-----------------------------------------------------------------------------
File : snp_annotation_singlehit_0-2-4-6_100010100010_pattern_mis.txt (83,603 SNPs)
Description :
     This file (commas separated file) contains the mismatches info of all
     the probes for each SNP and includes cases one or more probes of a given
     SNP do not have the expexted pattern : 100 010 100 010
       
     The file was generated running the program generate_table.pl using the main
     single hit file of the snpAnnotation probes (snp_annotation_singlehit.psl)

Note: the _0-2-4-6_ in the file name specify the probes# (probes *_0,*_2,*_4, and *_6 were used)

Fields:
   1. SNPID,
   2. P_0-0mm,
   3. P_0-1mm,
   4. P_0-2mm,
   5. P_2-0mm,
   6. P_2-1mm,
   7. P_2-2mm,
   8. P_4-0mm,
   9. P_4-1mm,
   10.P_4-2mm,
   11.P_6-0mm,
   12.P_6-1mm,
   13.P_6-2mm
