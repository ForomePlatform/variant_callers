Main problem for the package: provide 2nd stage (high precision) for
DeNovo variant detection.

Main call:

> cd <project directory>/denovo2
> python3 -m detect.detect2 
    -I <result of 1st stage>
    -U <path to directory with ad-libraries>
    -T <path to file with locations of trio BAM-files, proband should be last one
    -O <output file>

Special option -D <file> might be used to test original previous version of algorithm

======================================
Tools to build library 
> cd <project directory>/denovo2

1. Create ad-file <result.idx> for single person

> python3 -m adlib.collect_ad /data/hg19.fasta <input-sample.bam> <result.idx>

To repair peviously created incomplete <input.idx> add it as last (4th) parameter

2. Build library from personal ad-files:

> python3 -m adlib.collect_lib make [<input.ldx>] <list of idx files>

(Use parameter <input.ldx> in case of library extension by new idx files)

3. Information for ad-files in library:

> python3 -m adlib.collect_lib.py info <result.ldx>


======================================
Test/debug stuff

* Check if hg_19 works properly
> python3 -m adlib.hg19_rd /data/hg19.fasta

* Create partial hg19 blocks
> python3 -m adlib.hg19_portion /data/hg19.fasta 1/248605354 1/146055394

* Create partial ad-file
> python3 -m adlib.ad_miner create /wgs/data/old_cases/bgm0060_wgs/results/bgm0060a2.hg19.bam bgm0062a.short.idx *.part

* Read it
> python3 -m adlib.ad_miner read bgm0062a.short.idx

* Check it
> python3 -m adlib.ad_person bgm0062a.short.idx

* Create whole ad-file
> nohup python3 -m adlib.collect_ad /data/hg19.fasta /wgs/data/old_cases/bgm0060_wgs/results/bgm0060a2.hg19.bam bgm0062a.idx > res &

* Create library
> nohup python3 -m adlib.collect_lib make ll.ldx /data/model_output/denovo/data/bgm0004a1.idx  /data/model_output/denovo/data/bgm0005a1.idx > res &

* Append library
> nohup python3 -m adlib.collect_lib make ll1.ldx ll.ldx /data/model_output/denovo/data/bgm0004a1.idx  /data/model_output/denovo/data/bgm0005u1.idx > res & 

* Get info from library
> python3 -m adlib.collect_lib info /data/model_output/denovo/idxlib/l1.ldx

* Run detect2
> python3 -m detect.detect2 -I /wgs/data/old_cases/bgm0363_wgs/de_novo_stage1_output.txt -U /data/model_output/denovo/idxlib -T /wgs/data/old_cases/bgm0363_wgs/case_bams_.txt -O st2_0363.out -D lib_dump_0363.js

* Run prev/pysam_prac_.py (to compare with original previous algorithm):
> python2 -W ignore prev/pysam_prac_.py -I /wgs/data/old_cases/bgm0363_wgs/de_novo_stage1_output.txt -U ./lib_dump_0363.js -T /wgs/data/old_cases/bgm0363_wgs/case_bams_.txt -O st2_o_0363.out

* Build tiny test mdl-ad file
python3 -u -m detect.collect_mdl --test /data/model_output/denovo/idxlib qq.admd

* Build whole test mdl-ad file (too long with whole lib data) 
python3 -u -m detect.collect_mdl /data/model_output/denovo/idxlib result.admd
======================================
Unclean files/side stories

(A) Original code of previous variant of algorithm

prev/pysam_prac_.py & prev/pysam_prac.py

The original code of previous algorithm is in pysam_prac.py. Light modification of this code pysam_prac_.py represents programm for work in our environment for purposes to compare our version of project and previous one.
The programm uses unrelated AD-data from some special intermediate file, so to run pysam_prac_.py one needs to run detect2 with special option -D before. 
Modifications of original file care of the following:
 - special retrieval of AD-data from intermediate file
 - proper sort of results
 - got rid of uninformative prints
 - proper filtering of variants
 
(B) Intermediate version of algorithm simplification
detect/alg0.py


The previous algorithm was transformed into current form by (long) serie of simplification steps. The file alg0.py represents the latest stage of simplification equivalent to the original algorithm. The next after alg0.py key observation was as follows: return value of the method T_term_calc_for_rho(_, AD) is __almost__ the same as input value AD. This fact starts the long cascade of essential simplifications. Up to now, difference in results of algorithms looks unessential - but our version is essentially more understandable. However we prefer to keep intermediate code for possible future investigations. (Before run the code mights require adaptation to the current project environment)

(C) Bug in pysam module
adlib/check_lib.py

This code was used for deep check of ad-data: we care if the data in ad-files/libraries is the same to the variant of direct use of pysam module. Difference is: for forming ad-data in files we use wide diapason of positions in pysam call, in direct use diapason contains only one position. Results for a given position might differ, and that is the bug in pysam module. Using check_lib.py we made investigation on this problem and estimate the effect of this bug. It appeared that it happened with probalility near 1/1000 and usually it happened in dirty regions (full of repetitions) near zone of N-letters. (Before run the code mights require adaptation to the current project environment)


