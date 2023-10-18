#!/bin/tcsh

# python /home/user/orit/workspace/ElementCons/homer_rename_bamF.py //data//old_data//orit//ElementCons

# using RAMPAGE aligned BAM files perform peak calling with HOMER
python /home/user/orit/workspace/ElementCons/homer_findPeaksRAMPAGEbam.py /data/old_data/orit/ElementCons //private//Databases//Drosophila//Drosophila_melanogaster//UCSC//dm6

#create bed files based on the output of homer findcsRNATSS peak calling, download seq +-100 bp around TSS and run ElemeNT on it
python /home/user/orit/workspace/ElementCons/createBedFiles4ncRNAdowloadSeq.py /data/old_data/orit/ElementCons //private//Databases//Drosophila RAMPAGE

# normalize the ElemeNT output so that every motif appears in a separate line
python /home/user/orit/workspace/ElementCons/normalizeElemeNToutputGen.py /data/old_data/orit/ElementCons RAMPAGE

# prepare the input files for the R script that creates the graphs of elements position distribution and score
python /home/user/orit/workspace/ElementCons/prepMeanNpercent4graph-NascentAll.py /data/old_data/orit/ElementCons RAMPAGE
