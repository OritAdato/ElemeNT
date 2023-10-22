#!/bin/tcsh

# python homer_rename_bamF.py /ReplaceRootPath


#create bed files based on the output of homer findcsRNATSS peak calling, download seq +-100 bp around TSS and run ElemeNT on it
python createBedFiles4ncRNAdowloadSeq.py /ReplaceRootPath //genomes//Databases//Drosophila nascent

# normalize the ElemeNT output so that every motif appears in a separate line
python normalizeElemeNToutputGen.py /ReplaceRootPath nascent

# prepare the input files for the R script that creates the graphs of elements position distribution and score
python prepMeanNpercent4graph-NascentAll.py /ReplaceRootPath nascent
