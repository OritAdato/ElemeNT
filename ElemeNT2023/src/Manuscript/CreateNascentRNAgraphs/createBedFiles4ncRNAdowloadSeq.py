import os
import glob
from genProc import *
import pandas as pd
import pybedtools

# This program uses as input the TSS peak calling output files created by the analysis of nascent RNA data using
# using HOMER csRNA-seq Analysis utility. Based on the TSS peaks it created bed files for downloading +- bp around
# TSS. Using pybedtools it downloads the genomic sequences, adds them the gene annotations (annotated by HOMER)
# and send them as input to ElemeNT utility for searching corepromoter elements
#
def isnan(value):
    try:
        import math
        return math.isnan(float(value))
    except:
        return False
def get_specie_fasta_dir(specie_name):
    if specie_name == 'dp3':
        specie_fasta_dir = "//SpecieGenomePath//Databases//Drosophila//Drosophila_pseudoobscura//UCSC//dp3//WholeGenome//dp3.fa"
    else:
        if specie_name == 'dm6':
            specie_fasta_dir = "//SpecieGenomePath//Databases//Drosophila//Drosophila_melanogaster//UCSC//dm6//Sequence//WholeGenomeFasta//genome.fa"
        else:
            if specie_name == 'droAna3':
                specie_fasta_dir = "//SpecieGenomePath//Databases//Drosophila//Drosophila_ananassae//UCSC//droAna3//WholeGenomeFasta//droAna3.fa"
            else:
                if specie_name == 'droSim1':
                    specie_fasta_dir = "//SpecieGenomePath//Databases//Drosophila//Drosophila_simulans//UCSC//droSim1//WholeGenomeFasta//droSim1.fa"
                else:
                    if specie_name == 'droEre2':
                        specie_fasta_dir = "//SpecieGenomePath//Databases//Drosophila//Drosophila_erecta//UCSC//droEre2//WholeGenomeFasta//droEre2.fa"
                    else:
                        specie_fasta_dir = ' '


    return specie_fasta_dir

def create_CoordinatesBed(ncRNA_tssOnly_df,ncRNAinput_fn,targetBed_dir):

    bedfn = ncRNAinput_fn[:-4] + 'Only.bed'
    target_fn = os.path.join(targetBed_dir,bedfn)
    target_f = open(target_fn,'w')
    row_l = ncRNA_tssOnly_df.index.tolist()
    for rowInd in row_l:
        chrmsm_name = ncRNA_tssOnly_df.loc[rowInd,'chr']
        from_pos = ncRNA_tssOnly_df.loc[rowInd,'start'] -1 -25
        to_pos = ncRNA_tssOnly_df.loc[rowInd,'end'] +25
        entry_name = ncRNA_tssOnly_df.loc[rowInd,'#tssClusterID']
        score = ncRNA_tssOnly_df.loc[rowInd,'score']
        strand = ncRNA_tssOnly_df.loc[rowInd,'strand']
        line4print = (chrmsm_name + '\t' + str(from_pos) + '\t' + str(to_pos) + '\t' + entry_name + '\t' + str(score) +
                              '\t' + strand + '\n' )
        target_f.write(line4print)
    target_f.close()
    print('created',target_fn)
    return target_fn


def extract_seq4bed(target_fn,fasta_fn,targetFasta_dir):
    fasta_out_f = os.path.basename(target_fn)
    fasta_out_f = fasta_out_f[:-3] + 'fasta'
    fasta_out_fn = os.path.join(targetFasta_dir,fasta_out_f)
    bedfile_fn = pybedtools.BedTool(target_fn)       # wrapper for bedtools getfasta
    bedfile_fn = bedfile_fn.sequence(fi=fasta_fn, fo=fasta_out_fn, s=True)   # s=forces trandness removed name=True,
    return fasta_out_fn

def addGeneName2fasta(fasta_WithN_fn,ncRNA_df):
    fasta_WithGN_fn = fasta_WithN_fn[:-6] + '_wGN.fasta'
    fasta_WithGN_f = open(fasta_WithGN_fn,'w')
    with open (fasta_WithN_fn,'r') as fasta_WithN_f:
        for line in fasta_WithN_f:
            if line[0] == '>':
                chr, fromPos, toPos, strand = getGeneSearchFields(line)
                gene_name_l = ncRNA_df.loc[((ncRNA_df['chr']==chr) & (ncRNA_df['start']==int(fromPos) + 1 + 25) &
                              (ncRNA_df['end']==int(toPos) -25 ) & (ncRNA_df['strand']==strand)),'Closest TSS transcript'].tolist()
                if len(gene_name_l) !=1:
                    print ('len is != 1', chr,fromPos,toPos,strand,gene_name_l)
                    gene_name =gene_name_l[0]

                else:
                    gene_name = gene_name_l[0]
                    gn_nan = isnan(gene_name)
                    if gn_nan:
                        gene_name = ''
                line4print = line.rstrip() + gene_name + '\n'
                fasta_WithGN_f.write(line4print)
            else:
                fasta_WithGN_f.write(line)
    fasta_WithGN_f.close()

    return fasta_WithGN_fn

def getGeneSearchFields(line):
    # >chr3R:23026638-23026788(-)
    chr = line[1:line.find(':')]
    fromPos = line[line.find(':')+1:line.find('-')]
    toPos = line[line.find('-')+1:line.find('(')]
    strand = line[line.find('(')+1:line.find(')')]
    return chr, fromPos, toPos, strand

def prepare_ElemeNT_config_final(fasta_name,output_name,gc,ElemeNT_dir,startConstant,smooth):
    # Prepare the Config.txt file that will be use in the ElemeNT run. The file is saved under ElemeNT_dir
    os.chdir(ElemeNT_dir) # change current directory to ElemeNT binary dir

    dm_cutoff_fn = os.path.join(ElemeNT_dir, 'ConfigOrig.txt')
    config_fn = os.path.join(ElemeNT_dir,'Config.txt')
    gc = round(gc, 2)
    gc_content = 'GC: ' + str(gc) + ';\n'
    config_f = open(config_fn,'w')

    with open (dm_cutoff_fn,'r') as dm_cutoff_f:
        for line in dm_cutoff_f:
            if line[:13] == 'InputFileName':
                input_line = 'InputFileName: ' + fasta_name + ';\n'
                config_f.write(input_line)
            else:
                if line[:14] == 'OutputFileName':
                    output_line = 'OutputFileName: ' + output_name + ';\n'
                    config_f.write(output_line)
                else:
                    if line[:2]=='GC':
                        config_f.write(gc_content)
                    else:
                        if line[:13] == 'StartConstant':
                            output_line = 'StartConstant: ' + str(startConstant) + ';\n'
                            config_f.write(output_line)
                        else:
                            if line[:14] == 'SmoothConstant':
                                output_line = 'SmoothConstant: ' + str(smooth) + ';\n'
                                config_f.write(output_line)
                            else:
                                # Order: hInr, dInr, GAGA, BREu, TATA, BREd, XCPE1, XCPE2, Motif1, hTCT, dTCT, PB, BBCABW;
                                if line[:5] =='Order':
                                    output_line = 'Order: dInr, GAGA, BREu, TATA, BREd, XCPE1, XCPE2, Motif1, dTCT, PB' + ';\n'
                                    config_f.write(output_line)
                                else:
                                    config_f.write(line)
    config_f.close()
    return

dir_mainPath = "//Main//Path"
input_dir = os.path.join(dir_mainPath,"csRNAseq")
width = 150

specie_name = 'dm6'
gc = 42.0283    # this is the gc content of dm6

startConstant = 100
smooth = 10

targetFasta_dir = create_chippeak_dir(input_dir, 'ncRNAfasta_w' + str(width))
ElemeNT_out_dir = create_chippeak_dir(input_dir, 'ElemeNTout')
for ncRNAfn in glob.glob(os.path.join(input_dir, 'cas9*.tss.txt')):
    basefn = os.path.basename(ncRNAfn)
    targetBed_dir = create_chippeak_dir(input_dir, 'bed4fasta')
    ncRNA_df = pd.read_csv(ncRNAfn,header=0, sep='\t')
    ncRNA_tssOnly_df = ncRNA_df.loc[ncRNA_df['annotation']=='tss',:]
    fasta_fn = get_specie_fasta_dir(specie_name)
    target_fn = create_CoordinatesBed(ncRNA_tssOnly_df,basefn,targetBed_dir)
    #print (target_fn, '\n fasta:', fasta_fn)
    fasta_WithN_fn = extract_seq4bed(target_fn,fasta_fn,targetFasta_dir)
    fasta_WgeneName_fn = addGeneName2fasta(fasta_WithN_fn,ncRNA_df)

    basefn = os.path.basename(fasta_WithN_fn)

    output_fn = basefn[:-6] + '_start' + str(startConstant) + 'smooth' + str(smooth) + '.txt'
    output_name = os.path.join(ElemeNT_out_dir,output_fn)
    ElemeNT_dir = os.path.join(input_dir,'ElemeNT_V22')
    prepare_ElemeNT_config_final(fasta_WgeneName_fn, output_name, gc, ElemeNT_dir,startConstant,smooth)
    #print('current dir is: ', os.getcwd())
    os.system('ElemeNT_V2023_binary')
