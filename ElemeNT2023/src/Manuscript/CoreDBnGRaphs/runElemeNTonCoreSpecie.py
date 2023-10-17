import os
import glob
from genProc import *
import pandas as pd

def remove_seqWithNfromfasta (targetFasta_dir,fasta_WithN_f):
    fasta_WithN_fn = os.path.join(targetFasta_dir,fasta_WithN_f)
    #fasta_woN_f = fasta_WithN_f[:-6]+'_woN.fasta'
    fasta_woN_newfn = os.path.join(targetFasta_dir,fasta_WithN_f[:-7]+'_woN.fasta')
    removed_line =0
    fasta_woN_newf = open(fasta_woN_newfn,'w')
    seq_l_l =[]
    gene_names_l = []
    with open(fasta_WithN_fn,'r') as fasta_file:
        first_seq = True
        remove_seq = False
        for line in fasta_file:
            if line [0]== '>':
                line_l = line.split(' ')
                gene_name = line_l[1]     # skip the uuid of EPDnew and leave only gene name but make sure it is unique
                if gene_name[len(gene_name)-1:] == '_':
                    gene_name = gene_name + '1'
                if gene_name in gene_names_l:            # if the geneName already exist it means that this is another
                    gene_name = add_unique2gn(gene_name) # transcript - add seq num to its name so that it will be unique
                gene_names_l.append(gene_name)
                line = '>' +  gene_name
                for ind in range(2,len(line_l)):
                    line = line + ' ' + line_l[ind]
                if first_seq:
                    first_seq = False
                    title_l = line
                else:
                    if remove_seq:
                        print ('***removed: ', title_l[:20])
                        pass
                    else:
                        fasta_woN_newf.write(title_l)
                        for ind in range(len(seq_l_l)):
                            fasta_woN_newf.write(seq_l_l[ind])
                    title_l = line
                    remove_seq = False
                    seq_l_l = []
            else:
                seq_l = line
                if seq_l[0] != '>':
                    if (('N' in seq_l) or ('n' in seq_l)):
                        removed_line +=1
                        seq_l_l = []
                        remove_seq = True
                    else:
                        seq_l_l.append(seq_l)
    if remove_seq:
        print('***removed: ', title_l[:20])
        pass
    else:
        fasta_woN_newf.write(title_l)
        for ind in range(len(seq_l_l)):
            fasta_woN_newf.write(seq_l_l[ind])
    fasta_woN_newf.close()
    # if os.path.exists(fasta_WithN_fn):
    #     os.remove(fasta_WithN_fn)
    return fasta_woN_newfn

def add_unique2gn(gene_name):
    gene_suffix = gene_name[len(gene_name)-1:]
    gene_suffix_n = int(gene_suffix) + 1
    gene_name = gene_name[:len(gene_name)-1] + str(gene_suffix_n)
    return gene_name

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
                                config_f.write(line)
    config_f.close()
    return
#
# This program runs the ElemeNT utility on fasta sequences (+-100 bp around TSS) of 10 species that were downloaded
# from EPDNEW. The name pattern of sequences file for every specie is:
# specieName_EPDnew_minPlus100_seq.txt
# the program uses the table species_gc.csv that includes the GC content of every specie
#
dir_mainPath = "//Main//Path"
input_dir = os.path.join(dir_mainPath,"CoreSeq//sequences052023")
ElemeNT_dir = os.path.join(dir_mainPath,"CoreSeq//ElemeNT_V2023")
specie_gc_fn = os.path.join(input_dir,'species_gc.csv')
width = 100
startConstant=100
smooth = 10
# specie_name = 'dm6'
#gc = 42.0283    # this is the gc content of dm6
specie_gc_df = pd.read_csv(specie_gc_fn,header=0, sep=',')

ElemeNT_out_dir = create_chippeak_dir(input_dir, 'ElemeNTout')
for specieSeqfn in glob.glob(os.path.join(input_dir, '*EPDnew_minPlus100_seq.txt')):
    basefn = os.path.basename(specieSeqfn)
    specie_name = basefn[:basefn.find('_')]
    specie_gc = specie_gc_df.loc[specie_gc_df['EPDnew version']==specie_name,'%GC'].item()
    print ('specie name: ', specie_name, ' specie gc:', str(specie_gc) )
    fasta_woN_newfn = remove_seqWithNfromfasta(input_dir, basefn)
    output_fn = basefn[:-7] + '_start' + str(startConstant) + 'smooth' + str(smooth) +'.txt'
    output_name = os.path.join(ElemeNT_out_dir,output_fn)
    prepare_ElemeNT_config_final(fasta_woN_newfn, output_name, specie_gc, ElemeNT_dir,startConstant,smooth)
    #print('current dir is: ', os.getcwd())
    os.system('ElemeNT_V2023_binary')
