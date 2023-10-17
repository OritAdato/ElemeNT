
def get_specie_name(basefn):
    basefn_l = basefn.split('.')
    specie = basefn_l[1]
    return (specie)

def create_chippeak_dir(specie_sga_dir,dir2create):
    import os
    if os.path.isdir(os.path.join(specie_sga_dir, dir2create)):
        pass
   #     print('specie dir: ', specie_sga_dir,dir2create, 'Exists')
    else:
        os.makedirs(os.path.join(specie_sga_dir, dir2create),exist_ok=True)
    chippeak_dir=os.path.join(specie_sga_dir, dir2create)
    return (chippeak_dir)

def get_specie_fasta_dir(specie_name):
    if specie_name == 'dp3':
        specie_fasta_dir = "//private//Databases//Drosophila//Drosophila_pseudoobscura//UCSC//dp3//WholeGenome//dp3.fa"
    else:
        if specie_name == 'dm6':
            specie_fasta_dir = "//private//Databases//Drosophila//Drosophila_melanogaster//UCSC//dm6//Sequence//WholeGenomeFasta//genome.fa"
        else:
            if specie_name == 'droAna3':
                specie_fasta_dir = "//private//Databases//Drosophila//Drosophila_ananassae//UCSC//droAna3//WholeGenomeFasta//droAna3.fa"
            else:
                if specie_name == 'droSim1':
                    specie_fasta_dir = "//private//Databases//Drosophila//Drosophila_simulans//UCSC//droSim1//WholeGenomeFasta//droSim1.fa"
                else:
                    if specie_name == 'droEre2':
                        specie_fasta_dir = "//private//Databases//Drosophila//Drosophila_erecta//UCSC//droEre2//WholeGenomeFasta//droEre2.fa"
                    else:
                        specie_fasta_dir = ' '


    return specie_fasta_dir

#The following proc (split_longName, ) were copied from prepCommonGeneElemetsTable.py
def split_longName(elemeNTout_devStage_df):
    elemeNTout_devStage_df.rename(columns={'name':'longName'}, inplace=True)
    elemeNTout_devStage_df['splitPos'] = elemeNTout_devStage_df['longName'].str.find('::') + 2
    elemeNTout_devStage_df['name'] =  elemeNTout_devStage_df.apply(lambda x: x['longName'][x['splitPos']:],axis=1)
    return elemeNTout_devStage_df

def prepare_ElemeNT_config_final(fasta_name,output_name,gc,ElemeNT_dir,startConstant,smooth,order):
    # Prepare the Config.txt file that will be use in the ElemeNT run. The file is saved under ElemeNT_dir
    import os
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
                                if (line[:5] =='Order' and order !=''):
                                    output_line = 'Order: ' + order + '\n'
                                    config_f.write(output_line)
                                else:
                                    config_f.write(line)
    config_f.close()
    return

def handel_empty_pos(motifName):
    if (motifName == 'dInr' or motifName == 'DPE' or motifName == 'DPEhInr' or motifName == 'DPEdInr' or
            motifName == 'DPEBBCABW'):
        if motifName == 'dInr':
            pos_l = [-12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        else:
            pos_l = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39]
    else:
        if motifName =='TATA':
            pos_l = [-40, -39, -38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20]
        else:
            if motifName == 'Motif1':
                pos_l = [-17, -16,-15, -14, -13,-12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4]
            else:
                if motifName == 'dTCT':
                    pos_l = [-12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9]
                else:
                    print('4 motif:', motifName, 'no pos Exiting prog!!!!!!!!!!!!!!!')
                    exit()
    return pos_l

def unique(list1):
    import numpy as np
    x = np.array(list1)
    print(np.unique(x))
    ul = np.unique(x).tolist()
    return ul

def update_DPErealPos(elemeNTout_4motif_df):
    elemeNTout_4motif_rows = elemeNTout_4motif_df.index.tolist()
    for row_i in elemeNTout_4motif_rows:
        pos = elemeNTout_4motif_df.loc[row_i,'pos']
        if pos < 0:  # the difference is since there is no 0
            elemeNTout_4motif_df.at[row_i,'pos'] = pos + 30
        else:
            elemeNTout_4motif_df.at[row_i, 'pos'] = pos + 29
    return elemeNTout_4motif_df