import glob
import os
import pandas as pd
from genProc import *
import numpy as np

# This program prepares the input file for creation of Graphs that show for every motif for every specie the summary
# of % of detected motifs in every position and thier mean score. The files created by this program are used as
# input to an R script that creates the graph.
#
def createMeanNpercentCore_df():
    colnames = ['worm','honeybee','fly','zebrafish','chicken','dog','rat','mouse','macaque','human']

    meanScore_df = pd.DataFrame(columns=colnames)
    precent_df = pd.DataFrame(columns=colnames)
    return meanScore_df, precent_df

def prepFiles4graph(ElemeNT_out_df,motifName,devStage,motifMean_df, motifPercent_df,specie_name):
    import  numpy as np
    if 'DPE' in motifName:
        if len(motifName) > 3:
            dependent_inrName = motifName[3:]
            elemeNTout_4motif_df = ElemeNT_out_df.loc[((ElemeNT_out_df['motifName'] == 'DPE') &
                                              (ElemeNT_out_df['InrName'] == dependent_inrName)), :]
            print('DPE', dependent_inrName)
        else:
            elemeNTout_4motif_df = ElemeNT_out_df.loc[ElemeNT_out_df['motifName'] == motifName, :]

        elemeNTout_4motif_df = update_DPErealPos(elemeNTout_4motif_df)
    else:
        elemeNTout_4motif_df = ElemeNT_out_df.loc[ElemeNT_out_df['motifName']==motifName,:]
    pos_l = elemeNTout_4motif_df['pos'].values.tolist()
    unique_val = unique(pos_l)
    if len(unique_val) != 21:
        print('4 motif:', motifName, 'len of pos list is: ', str(len(unique_val)))
        unique_val = handel_empty_pos(motifName)
    else:
        print('motifName:', motifName, unique_val)
    if len(motifMean_df.index.tolist()) == 0:
        indTitle = motifName + '_meanScore'
        motifMean_df =add_emptyPosIndexes (unique_val,motifMean_df,indTitle)
        indTitle = motifName + '_percent'
        motifPercent_df = add_emptyPosIndexes (unique_val,motifPercent_df,indTitle)

    for pos in motifMean_df.index:
        scoreInpos_l = elemeNTout_4motif_df.loc[elemeNTout_4motif_df['pos']==pos,'score'].tolist()
        if len(scoreInpos_l) == 0:
            motifMean_df.at[pos, specie_name] = 0
        else:
            meanScoreInpos = np.mean(scoreInpos_l)
            motifMean_df.at[pos,specie_name] = round(meanScoreInpos,1)

        totNumOfElement = len(elemeNTout_4motif_df.index.tolist())
        if totNumOfElement == 0:
            motifPercent_df.at[pos, specie_name] = 0
        else:
            numOfElementInpos = len(scoreInpos_l)
            percentInpos = round((numOfElementInpos/totNumOfElement)*100,1)
            motifPercent_df.at[pos,specie_name] = percentInpos
    return motifMean_df, motifPercent_df

def add_emptyPosIndexes (pos_l,motifMean_df,indTitle):

    colnames = motifMean_df.columns.tolist()
    motifMeanWindexes_df = pd.DataFrame(columns=colnames, index=pos_l)
    motifMeanWindexes_df.index.name = indTitle

    return motifMeanWindexes_df


dir_mainPath = "//Main//path"
motifsInfo_dir = os.path.join(dir_mainPath,'scripts//Motifs')

width=100

specie_dict = {'ce6':'worm','amel5':'honeybee','dm6':'fly','danRer7':'zebrafish','galGal5':'chicken',
               'canFam3':'dog','rn6':'rat','mm10':'mouse','rheMac8':'macaque','hg38':'human'}

pipelineOut_dir = os.path.join(dir_mainPath,"scripts//pipelineOut")

motifName_l = ['dInr','DPE','DPEhInr','DPEdInr','DPEBBCABW', 'TATA','GAGA','PB','hInr','BBCABW','Motif1', 'dTCT']

specie_dir = os.path.join(dir_mainPath,"CoreSeq//sequences052023")
specie_ElemeNTannotNorm_dir = os.path.join(specie_dir,"ElemeNTout//normlized")


normalizedMotifs_dir = os.path.join(specie_dir,'ElemeNTout')
targetGraphInput_dir = create_chippeak_dir(normalizedMotifs_dir, 'graphInputFiles')

for motifName in motifName_l:
    print ('motifName:',motifName)
    motifMean_df, motifPercent_df = createMeanNpercentCore_df()
    # For every Motif and every specie reading the output of ElemeNT:
    for ElemenNT_outFile in glob.glob(os.path.join(specie_ElemeNTannotNorm_dir, '*_norm.csv')):
        ElemenNT_out_fn = os.path.basename(ElemenNT_outFile)
        specieShortName = ElemenNT_out_fn[:ElemenNT_out_fn.find('_')]
        specie_name = specie_dict[specieShortName]
        print(ElemenNT_out_fn,specieShortName,specie_name)
        ElemeNT_out_df = pd.read_csv(ElemenNT_outFile,sep=',',header=0)
        motifMean_df, motifPercent_df = prepFiles4graph(ElemeNT_out_df,motifName,specie_name,motifMean_df,
                                                        motifPercent_df,specie_name)
    meanScore_out_fn = os.path.join(targetGraphInput_dir,'all_Core_' + motifName + '_meanScore.txt')
    percent_out_fn = os.path.join(targetGraphInput_dir, 'all_Core_' + motifName + '_percent.txt')
    motifMean_df.to_csv(meanScore_out_fn,index=True,header=True,sep='\t')
    motifPercent_df.to_csv(percent_out_fn,index=True,header=True,sep='\t')







