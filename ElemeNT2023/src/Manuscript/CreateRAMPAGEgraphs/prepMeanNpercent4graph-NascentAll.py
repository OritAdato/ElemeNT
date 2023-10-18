import glob
import os
import pandas as pd
from genProc import *
import sys
# This program prepares the input file for creation of Graphs that show for every motif for all nascent RNA dev windows
# the summary of % of detected motifs in every position and their mean score. The files created by this program are
# used as input to an R script that creates the graph.
#
def readAllTimeDev(motifElemeNT_out_df, specie_ElemeNTannotNorm_dir, elemeNToutPref, motifName):
    for ElemenNT_outFile in glob.glob(os.path.join(specie_ElemeNTannotNorm_dir, elemeNToutPref +'*.csv')):
        ElemenNT_out_fn = os.path.basename(ElemenNT_outFile)
        dev_stage = ElemenNT_out_fn[ElemenNT_out_fn.find('-')-1:ElemenNT_out_fn.find('-')+2]
        # print(ElemenNT_out_fn,dev_stage)
        ElemeNT_out_df = pd.read_csv(ElemenNT_outFile,sep=',',header=0)
        if 'DPE' in motifName:
            if len(motifName) > 3:
                dependent_inrName = motifName[3:]
                elemeNTout_4motif_df = ElemeNT_out_df.loc[((ElemeNT_out_df['motifName'] == 'DPE') &
                                                           (ElemeNT_out_df['InrName'] == dependent_inrName)), :]
                # print('DPE', dependent_inrName)
            else:
                elemeNTout_4motif_df = ElemeNT_out_df.loc[ElemeNT_out_df['motifName'] == motifName, :]
            elemeNTout_4motif_df = update_DPErealPos(elemeNTout_4motif_df)
        else:
            elemeNTout_4motif_df = ElemeNT_out_df.loc[ElemeNT_out_df['motifName'] == motifName, :]
        motifElemeNT_out_df = pd.concat([motifElemeNT_out_df,elemeNTout_4motif_df])

    return  motifElemeNT_out_df


def createMeanNpercentNascent_df(motifName):
    #colnames = ['0-2','2-4','4-6','6-8']
    colnames = ["DevStages_0_to_8"]
    meanScore_df = pd.DataFrame(columns=colnames)
    precent_df = pd.DataFrame(columns=colnames)
    return meanScore_df, precent_df

def prepFiles4graph(elemeNTout_4motif_df,motifName,devStage,motifMean_df, motifPercent_df,specie_name):
    import  numpy as np
    pos_l = elemeNTout_4motif_df['pos'].values.tolist()
    unique_val = unique(pos_l)
    if len(unique_val) != 21:
        print('4 motif:', motifName, 'len of pos list is: ', str(len(unique_val)))
        unique_val = handel_empty_pos(motifName)
    else:
        print('motifName:', motifName, unique_val)
    if len(motifMean_df.index.tolist()) == 0:
        indTitle = specie_name + '-' + motifName + '_meanScore'
        motifMean_df =add_emptyPosIndexes (unique_val,motifMean_df,indTitle)
        indTitle = specie_name + '-' + motifName + '_percent'
        motifPercent_df = add_emptyPosIndexes (unique_val,motifPercent_df,indTitle)

    for pos in motifMean_df.index:
        scoreInpos_l = elemeNTout_4motif_df.loc[elemeNTout_4motif_df['pos']==pos,'score'].tolist()
        if len(scoreInpos_l) == 0:
            motifMean_df.at[pos, devStage] = 0
        else:
            meanScoreInpos = np.mean(scoreInpos_l)
            motifMean_df.at[pos,devStage] = round(meanScoreInpos,1)

        totNumOfElement = len(elemeNTout_4motif_df.index.tolist())
        if totNumOfElement == 0:
            motifPercent_df.at[pos, devStage] = 0
        else:
            numOfElementInpos = len(scoreInpos_l)
            percentInpos = round((numOfElementInpos/totNumOfElement)*100,1)
            motifPercent_df.at[pos,devStage] = percentInpos
    return motifMean_df, motifPercent_df


def add_emptyPosIndexes (pos_l,motifMean_df,indTitle):

    colnames = motifMean_df.columns.tolist()
    motifMeanWindexes_df = pd.DataFrame(columns=colnames, index=pos_l)
    motifMeanWindexes_df.index.name = indTitle

    return motifMeanWindexes_df

print(sys.argv[0])
for i in range(1, len(sys.argv)):
    print('argument:', i, 'value:', sys.argv[i])
    if i== 1:
        dir_mainPath = sys.argv[i]
    else:
        if i==2: # i=2
            inputFiles_type = sys.argv[i]

width=100

specie_name = 'dm6'

if inputFiles_type == 'RAMPAGE':
    dir_mainPath = os.path.join(dir_mainPath,"scripts/pipelineOut/GSE89299.dm6/HOMER")
    specie_dir = os.path.join(dir_mainPath, "findcsRNATSSoutput")
    middle_fn = "_HomerRampageAllDS_"
    elemeNToutPref = "tssAnalysis_GSM"
else:
    specie_dir = os.path.join(dir_mainPath, "csRNAseq")
    middle_fn = "_nascentRNAallDS_"
    elemeNToutPref = "cas9"

motifName_l = ['dInr','DPE', 'TATA']
motifName_l = ['dInr','DPE','TATA','PB','Motif1','dTCT','GAGA']
specie_name = 'dm6'

specie_ElemeNTannotNorm_dir = os.path.join(specie_dir,"ElemeNTout//normlized")

print ('specie dir: ', specie_name )

normalizedMotifs_dir = os.path.join(specie_dir,'ElemeNTout')
targetGraphInput_dir = create_chippeak_dir(normalizedMotifs_dir, 'graphInputFiles')
normalizedElemeNTout_colnames_l = ['name','motifName','pos','score','InrScore','InrName','motif_seq']

for motifName in motifName_l:
    print ('motifName:',motifName)
    motifMean_df, motifPercent_df = createMeanNpercentNascent_df(motifName)
    # For every specie and every dev stage reading the output of ElemeNT:
    motifElemeNT_out_df = pd.DataFrame(columns=normalizedElemeNTout_colnames_l)
    motifElemeNT_out_df = readAllTimeDev(motifElemeNT_out_df, specie_ElemeNTannotNorm_dir, elemeNToutPref, motifName)
    dev_stage = "DevStages_0_to_8"
    motifMean_df, motifPercent_df = prepFiles4graph(motifElemeNT_out_df,motifName,dev_stage,motifMean_df, motifPercent_df,specie_name)
    meanScore_out_fn = os.path.join(targetGraphInput_dir,specie_name + middle_fn + motifName + '_meanScore.txt')
    percent_out_fn = os.path.join(targetGraphInput_dir, specie_name + middle_fn + motifName + '_percent.txt')
    motifMean_df.to_csv(meanScore_out_fn,index=True,header=True,sep='\t')
    motifPercent_df.to_csv(percent_out_fn,index=True,header=True,sep='\t')








