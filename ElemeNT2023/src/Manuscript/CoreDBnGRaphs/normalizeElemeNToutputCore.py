import glob
import os
import pandas as pd
from genProc import *

def create_summary_df(specie_name):
    colnames =['fileN','devStage','smooth','GAGA','BREu','TATA','BREd','XCPE1','XCPE2','Motif1','dTCT','hTCT','BBCABW','hInr',
           'dInr','MTE','bridge','DPE','PB']
    sum_df = pd.DataFrame(columns=colnames)
    return sum_df

def summarize_ElemeNT_out(ElemeNT_out_df,ElemenNT_outFile,sum_df):
    basefn = os.path.basename(ElemenNT_outFile)
    sum_df.loc[len(sum_df.index), 'fileN'] = basefn
    dev_stage = basefn[5:9]
    sum_df.loc[sum_df['fileN'] == basefn, 'devStage'] = dev_stage    # nc dev_Stage
    sum_df.loc[sum_df['fileN'] == basefn, 'smooth'] = 0 # nc smooth
    #print("time is: ", str(timeOfDev), "dev: ", str(dev_stage))
    for coln in ElemeNT_out_df.columns:
        if ElemeNT_out_df[coln].isnull().values.all():
            print ("col name: %s is null" %(coln))
        else:
            sum_df.loc[sum_df['fileN']==basefn,coln] = len(ElemeNT_out_df[ElemeNT_out_df[coln]!= 'no'][coln])

    return sum_df

def summarize_Motifs_out(ElemeNT_out_col_l,ElemenNT_outFile,motifSum_df,normalize_motif_df):
    basefn = os.path.basename(ElemenNT_outFile)
    col_l = motifSum_df.columns.tolist()
    motifSum_df.loc[len(motifSum_df.index),'fileN'] = basefn
    dev_stage = basefn[5:9]
    motifSum_df.loc[motifSum_df['fileN'] == basefn, 'devStage'] = dev_stage
    motifSum_df.loc[motifSum_df['fileN'] == basefn, 'smooth'] = 0
    #print("time is: ", str(timeOfDev), "dev: ", str(dev_stage))
    for coln in ElemeNT_out_col_l:
        if coln =='name':
            motifSum_df.loc[motifSum_df['fileN'] == basefn, 'totalMotInFile'] = len(normalize_motif_df.index.values.tolist())

        else:
            len_df = len(normalize_motif_df[normalize_motif_df['motifName']== coln]['score'])
            motifSum_df.loc[motifSum_df['fileN']==basefn,coln] = len(normalize_motif_df[normalize_motif_df['motifName']== coln]['score'])

    return motifSum_df

def normalize_ElemeNT_out(ElemeNT_out_df):
    col_list = ['name','motifName','pos','score','InrScore','InrName']
    normalize_motif_df = pd.DataFrame(columns=col_list)
    for coln in ElemeNT_out_df.columns:
        if (coln == 'name' or  ElemeNT_out_df[coln].isnull().values.all()):
            continue
        motifValues_df = ElemeNT_out_df.loc[ElemeNT_out_df[coln]!='no']
        normalize_motif_df = getMotifPos(motifValues_df,coln,normalize_motif_df)

    return normalize_motif_df


def getMotifPos(motifValues_df,coln,normalized_motif_df):
    names_l = motifValues_df['name'].tolist()
    for name in names_l:
        check_dup = motifValues_df.loc[motifValues_df['name']==name,coln].tolist()
        if len(check_dup) == 1:
            motifPos = motifValues_df.loc[motifValues_df['name']==name,coln].item()
            normalized_motif_df = extractMotifPos(name,coln,motifPos,normalized_motif_df)
        else:
            print ("len of %s is not 1 len is: %d" %(name,len(check_dup)))

    return normalized_motif_df

def extractMotifPos(name,coln,motifPos,normalized_motif_df):
    if coln in ['MTE', 'DPE', 'bridge']:
        dependent_motif = True
        inr_dict = prep_inr_dict(name, coln, motifPos, ElemeNT_out_df)
    else:
        dependent_motif = False
        inr_dict ={}
    motifPos_l = motifPos.split(';')
    for pos in motifPos_l:
        if len(pos)==0:
            continue
        pos_score_l = pos.split(',')
        motifPos = pos_score_l[0]
        motif_score = pos_score_l[1]
        if dependent_motif:
            inrName_l = motif_score.split('(')
            inrName = inrName_l[1]
            inrName = inrName[:-1]
            inrScore_val = inr_dict[inrName][motifPos]
            motif_score = inrName_l[0]
        else:
            inrScore_val = ''
            inrName = ''
        normalized_motif_df.loc[len(normalized_motif_df.index)] = [name,coln,motifPos,motif_score,inrScore_val,inrName]
    return normalized_motif_df

def prep_inr_dict(name, coln, motifPos, ElemeNT_out_df):
    import collections
    dependentDict = collections.defaultdict(dict)
    motifPos_l = motifPos.split(';')
    for mPos in motifPos_l:
        if 'dInr' in mPos:
            motif_inr_pos = ElemeNT_out_df.loc[ElemeNT_out_df['name'] == name, 'dInr'].item()
            dependentDict = prep_motifPos_dict('dInr', dependentDict, motif_inr_pos)
        else:
            if 'hInr' in mPos:
                motif_inr_pos = ElemeNT_out_df.loc[ElemeNT_out_df['name'] == name, 'hInr'].item()
                dependentDict = prep_motifPos_dict('hInr', dependentDict, motif_inr_pos)
            else:
                if 'BBCABW' in mPos:
                    motif_inr_pos = ElemeNT_out_df.loc[ElemeNT_out_df['name'] == name, 'BBCABW'].item()
                    dependentDict = prep_motifPos_dict('BBCABW', dependentDict, motif_inr_pos)
    if len(motif_inr_pos) == 0:
        print('len of Inr for motif of gene name: %s motif_inr %s is 0 ' % (name, coln))
        exit()
    return dependentDict

def prep_motifPos_dict(motifName, dependentDict, motif_inr_pos):
    motif_inr_pos_l = motif_inr_pos.split(';')
    for inr_pos in motif_inr_pos_l:
        if len(inr_pos) == 0:
            continue
        inr_pos_l = inr_pos.split(',')
        pos = inr_pos_l[0]
        score = inr_pos_l[1]
        dependentDict[motifName][pos] = score

    return dependentDict

def save_summary_filesForSpecie(specie_name,pipelineOut_dir,sum_df,motifSum_df):
    sum_fname = specie_name + '_ElemeNTsum.csv'
    sum_fn = os.path.join(pipelineOut_dir,sum_fname)
    sum_df = sum_df.sort_values(by=['smooth', 'devStage'])
    sum_df.to_csv(sum_fn,index=False)

    sum_fname = specie_name + '_MotifSum.csv'
    sum_fn = os.path.join(pipelineOut_dir,sum_fname)
    motifSum_df = motifSum_df.sort_values(by=['smooth', 'devStage'])
    motifSum_df.to_csv(sum_fn,index=False)

    return

# This program 'normalizes' the ElemeNT output so that for every sequence, every motif in every position appears
# in a separate line
#

dir_mainPath = "//Main//Path"

motifsInfo_dir = os.path.join(dir_mainPath,'scripts//Motifs')
motifPosNlen_fn = os.path.join(motifsInfo_dir,'motifs_len.csv')
motifLen_df = pd.read_csv(motifPosNlen_fn,sep=',',header=0)

tss_start =100

specie_name = 'AllCoreSpecie'
core_dir = os.path.join(dir_mainPath,"CoreSeq/sequences052023")
specie_elemeNTout_dir = os.path.join(core_dir, 'ElemeNTout')

sum_df = create_summary_df(specie_name)
motifSum_df = create_summary_df(specie_name)

# For every specie and every dev stage reading the output of ElemeNT:
file_suffix = '*EPDnew_minPlus100__start' + str(tss_start) + 'smooth10.txt'
for ElemenNT_outFile in glob.glob(os.path.join(specie_elemeNTout_dir, file_suffix)):
    ElemenNT_out_fn = os.path.basename(ElemenNT_outFile)
    ElemenNTnormalized_out_fn = ElemenNT_out_fn[:-4] + '_norm.csv'
    print(ElemenNT_out_fn, ElemenNTnormalized_out_fn)
    ElemeNT_out_df = pd.read_csv(ElemenNT_outFile,sep='\t',header=0)
    print('shape',ElemeNT_out_df.shape)
    # sum_df includes a summary of how many lines in ElemeNT output included the motif (one time or more)
    sum_df = summarize_ElemeNT_out(ElemeNT_out_df,ElemenNT_outFile,sum_df)
    normalize_motif_df = normalize_ElemeNT_out(ElemeNT_out_df)
    targetElemeNT_dir = create_chippeak_dir(specie_elemeNTout_dir,'normlized')
    normalized_out_fn = os.path.join(specie_elemeNTout_dir,'normlized',ElemenNTnormalized_out_fn)
    normalize_motif_df.to_csv(normalized_out_fn,index=False)
    ElemeNT_out_col_l = ElemeNT_out_df.columns.tolist()
    # motif Sum includes a summary of how many times the motif was identified by ElemeNT in a specific ElemeNT output file
    motifSum_df = summarize_Motifs_out(ElemeNT_out_col_l,
                                       ElemenNT_outFile, motifSum_df, normalize_motif_df)
save_summary_filesForSpecie(specie_name,targetElemeNT_dir,sum_df,motifSum_df)








