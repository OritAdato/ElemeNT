import glob
import os
import pandas as pd
import sys

print(sys.argv[0])
for i in range(1, len(sys.argv)):
    print('argument:', i, 'value:', sys.argv[i])
    if i== 1:
        dir_mainPath = sys.argv[i]

pipelineOut_dir = os.path.join(dir_mainPath,"scripts//pipelineOut")

srrNames_matrix_dir = os.path.join(dir_mainPath,"scripts//matrix_files")
dm6SrrNames_matrix_fn = os.path.join(srrNames_matrix_dir,"samples_GSE89299")

homerBAM_dir = os.path.join(pipelineOut_dir,"GSE89299.dm6//HOMER//BAM")

srrMatrix_colNames = ['GSMdevStage', 'GSMonly', 'SRX', 'SRR', 'RAMPAGE','details']
dm6SrrNames_df = pd.read_csv(dm6SrrNames_matrix_fn,sep='\t',  header=None, index_col=False,names=srrMatrix_colNames,)

os.chdir(homerBAM_dir)
for SRRfile in glob.glob(os.path.join(homerBAM_dir, 'SRR*.bam')):
    basefn = os.path.basename(SRRfile)
    SRRid = basefn[:basefn.find(".")]
    GSMdevStage_id = dm6SrrNames_df.loc[dm6SrrNames_df["SRR"]==SRRid,'GSMdevStage'].item()
    print (SRRid,GSMdevStage_id)
    newBam_name = GSMdevStage_id + '.bam'
    renam_fn = ' mv ' + basefn + ' ' + newBam_name
    print(renam_fn)
    os.system(renam_fn)