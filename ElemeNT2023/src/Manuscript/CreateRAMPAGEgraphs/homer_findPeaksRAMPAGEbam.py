import glob
import os
import sys

def prepare_homerTagDir(BAMfile, dm6Genome_f):
    basefn = os.path.basename(BAMfile)
    tag_dirName = basefn[:-4] + '_l100/'
    makeTagDir_cmd = 'makeTagDirectory ' + tag_dirName + ' BAM/' + basefn + \
     ' -genome ' + dm6Genome_f + ' -checkGC -fragLength 100'
    # print(tag_dirName , '\n', makeTagDir_cmd)
    os.system(makeTagDir_cmd)
    return tag_dirName

def homer_findPeaks(tag_dirName,dm6Annotation_f,dm6Genome_f):
    
    findPeaks_cmd = 'findcsRNATSS.pl ' + tag_dirName + ' -o ' + 'findcsRNATSSoutput/tssAnalysis_' + tag_dirName[:-6] +\
                    ' -gtf ' +  dm6Annotation_f \
                    + ' -genome ' + dm6Genome_f

    print ( '==>' , findPeaks_cmd)
    os.system(findPeaks_cmd)
    return
print(sys.argv[0])
for i in range(1, len(sys.argv)):
    print('argument:', i, 'value:', sys.argv[i])
    if i== 1:
        dir_mainPath = sys.argv[i]
    else: # i=2
        dm6Genome_dir = sys.argv[i]

dm6Annotation_f = os.path.join(dm6Genome_dir,'Annotation//Genes','genes.gtf')
dm6Genome_f = os.path.join(dm6Genome_dir,'Sequence//WholeGenomeFasta','genome.fa')
pipelineOut_dir = os.path.join(dir_mainPath,"scripts//pipelineOut")

homer_dir = os.path.join(pipelineOut_dir,"GSE89299.dm6//HOMER")
BAM_dir = os.path.join(pipelineOut_dir,"GSE89299.dm6//HOMER//BAM")

dev_stage4analysis_l = ['01','02','03','04','05','06','07','08']

os.chdir(homer_dir)
for BAMfile in glob.glob(os.path.join(BAM_dir, 'GSM*.bam')):
    basefn = os.path.basename(BAMfile)
    dev_stage = basefn[-6:-4]
    if dev_stage in dev_stage4analysis_l:
        print(basefn, dev_stage)
        tag_dirName = prepare_homerTagDir(BAMfile, dm6Genome_f)
        homer_findPeaks(tag_dirName,dm6Annotation_f,dm6Genome_f)

