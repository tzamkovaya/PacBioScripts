#!/usr/bin/python

import cx_Oracle
import pandas as pd
import sys,os
import argparse
import fnmatch
from subprocess import Popen, PIPE
import pathlib

parser = argparse.ArgumentParser(description='create target,backbone, repCap fasta, sampleTable config, and summary config')
parser.add_argument('filepath', help='directory of new PacBio AAV data')
args = parser.parse_args()

connection = cx_Oracle.connect(user=user,password=pwd,dsn=dsn,encoding='UTF-8')
cursor=connection.cursor()
filePath=args.filePath
contaminantsPath = 'references/contaminants'
vectorPath = 'references/vectors'

file_list=[]
for path, folders, files in os.walk(filePath):
	for file in files:
		if fnmatch.fnmatch(file, '*Samples-and-Sequences.xlsx'):
			file_list.append(os.path.join(path, file))

origFileName = file_list[0]
origDf = pd.read_excel(origFileName)
ppbIdList = origDf['PPB ID'].tolist()
ppbIdList2 = [x.split("-")[1] for x in ppbIdList]

queryAllAAV = """select ppb.VISIBLE_ID PPB_VIS_ID, v.id, v.visible_id VEC_VIS_ID, v.alias, ci.LABEL, gs.content, nuf.name, nuf.RANGE_START, nuf.RANGE_STOP
FROM VECTOR v,
    GEN_SEQUENCE gs,
    GEN_NU_FEATURE nuf,
	protein_expression_batch peb,
	peb_vector pv,
	peb_ppb ep, 
	protein_purification_batch ppb,
	peb_ppt,
	ppt, 
	ppt_vector,
	chain_info ci
where
    v.id = pv.VECTOR_ID and
	pv.PROTEIN_EXPRESSION_BATCH_ID = peb.id and
	ep.PEB_ID = peb.id and
	ep.PPB_ID = ppb.ID and
	peb_ppt.PEB_ID = peb.id and
	peb_ppt.PPT_ID = ppt.id and
	ppt_vector.ppt_id = ppt.id and
	ppt_vector.VECTOR_ID = v.id and
	ci.id = ppt_vector.CHAIN_INFO_ID and 
	ppb.VISIBLE_ID = :ppbId and 
	ci.LABEL = 'Payload DNA' and
	gs.GES_ID = v.GES_ID and
	nuf.nus_id = v.GES_ID and
	(nuf.name like '%ITR' or 
	lower(nuf.name) ='payload_dna')
"""

queryRepCap = """select ppb.VISIBLE_ID PPB_VIS_ID, v.id, v.visible_id VEC_VIS_ID, v.alias, ci.LABEL, gs.content
FROM VECTOR v,
    GEN_SEQUENCE gs,
    GEN_NU_FEATURE nuf,
	protein_expression_batch peb,
	peb_vector pv,
	peb_ppb ep, 
	protein_purification_batch ppb,
	peb_ppt,
	ppt, 
	ppt_vector,
	chain_info ci
where
    v.id = pv.VECTOR_ID and
	pv.PROTEIN_EXPRESSION_BATCH_ID = peb.id and
	ep.PEB_ID = peb.id and
	ep.PPB_ID = ppb.ID and
	peb_ppt.PEB_ID = peb.id and
	peb_ppt.PPT_ID = ppt.id and
	ppt_vector.ppt_id = ppt.id and
	ppt_vector.VECTOR_ID = v.id and
	ci.id = ppt_vector.CHAIN_INFO_ID and 
	ppb.VISIBLE_ID = :ppbId and 
	ci.LABEL = 'Payload DNA' and
	gs.GES_ID = v.GES_ID and
	nuf.nus_id = v.GES_ID and
	ci.LABEL = 'Rep'
"""
dfL = []
vecIdList = []
targetVecL = []
backboneL = []
repCapL = []
targetSizeL = []

for ppbId in sorted(set(ppbIdList2), key=ppbIdList2.index):
    print('ppbId:', ppbId)
    allAavDf = pd.read_sql(queryAllAAV, connection, params=[ppbId])
    vecId = allAavDf['VEC_VIS_ID'].unique().tolist()
    vecId =  ''.join(map(str,vecId))
    vecIdList.append(vecId)
    
    repCapDf = pd.read_sql(queryRepCap, connection, params=[ppbId])
    repName = repCapDf['VEC_VIS_ID'].values[0]
	### change name to PPBID_VEC_VIS_ID.fa
    repCapName = os.path.join(contaminantsPath, "%s_%s.fa") % (ppbId, repName)

    df = allAavDf
    df = df.loc[-(df.duplicated(subset=['NAME', 'RANGE_START', 'RANGE_STOP'])), :]
    name = df['VEC_VIS_ID'].values[0]
	### calculate total length of sequence
    df['SeqLength'] = df['CONTENT'].str.count("") 
	### use VEC ID as name of target fasta file
    fileName = os.path.join(vectorPath, "%s_%s.fasta") % (ppbId, name)
	### create backbone file
    backboneName = os.path.join(vectorPath, "%s_%s_Backbone.fasta") % (ppbId, name)
	
    vector = str(ppbId) + "_" + str(name)
    targetVecL.append(vector)
    pBackbone = vector + "_Backbone"
    backboneL.append(pBackbone)
    repcap = str(ppbId) + "_" + str(repName)
    repCapL.append(repcap)
    helper = 'pHelper'
    
    ### create vector/backbone fasta files of unique vec id & nuc seq (CONTENT)
    with open(fileName, 'w') as f:
	    with open(backboneName, 'w') as b:
            f.write('>' + str(vector))
			f.write('\n')
            b.write('>' + str(vector) + '_Backbone')
            b.write('\n')
            df['NAME'] = df['NAME'].str.lower()
            for val in df['NAME']:
                if 'payload' in val:
                    df2 = df.loc[df['NAME'].str.contains('payload')]
					s = df2['RANGE_START'].item()
					e = df2['RANGE_STOP'].item()
					seq = ''.join(str(x) for x in df2['CONTENT'])
					payloadSeq = seq[s-1:e]
                    backboneSeq = seq[e:-1] + seq[:s]
					
            targetSizeL.append(len(payloadSeq))
            b.write(str(backboneSeq) + '\n')
            b.close()
            f.write(str(payloadSeq) + '\n')
            f.close()
	
    ### create repCapFa file
	with open(repCapName, 'w') as r:
        seq = ''.join(str(x) for x in repCapDf['CONTENT'])
        r.write('>' + str(repcap))
        r.write('\n')
		r.write(str(seq) + '\n')
		r.close()
		
    os.system("sbatch reference_masking_template.sh {} {} {}".format(str(vector), str(helper), str(repcap)))
	

sampleNameL = origDf['Sample'].tolist()
bcL = origDf['PacBio barcode index'].tolist()

dataPathL = origDf['Data Path'].tolist()
dataPathL2 = [os.path.basename(os.path.dirname(x)) for x in dataPathL]
fileName = [os.path.basename(os.path.normpath(x)) for x in dataPathL]
fileForBc = [s.split('.', 1)[0] for s in fileName]

refFastaL = [s + '_master.fasta' for s in targetVecL]
refmmiL = [s + '_master_ccs.fasta.mmi' for s in targetVecL]
extAlign = ['out_aligned_stripped_' + sampleName + '_hifi.csv' for sampleName in sampleNameL]
covFile = [sampleName + '.var.vcf' for sampleName in sampleNameL]
smallVariantFile = [sampleName + '.csv' for sampleName in sampleNameL]

outDir = [os.path.join(filePath, dataPath, 'out', s) for dataPath in dataPathL2 for s in sampleNameL]
readPath = [os.path.join(filePath, dataPath, 'ccs', f + '.demux' + bc + '_BAK8A_OA--' + bc + '_BAK8A_OA.bam') for dataPath,f,bc in zip(dataPathL2, fileForBc, bcL)]
config = ['CCS_HIFI_in.txt']*len(sampleNameL)

parentPath = [os.path.join(filePath, dataPath, 'out', s) for dataPath in dataPathL2 for s in sampleNameL]

r_itr_start = [x- 144 for x in targetSizeL]
r_50itr_start = [x- 72 for x in targetSizeL]
truncation_R = [x- 154 for x in targetSizeL]

### create new sampleTable
sampleTableDf = pd.DataFrame(list(
	zip(sampleNameL, 
	 readPath, 
	 refFastaL, 
	 refmmiL,
	 config,
	 targetVecL, 
	 outDir, 
	 backboneL,
	 ['pHelper']*len(sampleNameL),
	 repCapL,
	 ['Homo_sapiens']*len(sampleNameL),
	 ['E_coli']*len(sampleNameL),
	 [145]*len(sampleNameL),
	 r_itr_start,
	 [0]*len(sampleNameL),
	 [0]*len(sampleNameL))
),
columns = [
	'Sample_ID',
	'Reads_path',
	'Ref_fasta',
    'Ref_mmi',
	'Config',
	'TargetName',
	'OutputDirectory',
	'pBackbone',
	'pHelper',
	'pRepCap',
	'Host',
	'E_coli',
	'L_ITR_end',
	'R_ITR_start',
	'Amplicon_L',
	'Amplicon_R'
]
)

sampleTableFile = os.path.join(filePath, 'sampleTableConfigFile.csv')
sampleTableDf.to_csv(sampleTableFile, index=False)
print(sampleTableDf)