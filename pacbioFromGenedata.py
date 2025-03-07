#!/usr/bin/python

import cx_Oracle
import pandas as pd
import sys,os
import argparse
import re
import glob
import difflib
import fnmatch
import subprocess

parser = argparse.ArgumentParser(description='get back fasta files for vectors')
parser.add_argument('filepath', help='directory of PacBio data')
args = parser.parse_args()

connection = cx_Oracle.connect(user=user,password=pwd,dsn=dsn,encoding='UTF-8')
cursor=connection.cursor()
filePath=args.filePath
print('filePath:', filePath)

def readFasta(fileName):
	defline = ''
	sequence = ''
	fName = open(fileName)
	for line in fName:
		if line[0] == '>':
			defline = line.strip()
		else:
			sequence += line.strip()
	return(sequence)

file_list=[]
for path, folders, files in os.walk(filePath):
	for file in files:
		if fnmatch.fnmatch(file, '*Samples-and-Sequences.xlsx'):
			file_list.append(os.path.join(path, file))

origFileName = file_list[0]
origDf = pd.read_excel(origFileName)
vecIdList = origDf['Vector ID'].tolist()

dataL = origDf['Data Path'].tolist()
dataL2 = [s.replace('subreadset.xml', 'subreads.bam') for s in dataL]

bcL = origDf['PacBio Barcode Index'].tolist()
sampleL = origDf['Sample'].tolist()

vecQuery = """SELECT v.id, v.visible_id VEC_VIS_ID, v.ges_id VEC_GES_ID, v.alias, gs.content, nuf.name, nuf.RANGE_START, nuf.RANGE_STOP 
FROM VECTOR v, GEN_SEQUENCE gs, GEN_NU_FEATURE nuf
WHERE
v.visible_id = :vecId and
v.ges_id = gs.ges_id and
nuf.nus_id = v.ges_id
"""

vecFaQuery = """SELECT v.id, v.visible_id VEC_VIS_ID, v.ges_id VEC_GES_ID, v.alias, gs.content
FROM VECTOR v
JOIN GEN_SEQUENCE gs
ON v.ges_id = gs.ges_id
WHERE v.visible_id = :vecId
"""

### get back Genbank annotations and fasta sequence of each vector from Genedata database ###
for vecId in sorted(set(vecIdList), key=vecIdList.index):
	print('vecId:', vecId)
	testVecDf = pd.read_sql(vecQuery, connection, params=[vecId])
	vecFaName = os.path.join(filePath, 'VEC-', + "%s.fasta") % (vecId)
	vecDfName = os.path.join(filePath, '%s_vecDf.csv') % (vecId)
	testVecDf.to_csv(vecDfName)
	vecFaDf = pd.read_sql(vecFaQuery, connection, params=[vecId])
	vecFaDf['ALIAS'] = vecFaDf['ALIAS'].str.replace(" ", "_")

	with open(vecFaName, 'w') as f:
		f.write('>' + str(vecFaDf['VEC_VIS_ID'][0] + '_' + str(vecFaDf['ALIAS'][0]) + '\n' + str(vecFaDf['CONTENT'][0]) + '\n'))

### align reads, find structural variants, and create consensus  ###
for subreadsXML, bc, vec, sample in zip(dataL, bcL, vecIdList, sampleL):
	print(subreadsXML, bc, vec, sample)
	bamFile = subreadsXML.split('.')[0] + '.reads.bam'
	bamFile = os.path.basename(bamFile)
	demuxBam = bamFile.split('.')[0] + '.demux.' + bc + '_BAK8A_OA--' + bc + '_BAK8A_OA.bam'
	vecFasta = 'VEC-' + str(vec) + '.fasta'

	### generate index using pbmm2 ###
	os.system("ml smrtlink/10.1")
	os.system("pbmm2 index {} {}.mmi".format(vecFasta, str(vec)))
	
	### align using pbmm2 ###
	alignFile = str(vec) + '_' + str(bc)
	alignBam = alignFile + '.aligned.movie.bam'

	os.system("pbmm2 align {} {} --sort -j 8 -J 8".format(vecFasta, demuxBam, alignBam))
	print("alignment done")

	### find structural variants using freebayes, vcfstats ###
	fVar = alignFile + '.freebayes.var.vcf'
	fVarOut = alignFile + '.freebayes.stats.out'
	proc = subprocess.Popen(args['freebayes', '-f', vecFasta, alignBam], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = proc.communicate()
	
	with open(fVar, 'w') as f1:
		print(out.decode('utf-8'), file=f1)
	f1.close()

	proc2 = subprocess.Popen(args['vcfstats', fVar], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out2,err2 = proc2.communicate()
	
	with open(fVarOut, 'w') as f2:
		print(out2.decode('utf-8'), file=f2)
	f2.close()

	### create consensus using samtools consensus ###
	alignCons = alignFile + 'aligned.consensus.fasta'
	os.system('samtools-1.16/bin/samtools consensus {} --show-del yes --show-ins yes -o {} -@ 40'.format(alignBam, alignCons))
	print('consensus done')

	
	### pbsv variant analysis ###
	var1 = alignFile + '.svsig.gz'
	os.system('pbsv discover {} {}'.format(alignBam, var1))
	var1Out = alignFile + '.var.vcf'
	os.system('pbsv call --ccs -S 0 -m 10 --filter-near-contig-end 5 --cluster-max-ref-pos-diff 5 -x 500 {} {} {}'.format(vecFasta, var1,var1Out))	

	### mpileup for coverage info per read 
	mpileupFile = alignFile + '_pileup.out'
	### assume using sortedBam as alignBam
	proc = subprocess.Popen(args['samtools-1.16/bin/samtools', 'mpileup', '-f', vecFasta, alignBam], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = proc.communicate()
	
	with open(mpileupFile, 'w') as m1:
		print(out.decode('utf-8'), file=m1)
	m1.close()