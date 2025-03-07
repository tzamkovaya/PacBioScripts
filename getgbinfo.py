#!/usr/bin/env python

import os,sys
import gbparse_noquotes
from Bio import SeqIO
import pandas as pd
import numpy as np
import json
import utilities
import re
import operator
from Bio import SeqFeature
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

def qualityScore(pString):
 
    pString = pString.strip()
 
    phred = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ"
   
    total = 0
    counts = 0
   
    for i in range(len(pString)):
        pIdx = phred.find(pString[i])
        total += pIdx
        counts += 1
       
    return round(float(total) / float(counts), 2)


def readpileupPath(pileupPath):
    dataList = []
    lines = utilities.readlines(pileupPath)

    #print(len(lines))
    
    #new_consensus_seq = []
    for line in lines:
        items = line.split("\t")
        #print(items)
        #print(len(items))
        refPos = items[1]
        #print(refPos)
        refAa = items[2]
        #print(refAa)
        coverage = items[3]
        quality = qualityScore(items[5])
        #print("Coverage is: " + coverage)
        #print(quality)
        aaStr = items[4]
        #print('aaStr:',aaStr)
        aaCounts = {}
        for base in ["A", "C", "T", "G"]:
            if base != "G":
                aaCounts[base] = aaStr.count(base) + aaStr.count(base.lower())
            else:
                #aaCounts[base] = aaStr.count(base.lower()) - for Kristine AAV data
                aaCounts[base] = aaStr.count(base.lower())
        g_indices = [m.start() for m in re.finditer("[G]", aaStr)]
        if len(g_indices) != 0:
           for g_pos in g_indices:
               if aaStr[int(g_pos)-1] != "^":
                   aaCounts["G"] += 1
        aaCounts[refAa] = aaStr.count(".") + aaStr.count(",")
        #print(aaCounts[refAa])
        ins_indices = [m.start() for m in re.finditer("[+d]", aaStr)]
        ins_list = []
        if len(ins_indices) != 0:
            for ins_pos in ins_indices:
                if aaStr[int(ins_pos)+1].isdigit():
                    ins_list.append(int(aaStr[int(ins_pos)+1]))
        aaCounts["I"] = sum(ins_list)
        del_indices = [m.span() for m in re.finditer("-[0-9]+[ACGTNacgtn]+", aaStr)]
        del_list = []
        num_C_to_del = []
        num_A_to_del = []
        num_G_to_del = []
        num_T_to_del = []
        if len(del_indices) != 0:
            for del_pos in del_indices:
                if aaStr[int(del_pos[0])-1] != "^":
                    delStr = aaStr[del_pos[0]: del_pos[1]]
                    dList = []
                    for d in delStr:
                        if d.isdigit():
                            dList.append(d)
                    allD = ''.join(dList)
                    if delStr.count('C') != 0:
                        nc = delStr.count('C')
                        num_C_to_del.append(int(nc))
                    if delStr.count('G') != 0:
                        ng = delStr.count('G')
                        num_G_to_del.append(int(ng))
                    if delStr.count('A') != 0:
                        na = delStr.count('A')
                        num_A_to_del.append(int(na))
                    if delStr.count('T') != 0:
                        nt = delStr.count('T')
                        num_T_to_del.append(int(nt))
                    del_list.append(int(allD))
                    #del_list.append(int(aaStr[int(del_pos)+1]))
        aaCounts["D"] = sum(del_list)
        #print('C:', aaCounts["C"] , 'A:', aaCounts["A"], 'G:', aaCounts["G"], 'T:', aaCounts["T"])
        #print('total to delete:', sum(del_list))
        if len(num_C_to_del) != 0:
            totaldelC = sum(num_C_to_del)
            aaCounts["C"] = aaCounts["C"] - totaldelC
        if len(num_A_to_del) != 0:
            totaldelA = sum(num_A_to_del)
            aaCounts["A"] = aaCounts["A"] - totaldelA
        if len(num_G_to_del) != 0:
            totaldelG = sum(num_G_to_del)
            aaCounts["G"] = aaCounts["G"] - totaldelG
        if len(num_T_to_del) != 0:
            totaldelT = sum(num_T_to_del)
            aaCounts["T"] = aaCounts["T"] - totaldelT
        #print('newC:', aaCounts["C"] , 'newA:', aaCounts["A"], 'newG:', aaCounts["G"], 'newT:', aaCounts["T"])
        #print(aaCounts)
        #keymax = max(aaCounts.items(), key = operator.itemgetter(1))[0]
        #new_consensus_seq.append(keymax) 
        d = {"label": refPos, "refBase": refAa, "coverage": coverage, "quality": quality}

        for aa in aaCounts: d[aa] = aaCounts[aa]
        dataList.append(d)
    
    for i in range(len(dataList)):
        #print(dataList[i])
        refB = dataList[i]['refBase']
        #print("ref B is", refB)
        refBcov = dataList[i][refB]
        cov = dataList[i]['coverage']
    jsonStr = json.dumps(dataList)
    print(jsonStr)
    #print("".join(new_consensus_seq))
    return dataList
    
def jsontocsv(datajsonPath, outcsvPath):
    with open(datajsonPath, 'r') as datafile:
        data = json.load(datafile)
        #data = json.loads(json.dumps(datafile))
        df = pd.DataFrame(data)
        df.to_csv(outcsvPath)

def getConsSeqAS(newSeq, origGbPath, consGbPath):
    sec_rec2L = []
    for seq_rec in SeqIO.parse(newSeq, "genbank"):
        print('orig:',seq_rec.seq) 
        #sec_rec2 = seq_rec.seq.upper()
        #print('upper:',sec_rec2)
        sec_rec2L.append(str(seq_rec.seq))

    fullseq = "".join(sec_rec2L)
    fullseq2 = Seq(fullseq, generic_dna)
    print('this is fullseq:', fullseq)
    newSeqRecords = []
    for seq_record in SeqIO.parse(origGbPath, "genbank"):
        #print(seq_record.id)
        #print(seq_record.seq)
        #print(len(seq_record))
        #print(newSeq)
        seq_record.seq = fullseq2
        newSeqRecords.append(seq_record)
    SeqIO.write(newSeqRecords,consGbPath, "genbank")


def getConsSeq(pileupPath, origGbPath, consGbPath):
    dataList = []
    lines = utilities.readlines(pileupPath)
    first_line = lines[0]
    first_items = first_line.split("\t")
    first_ref_pos = first_items[1]
    #print(first_line)
    #print(first_ref_pos)
    new_consensus_seq = []
    for line in lines:
        items = line.split("\t") 

        refPos = items[1]
        #print(refPos)
        refAa = items[2]
        #print(refAa)
        coverage = items[3]
        quality = qualityScore(items[5])
        #print("Coverage is: " + coverage)
        #print(quality)
        aaStr = items[4]
        
        aaCounts = {}
        for base in ["A", "C", "T", "G"]:
            if base != "G":
                aaCounts[base] = aaStr.count(base) + aaStr.count(base.lower())
            else:
                aaCounts[base] = aaStr.count(base.lower()) 
        g_indices = [m.start() for m in re.finditer("[G]", aaStr)]
        if len(g_indices) != 0:
           for g_pos in g_indices:
               if aaStr[int(g_pos)-1] != "^":
                   aaCounts["G"] += 1       
        aaCounts[refAa] = aaStr.count(".") + aaStr.count(",")
        ins_indices = [m.start() for m in re.finditer("[+]", aaStr)]
        ins_list = []
        if len(ins_indices) != 0: 
            for ins_pos in ins_indices:
                ins_list.append(int(aaStr[int(ins_pos)+1]))
        aaCounts["I"] = sum(ins_list)
        del_indices = [m.start() for m in re.finditer("[-]", aaStr)]
        del_list = []     
        if len(del_indices) != 0:
            for del_pos in del_indices:
                if aaStr[int(del_pos)-1] != "^":
                    del_list.append(int(aaStr[int(del_pos)+1]))
        aaCounts["D"] = sum(del_list)
        #print(aaCounts)
        if(all(x==0 for x in aaCounts.values()) == True):
            #print("true")
            keymax = refAa
        else: 
            keymax = max(aaCounts.items(), key = operator.itemgetter(1))[0]
        
        new_consensus_seq.append(keymax)
        d = {"label": refPos, "refBase": refAa, "coverage": coverage, "quality": quality}

        for aa in aaCounts: d[aa] = aaCounts[aa]
        dataList.append(d)
    #print(dataList)
    new_cons_seq = ''.join(new_consensus_seq)
    #print(new_cons_seq)
    if int(first_ref_pos) > 1:
        amb_pos = 'N'
        num_amb = int(first_ref_pos) - 1
        n_to_add = amb_pos * num_amb 
        #print(len(n_to_add))
        new_cons_seq = n_to_add + new_cons_seq 
    else:
        new_cons_seq = new_cons_seq
    #print(new_cons_seq)
    new_cons_seq2 = Seq(new_cons_seq, generic_dna)
    #kpath = "/hpc/grid/wip_cmg_gbt/replication/users/debarj02/NGS_rawdata/kristine/phGBA22_pWTMYBPC3spacer_pWTMYBPC3/Sequences/p-cTNT-hMYBPC3-sPolyA-Kan.gb"
    newSeqRecords = []
    for seq_record in SeqIO.parse(origGbPath, "genbank"):
        #print(seq_record.id)
        #print(seq_record.seq)
        #print(len(seq_record))
        
        newSeq = new_cons_seq2
        #print(newSeq)
        seq_record.seq = newSeq
        newSeqRecords.append(seq_record)
 
        #start_pos = int(first_ref_pos)
        #end_pos = len(newSeq)
        #my_location = SeqFeature.FeatureLocation(start_pos, end_pos)

        #my_feat = SeqFeature(my_location, type = "misc_feature")
        #newSeqRecords.features.append(my_feat)
        #end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
        #my_location = SeqFeature.FeatureLocation(start_pos, end_pos)
 
    SeqIO.write(newSeqRecords,consGbPath, "genbank")

    #return dataList


dm = {"annotations": [],"data": []}

def gethg38info(pileupPath,countsPath):

    dm["data"] = readpileupPath(pileupPath)

    annot_counts = []

    test_pileup_data = dm["data"]
    tdf = pd.DataFrame(test_pileup_data, columns = ['label', 'refBase', 'coverage', 'quality', 'A', 'C', 'T', 'G', 'I', 'D'])
    tdf['coverage'] = pd.to_numeric(tdf['coverage'])
    total_coverage = sum(tdf.coverage)
    tdf['label'] = pd.to_numeric(tdf['label'])
    print("total coverage:", total_coverage)
    d = { "total_coverage": total_coverage}
    annot_counts.append(d)

    print(annot_counts)
    annot_counts_df = pd.DataFrame(annot_counts, columns = ['total_coverage'])
    annot_counts_df.to_csv(countsPath)


                

 
dm = {"annotations": [], "data": []}

def getjsoninfo(gbPath,pileupPath,countsPath,outputPath):
    annotations = gbparse_noquotes.getGenbankAnnotations(gbPath)

    dm["annotations"] = annotations
    dm["data"] = readpileupPath(pileupPath)

    #print(dm["annotations"])
    annot_counts = []
    
    for i in range(len(annotations)):
        #print(annotations[i])
        start_pos = int(annotations[i]['start'])
        #print(start_pos)
        end_pos = int(annotations[i]['end'])
        annot = annotations[i]['label']
        test_pileup_data = dm["data"]
        tdf = pd.DataFrame(test_pileup_data, columns = ['label', 'refBase', 'coverage', 'quality', 'A', 'C', 'T', 'G', 'I', 'D'])
        tdf['coverage'] = pd.to_numeric(tdf['coverage'])
        total_coverage = sum(tdf.coverage)
        tdf['label'] = pd.to_numeric(tdf['label'])
        tdf2 = tdf[(tdf.label >= start_pos) & (tdf.label <= end_pos)]
        annot_coverage = sum(tdf2.coverage)
        percent_mapped = annot_coverage / total_coverage * 100
        
        d = {"annotation": annot, "total_coverage": total_coverage, "annotation_coverage": annot_coverage, "percent_mapped": percent_mapped, "start": start_pos, "end": end_pos}
        annot_counts.append(d)
        
    print(annot_counts)
    dm["counts"] = annot_counts

    annot_counts_df = pd.DataFrame(annot_counts, columns = ['annotation', 'total_coverage', 'annotation_coverage', 'percent_mapped', 'start', 'end'])
    annot_counts_df[['percent_mapped', 'start', 'end']] = annot_counts_df[['percent_mapped', 'start', 'end']].apply(pd.to_numeric)
    annot_counts_df.to_csv(countsPath)
    
    utilities.writeJson(dm, outputPath) 
 

def getcountsAllPlasmids(helperCountsPath, repcapCountsPath, plasCountsPath, outPath):
    helper = pd.read_csv(helperCountsPath)
    repcap = pd.read_csv(repcapCountsPath)
    plas = pd.read_csv(plasCountsPath)
 
    total_cov_h = helper['total_coverage'].unique()
    total_cov_repcap = repcap['total_coverage'].unique()
    total_cov_p = plas['total_coverage'].unique()
    total_cov = int(total_cov_h) + int(total_cov_repcap) + int(total_cov_p)
    
    helper['new_percent_mapped'] = helper['annotation_coverage'] / total_cov * 100
    sum_helper_cov = sum(helper['new_percent_mapped'])
    repcap['new_percent_mapped'] = repcap['annotation_coverage'] / total_cov * 100
    sum_repcap_cov = sum(repcap['new_percent_mapped'])
    start_a = plas.loc[plas['annotation'] == '5 pri ITR', 'start'].iloc[0]
    end_a = plas.loc[plas['annotation'] == '3 pri ITR', 'end'].iloc[0]
    vec = plas[(plas.start >= start_a) & (plas.end <= end_a)]
    vec['new_percent_mapped'] = vec['annotation_coverage'] / total_cov * 100
    sum_vec = sum(vec['new_percent_mapped'])
    sum_bac = 100 - sum_vec - sum_helper_cov - sum_repcap_cov
    data = [['Helper', sum_helper_cov], ['RepCap', sum_repcap_cov], ['Vector', sum_vec], ['Plasmid', sum_bac]]
    df = pd.DataFrame(data, columns = ['refGenome', 'percMapping'])
    #print(df)
    #df.to_csv(outPath)
    full_plas_mapping = df.to_json(outPath,orient="records")
    print(full_plas_mapping)
    return(full_plas_mapping)
   # list_dict = []
   # for index,row in list(df.iterrows()):
   #     list_dict.append(dict(row))

   # with open(outPath, 'w') as f:
   #     f.write("\n".join(str(item) for item in list_dict))
  
def getAlignmentNum(bowtieLogPath):
    with open(bowtieLogPath) as f:
        for line in f:
            pass
        last_line = line

        num = re.findall( r'\d+\.*\d*', last_line)
        if len(num) == 1:
            num = float(num[0]) / 100 
    return(num) 


def getcountsAllPlasHG38(helperCountsPath, helperBowtie,repcapCountsPath, repcapBowtie,plasCountsPath, plasBowtie, hgCountsPath, hgBowtie,outPath):
    helper = pd.read_csv(helperCountsPath)
    repcap = pd.read_csv(repcapCountsPath)
    plas = pd.read_csv(plasCountsPath)
    hg38 = pd.read_csv(hgCountsPath)

    start_a = plas.loc[plas['annotation'] == '5 pri ITR', 'start'].iloc[0]
    end_a = plas.loc[plas['annotation'] == '3 pri ITR', 'end'].iloc[0]
    total_cov_h = helper['total_coverage'].unique()
    total_cov_repcap = repcap['total_coverage'].unique()
    total_cov_p = plas['total_coverage'].unique()
    total_cov_hg38 = hg38['total_coverage'].unique()
    
    helperNum = getAlignmentNum(helperBowtie)
    hgNum = getAlignmentNum(hgBowtie)
    repcapNum = getAlignmentNum(repcapBowtie)
    plasNum = getAlignmentNum(plasBowtie)
 
    total_cov_helper = int(total_cov_h) * helperNum
    total_cov_rep = int(total_cov_repcap) * repcapNum
    total_cov_plas = int(total_cov_p) * plasNum 
    total_cov_hg = int(total_cov_hg38) * hgNum
    
    #new way
    vec = plas[(plas.start >= start_a) & (plas.end <= end_a)]
    vec['new_percent_mapped'] = vec['annotation_coverage'] / total_cov_p * 100
    total_vec_cov = sum(vec['annotation_coverage']) * plasNum
    total_bac_cov = total_cov_p - (sum(vec['annotation_coverage'])) * plasNum

    total_cov = total_cov_helper + total_cov_rep + total_cov_plas + total_cov_hg
    #total_cov = total_cov_helper + total_cov_rep + int(total_vec_cov) + int(total_bac_cov)
    print(total_cov)
    sum_helper_cov = total_cov_helper / total_cov * 100
    sum_repcap_cov = total_cov_rep / total_cov * 100
    sum_hg38_cov = total_cov_hg / total_cov * 100
  
    vec = plas[(plas.start >= start_a) & (plas.end <= end_a)]
    vec['new_percent_mapped'] = vec['annotation_coverage'] 
    sum_vec = sum(vec['new_percent_mapped']) / total_cov * 100 
    sum_vec = sum_vec - sum_helper_cov - sum_repcap_cov - sum_hg38_cov
    sum_bac = 100 - sum_vec - sum_helper_cov - sum_repcap_cov - sum_hg38_cov
    
    data = [['Helper', sum_helper_cov], ['RepCap', sum_repcap_cov], ['HG38', sum_hg38_cov],['Vector', sum_vec], ['Plasmid', sum_bac]]
    df = pd.DataFrame(data, columns = ['refGenome', 'percMapping'])
    #print(df)
    #df.to_csv(outPath)
    full_plas_mapping = df.to_json(outPath,orient="records")
    #print(full_plas_mapping)
    return(full_plas_mapping)
   # list_dict = []
   # for index,row in list(df.iterrows()):
   #     list_dict.append(dict(row))

   # with open(outPath, 'w') as f:
   #     f.write("\n".join(str(item) for item in list_dict))



def getcountsVP(plasCountsPath, outPath):
    plas = pd.read_csv(plasCountsPath)

    start_a = plas.loc[plas['annotation'] == '5 pri ITR', 'start'].iloc[0]
    end_a = plas.loc[plas['annotation'] == '3 pri ITR', 'end'].iloc[0]
    vec = plas[(plas.start >= start_a) & (plas.end <= end_a)]
    total_cov_p = plas['total_coverage'].unique()
    vec['new_percent_mapped'] = vec['annotation_coverage'] / total_cov_p * 100
    sum_vec = sum(vec['new_percent_mapped'])
    sum_bac = 100 - sum_vec
    data = [['Helper', 0.00], ['RepCap', 0.00], ['Vector', sum_vec], ['Plasmid', sum_bac]]
    df = pd.DataFrame(data, columns = ['refGenome', 'percMapping'])
    full_plas_mapping = df.to_json(outPath, orient="records")
    

def getcountsVPhg38(plasCountsPath,plasBowtie,hgCountsPath,hgBowtie,outPath):
    plas = pd.read_csv(plasCountsPath)

    start_a = plas.loc[plas['annotation'] == '5 pri ITR', 'start'].iloc[0]
    end_a = plas.loc[plas['annotation'] == '3 pri ITR', 'end'].iloc[0]
    vec = plas[(plas.start >= start_a) & (plas.end <= end_a)]
    total_cov_p = plas['total_coverage'].unique()
    vec['new_percent_mapped'] = vec['annotation_coverage'] / total_cov_p * 100
    vec_cov = sum(vec['new_percent_mapped'])
    print("total vec % mapped:", vec_cov)
    bac_cov = 100 - vec_cov
    print("total bac % mapped:", bac_cov)
    plasNum = getAlignmentNum(plasBowtie)
    #vec['new_percent_mapped'] = vec['annotation_coverage'] / total_cov_p * 100
    total_vec_cov = sum(vec['annotation_coverage']) * plasNum
    total_bac_cov = total_cov_p - (sum(vec['annotation_coverage'])) * plasNum
    
    hg38 = pd.read_csv(hgCountsPath)
    total_cov_hg38 = hg38['total_coverage'].unique()
    print('total hg38 % mapped:', total_cov_hg38)
    hgNum = getAlignmentNum(hgBowtie)
#   plasNum = getAlignmentNum(plasBowtie)
    #total_cov_plas = int(total_cov_p) * plasNum
    total_cov_hg = int(total_cov_hg38) * hgNum
    total_cov = int(total_vec_cov) + int(total_bac_cov) + int(total_cov_hg)
    sum_hg38_cov = total_cov_hg / total_cov * 100
    #vec['new_percent_mapped'] = vec['annotation_coverage'] / total_cov * 100
    sum_vec = total_vec_cov / total_cov * 100
    sum_bac = int(total_bac_cov) / total_cov * 100

    data = [['Helper', 0.00], ['RepCap', 0.00],['HG38', sum_hg38_cov], ['Vector', sum_vec], ['Plasmid', sum_bac]]
    df = pd.DataFrame(data, columns = ['refGenome', 'percMapping'])
    full_plas_mapping = df.to_json(outPath, orient="records")


 
def createfinaljson(jsoninfoPath, plasmidcountsPath, outPath):
    with open(jsoninfoPath) as f:
        dm = json.load(f)
    with open(plasmidcountsPath) as f2:
        plas_counts = json.load(f2)
    dm["all_plas_counts"] = plas_counts
    
    with open(outPath, 'w') as outfile:
        json.dump(dm, outfile)

def addConfigJson(configJson,finalJson):
    with open(configJson) as f:
        dm = json.load(f)
    with open("/hpc/grid/wip_cmg_gbt/workspace/genedata/igSeqPipeline/ampliconAlignments/allDatasets.json") as origFile1:
        origInfo1 = json.load(origFile1)
        l1 = list(origInfo1.values())[0]
        l1.append(dm)
        data = {"data": l1}
    with open(finalJson, 'w') as out:
        json.dump(data, out)

    
args = sys.argv
if len(args) > 1:
    if args[1] == "-addConfigJson":
        addConfigJson(args[2], args[3])
    if args[1] == "-readpileupPath":
        readpileupPath(args[2])
    if args[1] == "-getcountsVP":
        getcountsVP(args[2], args[3])
    if args[1] == "-getcountsVPhg38":
        getcountsVPhg38(args[2], args[3], args[4], args[5], args[6])
    if args[1] == "-gethg38info":
        gethg38info(args[2], args[3])
    if args[1] == "-getcountsAllPlasHG38":
        getcountsAllPlasHG38(args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10])
    if args[1] == "-jsontocsv":
        jsontocsv(args[2], args[3])
    if args[1] == "-getConsSeq":
        getConsSeq(args[2], args[3], args[4])
    if args[1] == "-getConsSeqAS":
        getConsSeqAS(args[2], args[3], args[4])
    if args[1] == "-getjsoninfo":
        getjsoninfo(args[2], args[3], args[4], args[5])
    if args[1] == "-getcountsAllPlasmids":
        getcountsAllPlasmids(args[2], args[3], args[4], args[5])
    if args[1] == "-createfinaljson":
        createfinaljson(args[2], args[3], args[4])
