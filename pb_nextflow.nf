#!/usr/bin/env nextflow

/*
*load modules
*ml medsci smrtlink/10.1
*ml medsci anaconda3 oracle_instant_client/19.6
*/

params.input_bam='hifi_reads.bam'
params.excel='Samples-and-Sequences.xlsx'
params.vecOnly=''
params.diffOrient = true
params.align_from_fq = true
params.finalName = 'finalResults'

workflow {
    bam_ch = filterBamByRq(params.input_bam)
    convertToFq(filterBamByRq.out.filt.bam)
    makeDemuxSheet(params.excel)
    demux_ch = demuxFq(convertToFq.out.filt_fastq, makeDemuxSheet.out.demux_sample_sheet)
    (fasta_ch, gb_ch) = getVecInfo(params.excel)

    demux_ch.flatten().map { filepath -> [filepath.baseName.replace("_", "-"), filepath]}.set { demuxCh_tuple}
    fasta_ch.flatten().map { filepath -> [filepath.baseName.replace("_", "-"), filepath]}.set { fastaCh_tuple}
    gb_ch.flatten().map { filepath -> [filepath.baseName.replace("_", "-"), filepath]}.set { gbCh_tuple}

    fastaCh_tuple.join(demuxCh_tuple).set { joined_ch}
    align_ch  = alignFilesForCons(joined_ch)
    align_ch.collect().set{ cons_ch }

    getDistPlots(params.excel, params.finalName, fasta_ch, demux_ch, cons_ch, gb_ch)
    
    }

process filterBamByRq {
    cpus 24

    input:
    val params.input_bam

    output:
    path "filt_rq999.bam", emit: filt_bam

    shell:
    '''
    inBam="!{params.input_bam}"

    bamtools filter -in !{params.input_bam} -out "filt_rq999.bam" -tag "rq":">=0.999"
    '''
}

process convertToFq {
    input:
    val filt_bam

    output:
    path "filt_rq999.fastq", emit: filt_fastq

    shell:
    '''
    i="!{filt_bam}"
    samtools fastq !{filt_bam} > "filt_rq999.fastq"
    '''
}

process makeDemuxSheet {
    input:
    val params.excel

    output:
    path "*.tsv", emit: demux_sample_sheet

    shell:
    '''
    python makeDemuxSampleSheet.py !{params.excel}
    '''

}

process demuxFq {
    input:
    path filt_fastq
    val demux_sample_sheet

    output:
    path "*fastq", emit: vec_vb_fastq

    shell:
    if(params.diffOrient)
        '''
        python demuxSamples_nf.py -f !{filt_fastq} -b !{demux_sample_sheet}
        '''
    else
        '''
        python demuxSamples1_nf.py -f !{filt_fastq} -b !{demux_sample_sheet}
        '''
}

process getVecInfo {
    input:
    val params.excel
    
    output:
    path "*.fasta", emit:vec_fasta
    path "*_vecDf.csv", emit: vec_gb

    shell:
    if(params.vecOnly)
        '''
        python getVecInfo.py !{params.excel}
        '''
    else
        '''
        python getVecInfoFromVb.py !{params.excel}
        '''
}

process alignFilesForCons {
    memory '48 GB'
    cpus '15'
    executor 'slurm'
    time '4h'

    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample), path(vec_fasta), val(vec_vb_fastq)

    output:
    path "*_aligned.sorted.consensus.fasta", emit: cons

    shell:
    if(params.align_from_fq)

        '''
        ml medsci smrtlink/10.1

        origFile="!{vec_fasta}"
        newAlignSam=${origFile/.fasta/_aligned.sam}
        newAlignBam=${newAlignSam/.sam/.bam}
        echo ${newAlignBam}
        finalBam=${newAlignBam/.bam/.sorted.bam}
        consFile=${finalBam/.bam/.consensus.fasta}

        minimap2 -ax map-pb !{vec_fasta} | !{vec_vb_fastq} > ${newAlignSam}
        samtools view -bS ${newAlignSam} | samtools sort -o ${finalBam}

        ml medsci anaconda3
        samtools-1.16/bin/samtools consensus ${finalBam} --show-del yes --show-ins yes -o ${consFile} -@ 15
        '''
    
    else
        '''
        ml medsci smrtlink/10.1
        origFile="!{vec_fasta}"
        newAlignBam=${origFile/.fasta/_aligned.bam}
        echo ${newAlignBam}
        finalBam=${newAlignBam/.bam/.sorted.bam}
        consFile=${finalBam/.bam/.consensus.fasta}

        pbmm2 align !{vec_fasta} !{filt_bam} ${newAlignBam} --sort -j 8 -J 8
        samtools sort -o ${finalBam} ${newAlignBam}

        ml medsci anaconda3
        samtools-1.16/bin/samtools consensus ${finalBam} --show-del yes --show-ins yes -o ${consFile} -@ 15
        '''
}

process getDistPlots {
    publishDir 'finalMutPlotResults/', mode:'copy'

    input:
    val params.excel
    val params.finalName
    path(vec_fasta)
    path(vec_vb_fastq)
    path(cons)
    path(vec_gb)

    output:
    path "*_polyAdist.png", emit: polyA_dist_plots
    path "*_readLengthDist.png", emit: read_length_plots
    path "*_ist.csv", emit: dist_csv
    path "*allMutations.csv", emit: all_mut_csv
    path "*allPolyA_Info.csv", emit: all_polyA_csv

    shell:
    '''
    python getMutInfoFromVecOrVb.py -s !{params.excel} -o !{params.finalName} -ref !{vec.fasta} -fastq !{vec_vb_fastq} -cons !{cons} -gb !{vec_gb}
    '''
}