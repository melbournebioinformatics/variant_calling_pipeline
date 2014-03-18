#!/bin/env python

"""
GATK-based variant-calling pipeline, WGS version.

Authors: Bernie Pope, Clare Sloggett, Gayle Philip.
Thanks to Dmitri Mouradov and Maria Doyle for input on the initial 
analysis design.
Thanks to Matt Wakefield for contributions to Rubra 
(https://github.com/bjpop/rubra) during pipeline development.

Description:

This program implements a workflow pipeline for next generation
sequencing variant detection using the Broad Institute's GATK for
variant calling and using ENSEMBL for basic annotation.

It uses Rubra (https://github.com/bjpop/rubra) based on the 
Ruffus library.

It supports parallel evaluation of independent pipeline stages,
and can run stages on a cluster environment.

The pipeline is configured by an options file in a python file,
including the actual commands which are run at each stage.
"""


import sys
import re
import os.path
import os
from collections import defaultdict
from glob import *
import shutil
from ruffus import *
from rubra.utils import pipeline_options
from rubra.utils import (runStageCheck, mkLogFile, mkDir, mkForceLink)
from input_fastq import parse_and_link

def make_metadata_string(metadata):
    return r'-R"@RG\tID:%s\tSM:%s\tPL:%s"' % (metadata['ID'], metadata['SM'], metadata['PL'])

# Shorthand access to options
ref_files = pipeline_options.ref_files
working_files = pipeline_options.working_files
logDir = pipeline_options.pipeline['logDir']

# Data setup process and input organisation and metadata functions

#Metadata holding structures
fastq_metadata = defaultdict(dict)

original_fastq_files = []
for fastq_dir in working_files['fastq_dirs']:
    original_fastq_files += glob(os.path.join(fastq_dir, '*.fastq.gz'))

if len(original_fastq_files)==0:
    print "No input files found. Do the filenames follow the naming convention?"
    print "Directories searched:"
    print "\n".join(working_files['fastq_dirs'])
    sys.exit(1)
    
# Parse metadata out of input file names and construct symlinks
# Metadata is put into a dict (for the rest of ruffus) and some of it also into symlinks (for filename uniqueness)
# currently parsing by assuming AGRF naming structure and paired-end reads
mkDir(working_files['fastq_symlink_dir'])
all_fastq_files = []
for file in original_fastq_files:
    symlink = parse_and_link(file, working_files['fastq_symlink_dir'], fastq_metadata)
    all_fastq_files.append(symlink)

# Make a list of files we will actually use
if pipeline_options.pipeline['restrict_samples']:
    allowed_samples = set(pipeline_options.pipeline['allowed_samples'])
    fastq_files = [file for file in sorted(all_fastq_files) 
                        if (fastq_metadata[os.path.basename(file)]['sample'] in allowed_samples)]
else:
    fastq_files = sorted(all_fastq_files)

print "Symlinked files that will be used:"
for file in fastq_files:
    print file
print
print "Output dir is %s" % working_files['output_dir']
print "Log dir is %s" % logDir
print

# Create output subdirectories

output_dir = working_files['output_dir']

fastqc_dir = os.path.join(output_dir, "FastQC")
mkDir(fastqc_dir)

sambam_dir = os.path.join(output_dir, "alignments")
mkDir(sambam_dir)

variant_dir = os.path.join(output_dir, "variant_calls")
mkDir(variant_dir)

coverage_dir = os.path.join(output_dir, "coverage")
mkDir(coverage_dir)

ensembl_dir = os.path.join(output_dir, "ensembl")
mkDir(ensembl_dir)

# directory for final summary tables
results_dir = os.path.join(output_dir, "results")
mkDir(results_dir)

# Pipeline declarations

# Alignment and correction steps

@transform(fastq_files, regex('(.+\/)?(.+?)\.fastq\.gz'), 
        [r'%s/\2_fastqc' % fastqc_dir, r'%s/\2.fastqc.Success' % fastqc_dir])
def fastqc(inputs, outputs):
    """
    Run FastQC on each fastq file.
    """
    sequence = inputs
    fastqc_dest, flagFile = outputs
    runStageCheck('fastqc', flagFile, fastqc_dir, sequence)

@collate(fastq_files, regex(r".*?([^/]+)(_1|_2)\.fastq.gz"), 
        [r"%s/\1.sam" % sambam_dir, r"%s/\1.bwaPE.Success" % sambam_dir])
def bwaPE(inputs, outputs):
    """
    Aligns two paired-end fastq files to a reference genome to produce a sam file.
    """
    seq1, seq2 = sorted(inputs)
    output, flag_file = outputs
    fastq_name = os.path.basename(seq1)
    sample = fastq_metadata[fastq_name]['sample']
    runID = fastq_metadata[fastq_name]['run_id']
    lane = fastq_metadata[fastq_name]['lane']
    readgroup_metadata = { 'PL': 'ILLUMINA',
                           'SM': sample,
                           'ID': "%s_%s_Lane%d" % (sample, runID, lane) }
    metadata_str = make_metadata_string(readgroup_metadata)
    print "bwa-mem on %s and %s" % (os.path.basename(seq1), os.path.basename(seq2))
    runStageCheck('bwaMemPE', flag_file, metadata_str, ref_files['bwa_reference'], seq1, seq2, output)

@transform(bwaPE, suffix(".sam"),
            [".bam", ".samToBam.Success"])
def samToBam(inputs, outputs):
    """
    Convert sam to bam and sort, using Picard.
    """
    output, flag_file = outputs
    sam, _success = inputs
    print "converting to sorted bam: %s" % os.path.basename(sam)
    runStageCheck('samToSortedBam', flag_file, sam, output)

@collate(samToBam, regex(r'(.*?)([^/_]+)_([^/_]+_[^/_]+)\.bam'), 
            [r"\1\2.bam", r'\1\2.mergeBams.Success'])
def mergeBams(inputs, outputs):
    """
    Merge the sorted bams together for each sample.
    Picard should cope correctly if there is only one input.
    """
    bams = [bam for [bam, _success] in inputs]
    output, flag_file = outputs
    baminputs = ' '.join(["INPUT=%s" % bam for bam in bams])
    print "merging %s into %s" % (",".join([os.path.basename(bam) for bam in bams]), os.path.basename(output))
    runStageCheck('mergeBams', flag_file, baminputs, output)

@follows('indexMergedBams')
@transform(mergeBams, suffix('.bam'), 
            ['.dedup.bam', '.bam.dedup.Success'])
def dedup(inputs, outputs):
    """
    Remove apparent duplicates from merged bams using Picard MarkDuplicates.
    """
    input_bam, _success = inputs
    output_bam, flag_file = outputs
    logFile = mkLogFile(logDir, input_bam, '.dedup.log')
    print "de-duping %s" % os.path.basename(input_bam)
    runStageCheck('dedup', flag_file, input_bam, logFile, output_bam)

@follows('indexDedupedBams')  
@transform(dedup, suffix('.bam'), 
            ['.realigner.intervals', '.bam.realignIntervals.Success'])
def realignIntervals(inputs, outputs):
    """
    Run GATK RealignTargetCreator to find suspect intervals for realignment.
    """
    bam, _success = inputs
    output_intervals, flag_file = outputs
    logFile = mkLogFile(logDir, bam, '.realignIntervals.log')
    print "calculating realignment intervals for %s" % os.path.basename(bam)
    runStageCheck('realignIntervals', flag_file, ref_files['fasta_reference'], bam, ref_files['indels_realign_goldstandard'], ref_files['indels_realign_1000G'], logFile, output_intervals)

def remove_GATK_bai(bamfile):
    """
    A bug in some versions of GATK cause it to create an x.bai file, and this gets in the way of using the properly named x.bam.bai file. If the given file exists, delete it.
    """
    bad_bai = os.path.splitext(bamfile)[0] + ".bai"
    try:
        os.remove(bad_bai)
    except OSError, e:
        # Ignore error only if it is OSError #2, ie File Not Found
        if e.errno != 2:
            raise e

@transform(realignIntervals, regex(r"(.*?)([^/]+)\.realigner\.intervals"),
            add_inputs([r'\1\2.bam']), 
            [r'\1\2.realigned.bam', r'\1\2.bam.realign.Success'])
def realign(inputs, outputs):
    """
    Run GATK IndelRealigner for local realignment, using intervals found by realignIntervals.
    """
    [intervals, _success], [input_bam] = inputs
    output_bam, flag_file = outputs
    logFile = mkLogFile(logDir, input_bam, '.realign.log')
    print "realigning %s" % os.path.basename(input_bam)
    runStageCheck('realign', flag_file, ref_files['fasta_reference'], input_bam, intervals, logFile, output_bam)
    remove_GATK_bai(output_bam)

@follows('indexRealignedBams')
@transform(realign, suffix('.bam'),
            ['.recal_data.csv', '.baseQualRecalCount.Success'])
def baseQualRecalCount(inputs, outputs):
    """
    GATK CountCovariates, first step of base quality score recalibration.
    """
    bam, _success = inputs
    output_csv, flag_file = outputs
    logFile = mkLogFile(logDir, bam, '.baseQualRecalCount.log')
    print "count covariates using GATK for base quality score recalibration: %s" % os.path.basename(bam)
    runStageCheck('baseQualRecalCount', flag_file, bam, ref_files['fasta_reference'], ref_files['dbsnp'], logFile, output_csv)

@transform(baseQualRecalCount, regex(r'(.*?)([^/]+)\.recal_data\.csv'), 
            add_inputs([r'\1\2.bam']), 
            [r'\1\2.recal.bam', r'\1\2.baseQualRecalTabulate.Success'])
def baseQualRecalTabulate(inputs, outputs):
    """
    GATK TableRecalibration: recalibrate base quality scores using the output of CountCovariates.
    """
    [input_csv, _success], [input_bam] = inputs
    output_bam, flag_file = outputs
    logFile = mkLogFile(logDir, input_bam, '.baseQualRecalTabulate.log')
    print "recalibrate base quality scores using GATK on %s" % os.path.basename(input_bam)
    runStageCheck('baseQualRecalTabulate', flag_file, input_bam, ref_files['fasta_reference'], input_csv, logFile, output_bam)
    remove_GATK_bai(output_bam)

# Temporarily putting this indexing step here to work around bug
@transform(baseQualRecalTabulate, suffix('.bam'),
            ['.bam.bai', '.bam.indexRecalibratedBams.Success'])
def indexRecalibratedBams(inputs, outputs):
    """
    Index the recalibrated bams using samtools. 
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "samtools index on %s" % os.path.basename(bam)
    runStageCheck('indexBam', flag_file, bam)

# Variant calling steps

@follows(indexRecalibratedBams)
@transform(baseQualRecalTabulate, 
            regex(r'(.*?)([^/]+)\.recal\.bam'),
            [r'%s/\2.SNP.vcf' % variant_dir, 
             r'%s/\2.SNP.vcf.idx' % variant_dir, 
             r'%s/\2.callSNPs.Success' % variant_dir])
def callSNPs(inputs, outputs):
    """
    Use GATK UnifiedGenotyper to call SNPs from recalibrated bams.
    """
    bam, _success = inputs
    output_vcf, _idx, flag_file = outputs
    logFile = mkLogFile(logDir, bam, '.callSNPs.log')
    print "calling SNPs from %s" % bam
    runStageCheck('callSNPs', flag_file, ref_files['fasta_reference'], bam, ref_files['dbsnp'], logFile, output_vcf)

@follows(indexRecalibratedBams)
@transform(baseQualRecalTabulate, 
            regex(r'(.*?)([^/]+)\.recal\.bam'),
            [r'%s/\2.INDEL.vcf' % variant_dir, 
             r'%s/\2.INDEL.vcf.idx' % variant_dir, 
             r'%s/\2.callIndels.Success' % variant_dir])
def callIndels(inputs, outputs):
    """
    Use GATK UnifiedGenotyper to call indels from recalibrated bams.
    """
    bam, _success = inputs
    output_vcf, _idx, flag_file = outputs
    logFile = mkLogFile(logDir, bam, '.callIndels.log')
    print "calling Indels from %s" % bam
    runStageCheck('callIndels', flag_file, ref_files['fasta_reference'], bam, ref_files['dbsnp'], logFile, output_vcf)

@transform(callSNPs, suffix('.SNP.vcf'),
            ['.SNP.filtered.vcf', '.SNP.filtered.vcf.idx', '.filterSNPs.Success'])
def filterSNPs(inputs, outputs):
    """
    Use GATK VariantFiltration to filter raw SNP calls.
    """
    input_vcf, _idx, _success = inputs
    output_vcf, _idxout, flag_file = outputs
    logFile = mkLogFile(logDir, input_vcf, '.filterSNPs.log')
    print "filtering SNPs from %s" % input_vcf
    runStageCheck('filterSNPs', flag_file, ref_files['fasta_reference'], input_vcf, logFile, output_vcf)

@transform(callIndels, suffix('.INDEL.vcf'),
            ['.INDEL.filtered.vcf', '.INDEL.filtered.vcf.idx', '.filterIndels.Success'])
def filterIndels(inputs, outputs):
    """
    Use GATK VariantFiltration to filter raw INDEL calls.
    """
    input_vcf, _idx, _success = inputs
    output_vcf, _idxout, flag_file = outputs
    logFile = mkLogFile(logDir, input_vcf, '.filterIndels.log')
    print "filtering indels from %s" % input_vcf
    runStageCheck('filterIndels', flag_file, ref_files['fasta_reference'], input_vcf, logFile, output_vcf)


@transform([filterSNPs, filterIndels], regex(r'.*?([^/]+)\.vcf'), 
    [r'%s/\1.ensembl.vcf' % ensembl_dir,r'%s/\1.getEnsemblAnnotations.Success' % ensembl_dir])
def getEnsemblAnnotations(inputs, outputs):
    """
    Annotate vcf using ENSEMBL variant effect predictor.
    """
    vcf, _idx, _success = inputs
    output, flag_file = outputs
    logFile = mkLogFile(logDir, vcf, '.EnsemblAnnotation.log')
    print "Annotating %s with ENSEMBL variant effect predictor" % os.path.basename(vcf)
    runStageCheck('annotateEnsembl', flag_file, vcf, output, logFile)


# Indexing steps

@transform(mergeBams, suffix('.bam'),
            ['.bam.bai', '.bam.indexMergedBams.Success'])
def indexMergedBams(inputs, outputs):
    """
    Index the merged bams using samtools.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "samtools index on %s" % os.path.basename(bam)
    runStageCheck('indexBam', flag_file, bam)

@transform(dedup, suffix('.bam'),
            ['.bam.bai', '.bam.indexDedupedBams.Success'])
def indexDedupedBams(inputs, outputs):
    """
    Index the de-duplicated bams using samtools. Note that this actually goes from the fixMate-ed bams.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "samtools index on %s" % os.path.basename(bam)
    runStageCheck('indexBam', flag_file, bam)

@transform(realign, suffix('.bam'),
            ['.bam.bai', '.bam.indexRealignedBams.Success'])
def indexRealignedBams(inputs, outputs):
    """
    Index the locally realigned bams using samtools.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "samtools index on %s" % os.path.basename(bam)
    runStageCheck('indexBam', flag_file, bam)



@transform(mergeBams, suffix('.bam'),
            ['.bam.tdf', '.bam.igvcountMergedBams.Success'])
def igvcountMergedBams(inputs, outputs):
    """
    Use igvtools count to create a .tdf file for the merged bam files, to improve viewing of the bam coverage in igv.
    """
    bam, _success = inputs
    outfile, flag_file = outputs
    print "igvtools count on %s" % os.path.basename(bam)
    runStageCheck('igvcount', flag_file, bam, outfile)

@transform(realign, suffix('.bam'),
            ['.bam.tdf', '.bam.igvcountRealignedBams.Success'])
def igvcountRealignedBams(inputs, outputs):
    """
    Use igvtools count to create a .tdf file for the merged bam files, to improve viewing of the bam coverage in igv.
    """
    bam, _success = inputs
    outfile, flag_file = outputs
    print "igvtools count on %s" % os.path.basename(bam)
    runStageCheck('igvcount', flag_file, bam, outfile)

@transform(dedup, suffix('.bam'),
            ['.bam.tdf', '.bam.igvcountDedupedBams.Success'])
def igvcountDedupedBams(inputs, outputs):
    """
    Use igvtools count to create a .tdf file for the deduped bam files, to improve viewing of the bam coverage in igv. Note that this actually goes from the fixMate-ed bams.
    """
    bam, _success = inputs
    outfile, flag_file = outputs
    print "igvtools count on %s" % os.path.basename(bam)
    runStageCheck('igvcount', flag_file, bam, outfile)

@transform(baseQualRecalTabulate, suffix('.bam'),
            ['.bam.tdf', '.bam.igvcountRecalibratedBams.Success'])
def igvcountRecalibratedBams(inputs, outputs):
    """
    Use igvtools count to create a .tdf file for the recalibrated bam files, to improve viewing of the bam coverage in igv.
    """
    bam, _success = inputs
    outfile, flag_file = outputs
    print "igvtools count on %s" % os.path.basename(bam)
    runStageCheck('igvcount', flag_file, bam, outfile)

@transform(filterSNPs, suffix('.vcf'),
            ['.vcf.gz', '.vcf.gz.tbi', '.vcfindexSNPs.Success'])
def vcfIndexSNPs(inputs, outputs):
    """
    Use bgzip and tabix to prepare raw SNPs vcf for vcftools handling.
    """
    vcf, _idx, _success = inputs
    zipfile, tabix_index, flag_file = outputs
    print "bgzip and tabix (for vcftools) on %s" % vcf
    runStageCheck('indexVCF', flag_file, vcf)

@transform(filterIndels, suffix('.vcf'),
            ['.vcf.gz', '.vcf.gz.tbi', '.vcfindexIndels.Success'])
def vcfIndexIndels(inputs, outputs):
    """
    Use bgzip and tabix to prepare raw indels vcf for vcftools handling.
    """
    vcf, _idx, _success = inputs
    zipfile, tabix_index, flag_file = outputs
    print "bgzip and tabix (for vcftools) on %s" % vcf
    runStageCheck('indexVCF', flag_file, vcf)


# Coverage steps

@follows(indexMergedBams)
@transform(mergeBams, 
            regex(r'(.*?)([^/]+)\.bam'),
            [r'%s/\2.early.DepthOfCoverage.sample_cumulative_coverage_counts' % coverage_dir, 
            r'%s/\2.early.DepthOfCoverage.sample_cumulative_coverage_proportions' % coverage_dir, 
            r'%s/\2.early.DepthOfCoverage.sample_interval_statistics' % coverage_dir, 
            r'%s/\2.early.DepthOfCoverage.sample_interval_summary' % coverage_dir, 
            r'%s/\2.early.DepthOfCoverage.sample_statistics' % coverage_dir, 
            r'%s/\2.early.DepthOfCoverage.sample_summary' % coverage_dir, 
            r'%s/\2.earlyDepthOfCoverage.Success' % coverage_dir])
def earlyDepthOfCoverage(inputs, outputs):
    """
    Use GATK DepthOfCoverage to get a first pass at coverage statistics, after merging bams.
    """
    bam, _success = inputs
    flag_file = outputs[-1]
    output_example = outputs[0]
    output_base = os.path.splitext(output_example)[0]
    print "calculating coverage statistics using GATK DepthOfCoverage on %s" % bam
    runStageCheck('depthOfCoverage', flag_file, ref_files['fasta_reference'], bam, output_base)

@follows(indexDedupedBams)
@transform(dedup, 
        regex(r'(.*?)([^/]+)\.dedup\.bam'),
        [r'%s/\2.deduped.DepthOfCoverage.sample_cumulative_coverage_counts' % coverage_dir, 
         r'%s/\2.deduped.DepthOfCoverage.sample_cumulative_coverage_proportions' % coverage_dir, 
         r'%s/\2.deduped.DepthOfCoverage.sample_interval_statistics' % coverage_dir, 
         r'%s/\2.deduped.DepthOfCoverage.sample_interval_summary' % coverage_dir, 
         r'%s/\2.deduped.DepthOfCoverage.sample_statistics' % coverage_dir, 
         r'%s/\2.deduped.DepthOfCoverage.sample_summary' % coverage_dir, 
         r'%s/\2.dedupedDepthOfCoverage.Success' % coverage_dir])
def dedupedDepthOfCoverage(inputs, outputs):
    """
    Use GATK DepthOfCoverage to get a coverage statistics as soon as duplicates are removed.
    """
    bam, _success = inputs
    flag_file = outputs[-1]
    output_example = outputs[0]
    output_base = os.path.splitext(output_example)[0]
    print "calculating coverage statistics using GATK DepthOfCoverage on %s" % bam
    runStageCheck('depthOfCoverage', flag_file, ref_files['fasta_reference'], bam, output_base)

@follows(indexRecalibratedBams)
@transform(baseQualRecalTabulate, 
            regex(r'(.*?)([^/]+)\.recal\.bam'),
            [r'%s/\2.DepthOfCoverage.sample_cumulative_coverage_counts' % coverage_dir, 
            r'%s/\2.DepthOfCoverage.sample_cumulative_coverage_proportions' % coverage_dir, 
            r'%s/\2.DepthOfCoverage.sample_interval_statistics' % coverage_dir, 
            r'%s/\2.DepthOfCoverage.sample_interval_summary' % coverage_dir, 
            r'%s/\2.DepthOfCoverage.sample_statistics' % coverage_dir, 
            r'%s/\2.DepthOfCoverage.sample_summary' % coverage_dir, 
            r'%s/\2.depthOfCoverage.Success' % coverage_dir])
def finalDepthOfCoverage(inputs, outputs):
    """
    Use GATK DepthOfCoverage to get coverage statistics.
    """
    bam, _success = inputs
    flag_file = outputs[-1]
    output_example = outputs[0]
    output_base = os.path.splitext(output_example)[0]
    print "calculating coverage statistics using GATK DepthOfCoverage on %s" % bam
    runStageCheck('depthOfCoverage', flag_file, ref_files['fasta_reference'], bam, output_base)


# Read-counting steps

@transform(samToBam, suffix('.bam'), 
            ['.bam.flagstat', '.bam.countRunBam.Success'])
def countRunBam(inputs, outputs):
    """
    Run samtools flagstat on the initial per-lane, per-run bam file.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "Running samtools flagstat on %s" % bam
    runStageCheck('flagstat', flag_file, bam, output)

@transform(mergeBams, suffix('.bam'), 
            ['.bam.flagstat', '.bam.countRunBam.Success'])
def countMergedBam(inputs, outputs):
    """
    Run samtools flagstat on the merged bam file.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "Running samtools flagstat on %s" % bam
    runStageCheck('flagstat', flag_file, bam, output)

@transform(realign, suffix('.bam'), 
            ['.bam.flagstat', '.bam.countRealignedBam.Success'])
def countRealignedBam(inputs, outputs):
    """
    Run samtools flagstat on the realigned bam file.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "Running samtools flagstat on %s" % bam
    runStageCheck('flagstat', flag_file, bam, output)

@transform(dedup, suffix('.bam'), 
            ['.bam.flagstat', '.bam.countDedupedBam.Success'])
def countDedupedBam(inputs, outputs):
    """
    Run samtools flagstat on the deduped bam file.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "Running samtools flagstat on %s" % bam
    runStageCheck('flagstat', flag_file, bam, output)


# Data collation and plotting steps

@merge([countDedupedBam, countMergedBam],
        ["%s/readcounts.txt" % results_dir, "%s/readcount_fractions.txt" % results_dir, "%s/collateReadcounts.Success" % results_dir])
def collateReadCounts(inputs, outputs):
    """
    Collate read counts from samtools flagstat output into a table.
    """
    # Note expected input and output directories are effectively hard-coded
    in_dir =  sambam_dir
    out_dir = results_dir  
    flag_file = outputs[-1]
    print "Collating read counts"
    runStageCheck('collateReadcounts', flag_file, in_dir, out_dir)


