
# This section is used by the variant calling pipeline.py to specify input data and 
# working directories.
#
# Note that if you have downloaded the pipeline the directory names below are examples 
# only and you will need to edit them to suit your needs.
#
# Required variables:
#  - fastq_dirs: a list of directories where the raw input data is found. Currently this
#      data is expected to be paired-end gzipped fastq and to follow a specific naming
#      convention (see below).
#  - fastq_symlink_dir: symlinks to all raw fastq files will be written to this directory
#      and used by the rest of the pipeline. These symlinks have standardised names and
#      are a useful flattened summary of all known input data.
#  - output_dir: the directory used by the pipeline for output and intermediate files.
#      A directory structure will be created under this directory by pipeline.py.
#
# Input data naming convention:
# The input fastq files must follow a naming convention so that the pipeline can determine
# the metadata fields. This convention in the default script is to use the regex
#     ([a-zA-Z0-9-.]+)_([^_/]+)_[CAGTN]+_L([0-9]+)_R(1|2).fastq.gz
# This corresponds to metadata fields
#    SAMPLE_RUN_TAG_LANE_READPAIR.fastq.gz
# where
#    SAMPLE is a unique identifier for the sample sequenced
#    RUN is a unique identifier for the experiment (e.g. run or flowcell ID)
#    TAG is the barcode sequence used for multiplexing (NA if none)
#    LANE is the flowcell lane identifier, written like L001
#    READPAIR identifies whether the file contains forward or reverse reads, R1 or R2
#
# For example: Sample395_C0WK7ACXX_ACTTGA_L007_R1.fastq.gz
#
# This file naming convention follows that returned by many sequencing centres for
# Illumina data.    
#
working_files = {
    'fastq_dirs': [
         './example_data/input_data_wgs'
    ],
    'fastq_symlink_dir': './example_data/output_wgs/fastq_symlinks',
    'output_dir': './example_data/output_wgs'
}

# This section is used by the variant calling pipeline.py to specify reference data files.
#
# Note that if you have downloaded the pipeline the filenames below are examples only and
# you will need to get the relevant reference files for your data. Exactly which files
# you need depend on your data. At time of writing reference data can be obtained from:
#  - Reference genome: many sources depending on data. For our human data we used the 
#     1000 genomes version of the b37 (hg19) genome build, found at
#     ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
#     Note that your genome must use the same chromosome naming convention as any other
#     reference files (such as dbSNP); if you use hg19 (chr1,chr2) instead of b37 (1,2)
#     you may need to convert the files suggested below.
#  - dbSNP variants: dbSNP is at http://www.ncbi.nlm.nih.gov/projects/SNP/
#     A useful release summary is at http://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi
#     We used human variants which were obtained in VCF format from
#     ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
#  - known indels for local realignment: these follow the Broad recommendations for the
#     GATK tool suite and come from the GATK resource bundle. See
#     http://gatkforums.broadinstitute.org/discussion/1213/what-s-in-the-resource-bundle-and-how-can-i-get-it
#
# Expected variables (if you use the relevant pipeline steps):
#  - fasta_reference: the reference genome fasta. Should be in the same location as the
#      .fai files produced by samtools faidx. 
#      TODO: do this indexing as part of the pipeline and check for index files.
#  - bwa_reference: the reference genome fasta. Should be in the same location as the
#      index files produced by bwa index.
#      TODO: do this indexing as part of the pipeline and check for index files.
#  - dbsnp: the dbSNP variants file in VCF format, for annotating variants and for 
#       GATK base quality recalibration.
#  - indels_realign_goldstandard and
#  - indels_realign_1000G: files of known indels for use in GATK local realignment.
#       Currently the Broad Institute recommends using these two files (see above).
ref_files = {
    'fasta_reference': '/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/bwa_Indexed/human_g1k_v37.fasta',
    'bwa_reference': '/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/bwa_Indexed/human_g1k_v37.fasta',
    
    'dbsnp': '/vlsci/VR0002/shared/Reference_Files/SNP_db/dbSNP137.vcf',  

    'indels_realign_goldstandard': '/vlsci/VR0002/shared/Reference_Files/Indels_for_realignment/Mills_and_1000G_gold_standard.indels.b37.vcf',
    'indels_realign_1000G': '/vlsci/VR0002/shared/Reference_Files/Indels_for_realignment/1000G_phase1.indels.b37.vcf'
}

# pipeline should hold configuration options for Rubra and for the pipeline.
# This section is required for every Rubra pipeline,
# but restrict_samples and allowed_samples are specific to the variant-calling pipeline.
#
# Rubra variables:
#  - logDir: the directory where batch queue scripts, stdout and sterr dumps are stored.
#  - logFile: the file used to log all jobs that are run.
#  - style: the default style, one of 'flowchart', 'print', 'run', 'touchfiles'. Can be 
#      overridden by specifying --style on the command line.
#  - procs: the number of python processes to run simultaneously. This determines the
#      maximum parallelism of the pipeline. For distributed jobs it also constrains the
#      maximum total jobs submitted to the queue at any one time.
#  - verbosity: one of 0 (quiet), 1 (normal), 2 (chatty). Can be overridden by specifying
#      --verbose on the command line.
#  - end: the desired tasks to be run. Rubra will also run all tasks which are dependencies 
#      of these tasks. Can be overridden by specifying --end on the command line.
#  - force: tasks which will be forced to run, regardless of timestamps. Can be overridden
#      by supplying --force on the command line.
#  - rebuild: one of 'fromstart','fromend'. Whether to calculate which dependencies will
#      be rerun by working back from an end task to the latest up-to-date task, or forward
#      from the earliest out-of-date task. 'fromstart' is the most conservative and 
#      commonly used as it brings all intermediate tasks up to date.
#
# Variant-calling pipeline variables:  (TODO: move to a separate section)
#  - restrict_samples: whether to restrict input files to those specified by allowd_samples
#  - allowed_samples: sample names that will be run of restrict_samples is True
pipeline = {
    'logDir': 'log_example_wgs',    
    'logFile': 'pipeline.log',
    'style': 'print',
    'procs': 30,
    'verbose': 1,
    'end': ['earlyDepthOfCoverage', 'finalDepthOfCoverage', 
            'fastqc',
            'igvcountMergedBams', 'countDedupedBam', 'countRunBam', 'countMergedBam',
            'getEnsemblAnnotations',
            'collateReadCounts',
            'vcfIndexSNPs', 'vcfIndexIndels'
            ],
    'force': [],
    'rebuild' : "fromstart",

    'restrict_samples': False,
    'allowed_samples': []
}
