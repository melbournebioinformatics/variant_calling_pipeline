working_files = {
    'output_dir': '/vlsci/VR0244/shared/test_ruffus_wgs/output',
    'fastq_dirs': [
         '/vlsci/VR0244/shared/meg/EXOME/FASTQ'
    ],
    'fastq_symlink_dir': '/vlsci/VR0244/shared/test_ruffus_wgs/fastq_symlinks'
}
ref_files = {
    'fasta_reference': '/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/bwa_Indexed/human_g1k_v37.fasta',
    'bwa_reference': '/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/bwa_Indexed/human_g1k_v37.fasta',
    'dbsnp': '/vlsci/VR0002/shared/Reference_Files/SNP_db/dbSNP137.vcf',  #note currently contains somatic variants too, no SAO filter
    'indels_realign_goldstandard': '/vlsci/VR0002/shared/Reference_Files/Indels_for_realignment/Mills_and_1000G_gold_standard.indels.b37.vcf',
    'indels_realign_1000G': '/vlsci/VR0002/shared/Reference_Files/Indels_for_realignment/1000G_phase1.indels.b37.vcf'
}
pipeline = {
    'logDir': 'log_varcall_test',
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
    'restrict_samples': True,
    'allowed_samples': ['mw4241']
}
