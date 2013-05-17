working_files = {
    'output_dir': '/vlsci/VR0244/shared/test_ruffus_exome/output',
    'fastq_dirs': [
         '/vlsci/VR0244/shared/meg/EXOME/FASTQ'
    ],
    'fastq_symlink_dir': '/vlsci/VR0244/shared/test_ruffus_exome/fastq_symlinks'
}
ref_files = {
    'fasta_reference': '/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/bwa_Indexed/human_g1k_v37.fasta',
    'bwa_reference': '/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/bwa_Indexed/human_g1k_v37.fasta',
    'exon_bed': '/vlsci/VR0244/shared/Ruth_Exome/SureSelect_XT_Human_All_Exon_v4plus_UTR/TruSeqCustom_WOAgilentSSXTHAEv4TR_targeted_regions.bed',
    'exon_bed_extended': '/vlsci/VR0244/shared/Ruth_Exome/SureSelect_XT_Human_All_Exon_v4plus_UTR/TruSeqCustom_WOAgilentSSXTHAEv4TR_targeted_regions_plusMinus150bp.bed',
    'dbsnp': '/vlsci/VR0002/shared/Reference_Files/SNP_db/dbSNP137.vcf',  #note currently contains somatic variants too, no SAO filter
    'indels_realign_goldstandard': '/vlsci/VR0002/shared/Reference_Files/Indels_for_realignment/Mills_and_1000G_gold_standard.indels.b37.vcf',
    'indels_realign_1000G': '/vlsci/VR0002/shared/Reference_Files/Indels_for_realignment/1000G_phase1.indels.b37.vcf'
}
pipeline = {
    'logDir': 'log_exome_test',
    'logFile': 'pipeline.log',
    'style': 'print',
    'procs': 10,
    'verbose': 1,
    'end': ['earlyDepthOfCoverage', 'dedupedDepthOfCoverage', 'finalDepthOfCoverage',
            'fastqc', 
            'igvcountMergedBams', 'countRunBam', 
            'collateReadCounts',
            'vcfIndexSNPs', 'vcfIndexIndels',
            'getEnsemblAnnotations',
            'exonCoverage'],
    'force': [],
    'rebuild' : "fromstart",
    'restrict_samples': True,
    'allowed_samples': ['mw4241']
}
