
# stageDefaults contains the default options which are applied to each stage (command).
# This section is required for every Rubra pipeline.
# These can be overridden by options defined for individual stages, below.
# Stage options which Rubra will recognise are: 
#  - distributed: a boolean determining whether the task should be submitted to a cluster
#      job scheduling system (True) or run on the system local to Rubra (False). 
#  - walltime: for a distributed PBS job, gives the walltime requested from the job
#      queue system; the maximum allowed runtime. For local jobs has no effect.
#  - memInGB: for a distributed PBS job, gives the memory in Gigabytes requested from the 
#      job queue system. For local jobs has no effect.
#  - queue: for a distributed PBS job, this is the name of the queue to submit the
#      job to. For local jobs has no effect. This is currently a mandatory field for
#      distributed jobs, but can be set to None.
#  - modules: the modules to be loaded before running the task. This is intended for  
#      systems with environment modules installed. Rubra will call module load on each 
#      required module before running the task. Note that defining modules for individual 
#      stages will override (not add to) any modules listed here. This currently only
#      works for distributed jobs.
stageDefaults = {
    'distributed': True,
    'queue': None,
    'walltime': "01:00:00",
    'memInGB': 8,
    'modules': [
        "bwa-intel/0.6.2",
        "samtools-intel/0.1.19",
        "picard/1.53",
        "python-gcc/2.7.5",
        "R-gcc/3.0.2",
        "gatk/1.6-7"
    ]
}

# stages should hold the details of each stage which can be called by runStageCheck.
# This section is required for every Rubra pipeline.
# Calling a stage in this way carries out checkpointing and, if desired, batch job
# submission. 
# Each stage must contain a 'command' definition. See stageDefaults above for other 
# allowable options.
stages = {
    "fastqc": {
        "command": "fastqc --quiet -o %outdir %seq",
        'modules': [ "fastqc/0.10.1" ]
    },
    'alignBWA': {
        'command': "bwa aln -t 8 %encodingflag %ref %seq > %out",
        'walltime': "3:00:00",
        'queue': 'smp',
        'memInGB': 23
    },
    'alignToSamSE': {
        'command': "bwa samse %ref %meta %align %seq > %out"
    },
    'alignToSamPE': {
        'command': "bwa sampe %ref %meta %align1 %align2 %seq1 %seq2 > %out"
    },
    'samToSortedBam': {
        'command': "./SortSam 6 VALIDATION_STRINGENCY=LENIENT INPUT=%seq OUTPUT=%out SORT_ORDER=coordinate",
        'walltime': "5:00:00",
    },
    'mergeBams': {
        'command': "./PicardMerge 6 %baminputs USE_THREADING=true VALIDATION_STRINGENCY=LENIENT AS=true OUTPUT=%out",
        'walltime': "5:00:00"
    },
    'indexBam': {
        'command': "samtools index %bam"
    },
    'flagstat': {
        'command': "samtools flagstat %bam > %out",
        'walltime': "00:10:00"
    },
    'igvcount': {
        'command': "igvtools count %bam %out hg19",
        'modules': [ "igvtools/1.5.15" ]
    },
    'indexVCF': {
        'command': "./vcftools_prepare.sh %vcf",
        'modules': [ "tabix/0.2.5" ]
    },
    'realignIntervals': {
        # Hard-coded to take 2 known indels files right now
        'command': "./GenomeAnalysisTK 1 -T RealignerTargetCreator -R %ref -I %bam --known %indels_goldstandard --known %indels_1000G -log %log -o %out",
        'memInGB': 23,
        'walltime': "5:00:00"
    },
    'realign': {
        'command': "./GenomeAnalysisTK 22 -T IndelRealigner -R %ref -I %bam -targetIntervals %intervals -log %log -o %out",
        'memInGB': 23,
        'walltime': "5:00:00"
    },
    'dedup': {
        'command': "./MarkDuplicates 6 INPUT=%bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT AS=true METRICS_FILE=%log OUTPUT=%out",
        'walltime': '5:00:00'
    },
    'baseQualRecalCount': {
        'command': "./GenomeAnalysisTK 12 -T CountCovariates -I %bam -R %ref --knownSites %dbsnp -nt 8 -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -log %log -recalFile %out",
        'queue': 'smp',
        'memInGB': 23,
        'walltime': "5:00:00"
    },
    'baseQualRecalTabulate': {
        'command': "./GenomeAnalysisTK 4 -T TableRecalibration -I %bam -R %ref -recalFile %csvfile -l INFO -log %log -o %out",
        'walltime': "5:00:00"
    },
    'callSNPs': {
        'command': "./GenomeAnalysisTK 12 -T UnifiedGenotyper -nt 8 -R %ref -I %bam --dbsnp %dbsnp -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1600 -l INFO -A AlleleBalance -A DepthOfCoverage -A FisherStrand -glm SNP -log %log -o %out",
        'queue': 'smp',
        'memInGB': 23,
        'walltime': "3:00:00"
    },
    'callIndels': {
        'command': "./GenomeAnalysisTK 12 -T UnifiedGenotyper -nt 8 -R %ref -I %bam --dbsnp %dbsnp -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1600 -l INFO -A AlleleBalance -A DepthOfCoverage -A FisherStrand -glm INDEL -log %log -o %out",
        'queue': 'smp',
        'memInGB': 23,
        'walltime': "3:00:00"
    },
    'filterSNPs': {
        # Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
        'command': "./GenomeAnalysisTK 4 -T VariantFiltration -R %ref --variant %vcf --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'GATK_MINIMAL_FILTER' -log %log -o %out",
    },
    'filterIndels': {
        # Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
        # If you have 10 or more samples GATK also recommends the filter InbreedingCoeff < -0.8
        'command': "./GenomeAnalysisTK 4 -T VariantFiltration -R %ref --variant %vcf --filterExpression 'QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0' --filterName 'GATK_MINIMAL_FILTER' -log %log -o %out",
    },
    'annotateEnsembl': {
        # This command as written assumes that VEP and its cache have been
        # downloaded in respective locations
        # ./variant_effect_predictor_2.5
        # ./variant_effect_predictor_2.5/vep_cache
        'command': "perl variant_effect_predictor_2.5/variant_effect_predictor.pl --cache --dir variant_effect_predictor_2.5/vep_cache -i %vcf --vcf -o %out -species human --canonical --gene --protein --sift=b --polyphen=b > %log",
        'modules': [ "perl/5.10.1", "ensembl/67" ]
    },
    'depthOfCoverage': {
        'command': "./GenomeAnalysisTK 4 -T DepthOfCoverage -R %ref -I %bam -omitBaseOutput -ct 1 -ct 10 -ct 20 -ct 30 -o %out",
    },
    'collateReadcounts': {
        'command': 'python count_flagstat_wgs.py %dir %outdir',
        'walltime': "00:10:00"
    }
}
