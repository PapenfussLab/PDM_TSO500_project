#!/usr/bin/env nextflow

projectDir="/vast/ai_projects/melanoma_depigmentation_biomarker_discovery/peinan/PDM_SCM"
params.metadata="$projectDir/data/bam_meta.csv"
params.refDir="$projectDir/bin/TSO500_RUO_LocalApp_v2/resources/genomes/hg19_hardPAR"
params.gnomAD="$projectDir/data/ref/gnomAD_hg19"
params.mapbility="$projectDir/data/ref/wgEncodeCrgMapabilityAlign100mer.bigWig"
params.TSObed="$projectDir/data/ref/resources/TST500C_manifest.bed"
params.Rlib="$projectDir/renv/library/linux-rhel-9.3/R-4.4/x86_64-pc-linux-gnu"
params.gatkPATH="$projectDir/data/GATK_CNA_2025"
params.intervalGATK="$projectDir/data/ref/resources/preprocessed_hg19.interval_list"
params.gvcfDir="$projectDir/data/normals/gvcf"
params.ponvcfDir="$projectDir/data/normals/ponvcf"

process generateInterval{
    publishDir path: "$projectDir/PureCN/output/ref", mode: 'copy'
    executor 'slurm'
    cpus = 2
    memory = 16.GB
    time = 1.hour

input:
    path (Rlib)
    path (refDir)
    path (mapbility)
    path (TSObed)
    

output: 
    path ("TSO500_intervals.txt")

script:
	"""
	runGenerateInterval.sh $Rlib $refDir $mapbility $TSObed
    """
}

process checksex {     
    executor 'slurm'
    cpus = 1
    memory = 4.GB
    time = 1.hour
input:
    tuple val (sample), path(bam), val(condition), path(bai)
    path (regionFile) 

output:
	tuple val (sample), path(bam), val(condition), path(bai), path("${sample}_sex.txt")
       
script:
	"""
	runSexCheck.sh $bam $sample $regionFile
    """
}

 process mutect2 {     
    publishDir path: "$projectDir/PureCN/output/mutect2", mode: 'copy'
    executor 'slurm'
    cpus = 20
    memory = 200.GB
    time = 48.hour

input:
    tuple val (sample), path(bam), val(condition), path(bai), val(sex)
    path (refDir)
    path (intervalGATK)
    path (ponvcfDir)
    path (gnomAD)


output:
    tuple val (sample), val(sex), path("${sample}_mutect2")
       
script:
    """
    runMutect2.sh $sample $bam $sex $refDir $intervalGATK $ponvcfDir $gnomAD
    """
}

process normalDB {
    publishDir path: "$projectDir/PureCN/output/ref", mode: 'copy'
    executor 'slurm'
    cpus = 8
    memory = 80.GB
    time = 48.hour

input:
    tuple path(files), val(sex)
    path (Rlib)
    path (intervalPureCN)
    path (gvcfDir)

output:
    path ("normal_${sex}")

script:
    """
    runNormalDB.sh $sex $Rlib $intervalPureCN $gvcfDir
    """
}

process PureCN_GATK {
    publishDir path: "$projectDir/PureCN/output/results", mode: 'copy'
    executor 'slurm'
    cpus = 8
    memory = 80.GB
    time = 48.hour

input:
    tuple val(sample), val(sex), path(mutect2PATH)
    path (gatkPATH)
    path (R_LIBS)
    path (normaldb)

output: 
    path ("${sample}_GATK_PureCN")

script:
    """
    runPureCN_GATK.sh $sample $sex $mutect2PATH $gatkPATH $R_LIBS
    """
}

workflow  {
    // Pre-process TSO500 interval for PureCN analysis
    generateInterval(params.Rlib, params.refDir, params.mapbility, params.TSObed)

    // Load input channel from metadata file
    input_ch=Channel
    .fromPath(params.metadata)
	.splitCsv( header:true )
	.map { row -> 
		[row.sample, row.bam, row.condition]
        row.bai = row.bam+".bai"
        return [row.sample, row.bam, row.condition, row.bai]
	}
    
    // check sex based on X/Y chromosome coverage and validated with clinical annotation
    meta_ch=checksex(input_ch, params.TSObed)
    .map { tuple ->
        def (sample, bam, condition, bai, sexfile) = tuple
        def sex = sexfile.text.trim()
        return [sample, bam, condition, bai, sex]
    }
    .branch {
		normal: it[2] == "normal"
		tumour: it[2] == "tumour"
	}
    
    tumour_ch=meta_ch.tumour

    // Run mutect2 on tumour samples, keep all germline variants
    tumour_mutect=mutect2(tumour_ch, params.refDir, params.intervalGATK, params.ponvcfDir, params.gnomAD)


    normal_ch=meta_ch.normal
    .branch{
        female: it[4] == "F"
		male: it[4] == "M"
    }

    normal_female_bam=normal_ch.female
    .map { [ it[1], it[3] ] }
    .collect()
    .map { [it, "F"] }
    
    normal_male_bam=normal_ch.male
    .map { [ it[1], it[3] ] }
    .collect()
    .map { [it, "M"] }

    normal_merged_bam=normal_female_bam
    .merge(normal_male_bam)
    .map { tuple ->
        def (file_F, sex_F, file_M, sex_M) = tuple
        [[file_F, sex_F], [file_M, sex_M]]
    }
    .flatMap()

    normal_coverage=normalDB(normal_merged_bam, params.Rlib, generateInterval.out, params.gvcfDir)
    .collect()

    PureCN_GATK(tumour_mutect, params.gatkPATH, params.Rlib, normal_coverage)
}
