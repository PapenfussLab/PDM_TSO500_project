#!/usr/bin/env nextflow

projectDir="/vast/ai_projects/melanoma_depigmentation_biomarker_discovery/peinan/PDM_SCM"
params.metadata="$projectDir/NM/NM_meta.csv"
params.refDir="$projectDir/bin/TSO500_RUO_LocalApp_v2/resources/genomes/hg19_hardPAR"
params.bamDir="$projectDir/NM/bam"
params.TSObamDir="$projectDir/NM/TSO_subset/bam"
params.Exomebed="$projectDir/NM/resource/SeqCap_EZ_Exome_v2_UCSC.bed"
params.TSObed="$projectDir/data/ref/resources/TST500C_manifest.bed"
params.gnomadDir="$projectDir/data/ref/gnomAD_hg19"

process preprocessInterval {
	publishDir path: "$projectDir/NM/TSO_CN_analysis/result/ref", mode: 'copy'
    executor 'slurm'
    cpus = 4
    memory = 32.GB
    time = 1.hour

input:
	path refDir
    path regionFile

output:
	path ("preprocessed_hg19.interval_list")
       
script:
    
	"""
	preprocessIntervals.sh $refDir $regionFile
    """
}

process checksex {    
    publishDir path: "$projectDir/NM/TSO_CN_analysis/result/temp", mode: 'copy'
    executor 'slurm'
    cpus = 2
    memory = 16.GB
    time = 1.hour
input:
    tuple val (sample), val(type), path(bam), path(bai)
    path (regionFile)
    
output:
	tuple val (sample), path("${sample}_sex.txt")
       
script:
	"""
	checkSex.sh $bam $sample $regionFile
    """
}

process normal_read_count {    
    publishDir path: "$projectDir/NM/TSO_CN_analysis/result/normal/readCount", mode: 'copy'
    executor 'slurm'
    cpus = 2
    memory = 32.GB
    time = 1.hour
input:
    tuple val (sample), val(type), path(bam), path(bai), val(sex)
    path (interval)
    
output:
	tuple val (sample), val(sex), path("${sample}_normal_counts.tsv")
       
script:
	"""
	runCollectReadCounts.sh $sample $bam $interval
    """
}

process read_count_PON {    
    publishDir path: "$projectDir/NM/TSO_CN_analysis/result/normal/readCountPON", mode: 'copy'
    executor 'local'
    cpus = 2
    memory = 16.GB
    time = 1.hour
input:
    tuple val (sex), file(readcounts)
    
output:
	path("read_count_PON_${sex}.hdf5")
       
script:
    
    readcount_files = readcounts.join(",")

    """
	createReadCountPON.sh $sex $readcount_files
    """
}

process GATK_CNV {    
    publishDir path: "$projectDir/NM/TSO_CN_analysis/result/GATK_CNV", mode: 'copy'
    executor 'slurm'
    cpus = 4
    memory = 32.GB
    time = 8.hour
input:
    tuple val (sample), val(type), path(bam), path(bai), val(sex)
    path (readCountPON)
    path (interval)
    path (refDir)
    
output:
	tuple val (sample), path ("${sample}_GATK")
       
script:
	"""
	runGATK_CNV.sh $sample $bam $sex $interval $refDir
    """
}

process GATK_CNV_plot {    
    publishDir path: "$projectDir/NM/TSO_CN_analysis/result/GATK_plot", mode: 'copy'
    executor 'slurm'
    cpus = 1
    memory = 8.GB
    time = 1.hour
input:
    tuple val (sample), path (GATKDir)
    path (refDir)
    
output:
	path ("${sample}.modeled.png")
       
script:
	"""
	plotGATK.sh $sample $GATKDir $refDir
    """
}

process mutect2_germline {     
    publishDir path: "$projectDir/NM/TSO_CN_analysis/result/normal/mutect2", mode: 'copy'
    executor 'slurm'
    cpus = 8
    memory = 80.GB
    time = 48.hour

input:
    tuple val (sample), val(type), path(bam), path(bai), val(sex)
    path (interval)
    path (refDir)
    path (gnomadDir)

output:
    tuple val (sample), val(sex), path("${sample}_normal_mutect2.vcf*") 
       
script:
    """
    mutect2_germline.sh $sample $bam $interval $refDir $gnomadDir
    """
}

process genomicsDB {     
    publishDir path: "$projectDir/NM/TSO_CN_analysis/result/normal/genomicsDB", mode: 'copy'
    executor 'slurm'
    cpus = 4
    memory = 80.GB
    time = 48.hour

input:
    tuple val (sex), path (mutect2Files)
    path (interval)
    path (refDir)
    

output:
    tuple val(sex), path("pondb_${sex}") 
       
script:
    """
    runGenomicsDBImport.sh $sex $interval $refDir
    """
}

process createPONVCF {     
    publishDir path: "$projectDir/NM/TSO_CN_analysis/result/normal/ponVCF", mode: 'copy'
    executor 'slurm'
    cpus = 4
    memory = 80.GB
    time = 48.hour

input:
    tuple val (sex), path (genomicsDB)
    path (refDir)
    

output:
    tuple val(sex), path("${sex}_pon.vcf.gz*") 
       
script:
    """
    createPonVCF.sh $sex $genomicsDB $refDir
    """
}

process createGVCF {     
    publishDir path: "$projectDir/NM/TSO_CN_analysis/result/normal/gVCF", mode: 'copy'
    executor 'slurm'
    cpus = 4
    memory = 80.GB
    time = 48.hour

input:
    tuple val (sex), path (genomicsDB)
    path (refDir)
    

output:
    tuple val(sex), path("${sex}.gvcf.gz*") 
       
script:
    """
    createGVCF.sh $sex $genomicsDB $refDir
    """
}



workflow  {
meta_ch=Channel.fromPath(params.metadata)
    .splitCsv( header:true )
    .map{ row -> [row.sample, row.type, row.capture] }
    .distinct()
    .branch {
		nimblegen: it[2] == "nimblegen"
		agilent: it[2] == "agilent" 
	}

WES_bam_ch=meta_ch.nimblegen
    .map{ sample, type, capture -> 
        def bam = params.bamDir + "/" + sample + "_" + type + ".bam" 
        def bai = params.bamDir + "/" + sample + "_" + type + ".bam.bai" 
        [ sample, type, file(bam), file(bai) ]
        }
    .branch {
		normal: it[1] == "normal"
		tumour: it[1] == "tumour"
	}

TSO_bam_ch=meta_ch.nimblegen
    .map{ sample, type, capture -> 
        def bam = params.TSObamDir + "/" + sample + "_" + type + "_TSO500.bam" 
        def bai = params.TSObamDir + "/" + sample + "_" + type + "_TSO500.bam.bai" 
        [ sample, type, file(bam), file(bai) ]
        }
    .branch {
		normal: it[1] == "normal"
		tumour: it[1] == "tumour"
	}

sample_sex=checksex(WES_bam_ch.normal, params.Exomebed)
    .map { tuple ->
        def (sample, sexfile) = tuple
        def sex = sexfile.text.trim()
        return [sample, sex]
    }

TSO_normal_ch=TSO_bam_ch.normal
    .join(sample_sex)

TSO_tumour_ch=TSO_bam_ch.tumour
    .join(sample_sex)

preprocessInterval(params.refDir, params.TSObed)
normal_readCount=normal_read_count(TSO_normal_ch, preprocessInterval.out)
    .branch{
        female: it[1] == "F"
		male: it[1] == "M"
    }
normal_readCount_M=normal_readCount.male
    .map{ it[2] }
    .collect()
    .map{ ["M", it ] }
normal_readCount_F=normal_readCount.female
    .map{ it[2] }
    .collect()
    .map{ ["F", it ] }
normal_readCount_combine=normal_readCount_M
    .mix(normal_readCount_F)

readCountPON=read_count_PON(normal_readCount_combine)
    .collect()

GATK_CNV(TSO_tumour_ch, readCountPON, preprocessInterval.out, params.refDir)
GATK_CNV_plot(GATK_CNV.out, params.refDir)

normal_VCF=mutect2_germline(TSO_normal_ch, preprocessInterval.out, params.refDir, params.gnomadDir)
    .branch{
        female: it[1] == "F"
		male: it[1] == "M"
    }
normal_VCF_F=normal_VCF.female
    .map { it[2] }
    .collect()
    .map{ ["F", it] }
normal_VCF_M=normal_VCF.male
    .map { it[2] }
    .collect()
    .map{ ["M", it] }
normal_VCF_combine=normal_VCF_F
    .mix(normal_VCF_M)

genomicsDB(normal_VCF_combine, preprocessInterval.out, params.refDir)
createPONVCF(genomicsDB.out, params.refDir)
createGVCF(genomicsDB.out, params.refDir)
}
