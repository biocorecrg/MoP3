#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* 
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '2.0'

params.help            = false
params.resume          = false

log.info """

╔╦╗╔═╗╔═╗  ╔═╗┬─┐┌─┐┌─┐┬─┐┌─┐┌─┐┌─┐┌─┐┌─┐
║║║║ ║╠═╝  ╠═╝├┬┘├┤ ├─┘├┬┘│ ││  ├┤ └─┐└─┐
╩ ╩╚═╝╩    ╩  ┴└─└─┘┴  ┴└─└─┘└─┘└─┘└─┘└─┘
                                                                                       
====================================================
BIOCORE@CRG Master of Pores 2. Preprocessing - N F  ~  version ${version}
====================================================

conffile                  : ${params.conffile}

fast5                     : ${params.fast5}
fastq                     : ${params.fastq}

reference                 : ${params.reference}
annotation                : ${params.annotation}

granularity               : ${params.granularity}

ref_type                  : ${params.ref_type}
pars_tools                : ${params.pars_tools}
barcodes                  : ${params.barcodes}

output                    : ${params.output}

GPU                       : ${params.GPU}

basecalling               : ${params.basecalling} 
demultiplexing            : ${params.demultiplexing} 
demulti_fast5             : ${params.demulti_fast5}

filtering                 : ${params.filtering}
mapping                   : ${params.mapping}

counting                  : ${params.counting}
discovery                 : ${params.discovery}

cram_conv                 : ${params.cram_conv}
subsampling_cram          : ${params.subsampling_cram}


saveSpace                 : ${params.saveSpace}

email                     : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// check multi5 and GPU usage. GPU maybe can be removed as param if there is a way to detect it
if (params.GPU != "cuda11" && params.GPU != "cuda10" && params.GPU != "OFF" && params.GPU != "ON") exit 1, "Please specify cuda11, cuda10, ON or OFF if GPU processors are available. ON is legacy for cuda10"

// include functions, outdirs from other files
evaluate(new File("../outdirs.nf"))
def local_modules = file("$baseDir/../local_modules.nf")
def subworkflowsDir = "${baseDir}/../BioNextflow/subworkflows"
joinScript = file("$baseDir/bin/join.r")

// check input files
if (params.mapping != "NO") {
    reference = file(params.reference)
    if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
} else {
    reference = ""
}
config_report = file("$baseDir/config.yaml")
if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("$baseDir/../img/logo_small.png")
Channel.fromPath( "$baseDir/deeplexicon/*.h5").set{deepmodels}
Channel.fromPath( "$baseDir/seqtagger_models").set{seqtagger_models}



Channel
    .from( config_report, logo )
    .collect().set{multiqc_info}


def gpu = (params.GPU != 'OFF' ? 'ON' : 'OFF')
def cuda_cont = (params.GPU == 'cuda11' ? 'biocorecrg/mopbasecallc11:0.3' : 'biocorecrg/mopbasecall:0.3')

def tools = [:]
tools["basecalling"] = params.basecalling
tools["demultiplexing"] = params.demultiplexing
tools["mapping"] = params.mapping
tools["filtering"] = params.filtering
tools["counting"] = params.counting
tools["discovery"] = params.discovery
//tools["variantcall"] = params.variantcall

// Output files
outputReport   = file("${outputMultiQC}/multiqc_report.html")

/*
* move old multiQCreport
*/
if( outputReport.exists() ) {
  log.info "Moving old report to multiqc_report.html multiqc_report.html.old"
  outputReport.moveTo("${outputMultiQC}/multiqc_report.html.old")
}

/*
* This is default value in case guppy will be used for RNA demultiplexing
*/
params.barcodekit = ""

if (params.ref_type == "genome") {
    if (params.annotation != "") {
        annotation = file(params.annotation)
        if( !annotation.exists() ) exit 1, "Missing annotation file: ${params.annotation}!"
    }
}
def demulti_fast5_opt = "OFF"
if (params.demulti_fast5 == "ON" || params.demulti_fast5 == "YES" ) {
    demulti_fast5_opt = "ON"
}

def guppy_basecall_label = (params.GPU != 'OFF' ? 'basecall_gpus' : 'big_cpus')
def deeplexi_basecall_label = (params.GPU != 'OFF' ? 'demulti_gpus' : '')
def output_bc = (demulti_fast5_opt == 'ON' ? '' : outputFast5)
def outputMinionQC = (demulti_fast5_opt == 'ON' ? '': outputQual)


if (params.saveSpace == "YES") outmode = "move"
else outmode = "copy"

// Check the outmode with guppy and read-ducks!


include { RNA2DNA; preparing_demultiplexing_fast5_seqtagger; preparing_demultiplexing_fast5_deeplexicon; extract_seqtagger_fastq; extract_deeplexicon_fastq; parseFinalSummary; checkTools; reshapeSamples; reshapeDemuxSamples; checkRef; getParameters } from "${local_modules}" 
include { extract_demultiplexed_fast5 as extract_demultiplexed_fast5_seqtagger;  extract_demultiplexed_fast5 as extract_demultiplexed_fast5_deeplexicon } from "${local_modules}" addParams(OUTPUTF5: outputFast5, OUTPUTST: outputQual, LABEL: 'big_cpus')
include { extract_demultiplexed_fast5_guppy } from "${local_modules}" addParams(OUTPUT: outputFast5, LABEL: 'big_cpus')

def guppypars = parseFinalSummary(params.conffile)

// Create a channel for tool options
if (workflow.profile == "awsbatch") guppypars = guppypars + " --data_path /nextflow-bin/ont-guppy/data"
progPars = getParameters(params.pars_tools)
def guppy_basecall_pars = guppypars + " " + progPars["basecalling--guppy"]

// Create a channel for excluded ids
if (params.barcodes != "") {
    barcodes_to_include = file(params.barcodes)
    if( !barcodes_to_include.exists() ) exit 1, "Missing barcodes_to_include file: ${params.barcodes}!"
    Channel.from(barcodes_to_include.readLines())
    .map { line ->
        [ line ]
    }.set{ barcodes_to_include}
} else {
    barcodes_to_include = Channel.empty()
}


include { GET_WORKFLOWS; BASECALL as GUPPY_BASECALL; BASECALL_DEMULTI as GUPPY_BASECALL_DEMULTI } from "${subworkflowsDir}/basecalling/guppy" addParams(EXTRAPARS_BC: guppy_basecall_pars, EXTRAPARS_DEM: progPars["demultiplexing--guppy"], LABEL: guppy_basecall_label, GPU: gpu, MOP: "YES", OUTPUT: output_bc, CONTAINER: cuda_cont, OUTPUTMODE: outmode)
include { DEMULTIPLEX as READUCKS_DEMULTIPLEX } from "${subworkflowsDir}/demultiplexing/readucks" addParams(EXTRAPARS: progPars["demultiplexing--readucks"], LABEL: 'big_cpus', OUTPUT: output_bc, OUTPUTMODE: outmode)
include { DEMULTIPLEX as SEQTAGGER_DEMULTIPLEX } from "${subworkflowsDir}/demultiplexing/seq_tagger" addParams(EXTRAPARS: progPars["demultiplexing--seqtagger"], LABEL: 'demulti_gpus', OUTPUT: output_bc)

include { GET_VERSION as DEMULTIPLEX_VER; DEMULTIPLEX as DEMULTIPLEX_DEEPLEXICON } from "${subworkflowsDir}/demultiplexing/deeplexicon" addParams(EXTRAPARS: progPars["demultiplexing--deeplexicon"], LABEL:deeplexi_basecall_label, GPU: gpu)
include { GET_VERSION as NANOFILT_VER; FILTER as NANOFILT_FILTER} from "${subworkflowsDir}/trimming/nanofilt" addParams(EXTRAPARS: progPars["filtering--nanofilt"])
include { GET_VERSION as NANOQ_VER; FILTER as NANOQ_FILTER} from "${subworkflowsDir}/trimming/nanoq" addParams(EXTRAPARS: progPars["filtering--nanoq"])
include { MAP as GRAPHMAP} from "${subworkflowsDir}/alignment/graphmap" addParams(EXTRAPARS: progPars["mapping--graphmap"], LABEL:'big_mem_cpus')
include { MAP as GRAPHMAP2} from "${subworkflowsDir}/alignment/graphmap2" addParams(EXTRAPARS: progPars["mapping--graphmap2"], LABEL:'big_mem_cpus')
include { MAP as MINIMAP2} from "${subworkflowsDir}/alignment/minimap2" addParams(EXTRAPARS: progPars["mapping--minimap2"], LABEL:'big_mem_cpus')
include { ALL as BWA} from "${subworkflowsDir}/alignment/bwa" addParams(EXTRAPARS: progPars["mapping--bwa"], LABEL:'big_mem_cpus')
include { GET_VERSION as BWA_VER} from "${subworkflowsDir}/alignment/bwa" 
include { GET_VERSION as GRAPHMAP_VER} from "${subworkflowsDir}/alignment/graphmap" 
include { GET_VERSION as GRAPHMAP2_VER} from "${subworkflowsDir}/alignment/graphmap2" 
include { GET_VERSION as MINIMAP2_VER} from "${subworkflowsDir}/alignment/minimap2" 
include { FASTQCP as FASTQC} from "${subworkflowsDir}/qc/fastqc" addParams(LABEL: 'big_cpus')
include { GET_VERSION as FASTQC_VER} from "${subworkflowsDir}/qc/fastqc"
include { SORT as SAMTOOLS_SORT } from "${subworkflowsDir}/misc/samtools" addParams(LABEL: 'big_cpus', OUTPUT:outputMapping)
include { INDEX as SAMTOOLS_INDEX } from "${subworkflowsDir}/misc/samtools" addParams(OUTPUT:outputMapping)
include { GET_VERSION as SAMTOOLS_VERSION; CAT as SAMTOOLS_CAT } from "${subworkflowsDir}/misc/samtools"
include { MOP_QC as NANOPLOT_QC } from "${subworkflowsDir}/qc/nanoplot" addParams(LABEL: 'big_cpus_ignore')
include { GET_VERSION as NANOPLOT_VER } from "${subworkflowsDir}/qc/nanoplot" 
include { GET_VERSION as NANOCOUNT_VER } from "${subworkflowsDir}/read_count/nanocount"
include { COUNT as NANOCOUNT } from "${subworkflowsDir}/read_count/nanocount" addParams(LABEL: 'big_mem', EXTRAPARS: progPars["counting--nanocount"], OUTPUT:outputCounts)
include { COUNT_AND_ANNO as HTSEQ_COUNT } from "${subworkflowsDir}/read_count/htseq" addParams(CONTAINER:"biocorecrg/htseq:30e9e9c", EXTRAPARS: progPars["counting--htseq"], OUTPUT:outputCounts, LABEL:'big_cpus')
include { GET_VERSION as HTSEQ_VER } from "${subworkflowsDir}/read_count/htseq" addParams(CONTAINER:"biocorecrg/htseq:30e9e9c")

include { GET_VERSION as BAMBU_VER } from "${subworkflowsDir}/assembly/bambu" 
include { ASSEMBLE as BAMBU_ASSEMBLE } from "${subworkflowsDir}/assembly/bambu" addParams(EXTRAPARS: progPars["discovery--bambu"], OUTPUT:outputAssembly, LABEL:'big_mem_cpus')

include { GET_VERSION as ISOQUANT_VER } from "${subworkflowsDir}/assembly/isoquant" 
include { ASSEMBLE as ISOQUANT_ASSEMBLE } from "${subworkflowsDir}/assembly/isoquant" addParams(EXTRAPARS: progPars["discovery--isoquant"], OUTPUT:outputAssembly, LABEL:'big_mem_cpus', CONTAINER:'quay.io/biocontainers/isoquant:3.2.0--hdfd78af_0')

include { REPORT as MULTIQC; GET_VERSION as MULTIQC_VER } from "${subworkflowsDir}/reporting/multiqc" addParams(EXTRAPARS: "-c ${config_report.getName()}", OUTPUT:outputMultiQC)
include { concatenateFastQFiles} from "${local_modules}" addParams(OUTPUT:outputFastq)
include { MinIONQC} from "${local_modules}" addParams(OUTPUT:outputMinionQC, LABEL: 'big_mem_cpus')
include { bam2stats; countStats; joinCountStats; joinAlnStats} from "${local_modules}" 
include { cleanFile as fastqCleanFile; cleanFile as bamCleanFile; cleanFile as fast5CleanFile} from "${local_modules}"
include { AssignReads} from "${local_modules}" addParams(OUTPUT:outputAssigned)
include { bam2Cram } from "${local_modules}" addParams(OUTPUT:outputCRAM, LABEL: 'big_cpus_ignore')

/*
* Simple flow of basecalling
*/
workflow flow1 {
    take: 
        fast5_4_analysis
    main:
        outbc = GUPPY_BASECALL (fast5_4_analysis)
        basecalled_fastq = outbc.basecalled_fastq

        // Optional fastq filtering 
        if (params.filtering == "nanofilt") {
            basecalled_fastq = NANOFILT_FILTER(outbc.basecalled_fastq)
            //basecalled_fastq = reshapeSamples(nanofilt.out)
        } else if (params.filtering == "nanoq") {
            basecalled_fastq = NANOQ_FILTER(outbc.basecalled_fastq)
            //basecalled_fastq = reshapeSamples(nanofilt.out)
        }
        
        //bc_fastq = reshapeSamples(basecalled_fastq)
        bc_fast5 = reshapeSamples(outbc.basecalled_fast5)
        bc_stats = reshapeSamples(outbc.basecalling_stats)

    emit:
        basecalled_fast5 = bc_fast5
        basecalled_fastq = basecalled_fastq
        basecalled_stats = bc_stats

}

/*
*  Basecalling and Demultiplexing
*/

workflow flow2 {        
    take: 
        fast5_4_analysis
        barcodes_to_include
    main:
        // IF DEMULTIPLEXING IS DEEPLEXICON 
        if(params.demultiplexing == "deeplexicon") {
            outbc = GUPPY_BASECALL(fast5_4_analysis)

            demux = DEMULTIPLEX_DEEPLEXICON(deepmodels, fast5_4_analysis)
            fast5_res = outbc.basecalled_fast5
        
            // Optional demultiplex fast5       
            if (demulti_fast5_opt == "ON") {
                basecalledbc = reshapeSamples(outbc.basecalled_fast5)
                alldemux = reshapeSamples(demux)                
                
                prep_demux = preparing_demultiplexing_fast5_deeplexicon(alldemux.groupTuple()).transpose()
                data_for_demux = prep_demux.combine(basecalledbc.transpose().groupTuple(),  by: 0)              
                extract_demultiplexed_fast5_deeplexicon(data_for_demux)
                
                // OPTIONAL CLEANING FASTQ5 FILES
                if (params.saveSpace == "YES") {
                    fast5CleanFile(basecalledbc.transpose().groupTuple(), fast5_res.map{it[1]}.collect(), ".fast5")
                }
            }
        // Demultiplex fastq    
            demufq = extract_deeplexicon_fastq(demux.join(outbc.basecalled_fastq))

        } else if (params.demultiplexing == "guppy" || params.demultiplexing == "readucks") {
            // IF DEMULTIPLEXING IS GUPPY   
            outbc = GUPPY_BASECALL_DEMULTI (fast5_4_analysis)
            demufq = outbc.basecalled_fastq
            fast5_res = outbc.basecalled_fast5

            // In case of readucks demultiplexing       
            if (params.demultiplexing == "readucks" ) {
                read_demufq = READUCKS_DEMULTIPLEX(demufq)
                demufq = read_demufq
            }
            

            // Optional demultiplex fast5       
            if (demulti_fast5_opt == "ON" ) {
                basecalledbc = reshapeSamples(outbc.basecalled_fast5)
                alldemux = reshapeSamples(outbc.basecalling_stats)                              
                if (params.demultiplexing == "readucks" || params.demultiplexing == "seqtagger") {
                    println "################################################################"
                    println "WARNING FAST5 DEMULTIPLEXING NOT IMPLEMENTED WHEN USING READUCKS"
                    println "################################################################"
                } else {
                    fast5_res = extract_demultiplexed_fast5_guppy(alldemux.groupTuple().join(basecalledbc.transpose().groupTuple()))
                }
                // OPTIONAL CLEANING FASTQ5 FILES
                fast5CleanFile(basecalledbc.transpose().groupTuple(), fast5_res.map{it[1]}.collect(), ".fast5")
            }
        } else if (params.demultiplexing == "seqtagger" ) {
            outbc = GUPPY_BASECALL(fast5_4_analysis)
            fast5_res = outbc.basecalled_fast5

            demux = SEQTAGGER_DEMULTIPLEX(fast5_4_analysis, seqtagger_models)
            demufq = extract_seqtagger_fastq(demux.join(outbc.basecalled_fastq))

            // Optional demultiplex fast5       
            if (demulti_fast5_opt == "ON") {

                basecalledbc = reshapeSamples(outbc.basecalled_fast5)
                alldemux = reshapeSamples(demux)                
                prep_demux = preparing_demultiplexing_fast5_seqtagger(alldemux.groupTuple()).transpose()
                data_for_demux = prep_demux.combine(basecalledbc.transpose().groupTuple(),  by: 0)
                extract_demultiplexed_fast5_seqtagger(data_for_demux)
                
            }

        }
        
        reshapedPrefiltDemufq = demufq.transpose().map{
            [it[1].name.replace(".fastq.gz", "").replace(".fq.gz", ""), it[1] ]
        }
        
        if (params.barcodes != "") {    
            reshapedPrefiltDemufq.map{
                def id_raw = it[0].split("---")
                def id_raw2 = id_raw[1].split("\\.")
                def ori_id = "${id_raw[0]}---${id_raw2[1]}"
                [ori_id, it]
            }.join(barcodes_to_include).map{
                it[1]
            }.set{reshapedDemufq}
        } else {
            reshapedDemufq = reshapedPrefiltDemufq
        }
                
        basecalled_stats = reshapeSamples(outbc.basecalling_stats)

        // Optional fastq filtering 
        if (params.filtering == "nanofilt") {
            nanofilt = NANOFILT_FILTER(reshapedDemufq)
            reshapedDemufq = nanofilt
        } else if (params.filtering == "nanoq") {
            //nanofilt = NANOQ_FILTER(outbc.basecalled_fastq)
            //basecalled_fastq = reshapeSamples(nanofilt.out)
            nanofilt = NANOQ_FILTER(reshapedDemufq)
//          reshapedDemufq = nanofilt
        }
    emit:
        basecalled_fast5 =  fast5_res
        //basecalled_fastq = basecalled_fastq_res
        basecalled_fastq = reshapedDemufq
        basecalled_stats = basecalled_stats
            
}


workflow preprocess_flow {
    take:
        bc_fast5
        bc_fastq
        basecalled_stats
        
    main:   
    // Perform MinIONQC on basecalling stats
    basecall_qc = MinIONQC(basecalled_stats.groupTuple())   
    multiqc_data = basecall_qc.QC_folder.map{it[1]}.mix(multiqc_info)

    // Perform mapping on fastq files
    if (params.mapping == "NO") {
        stats_aln = Channel.value() 
        sorted_alns = Channel.value()   
        nanoplot_qcs = Channel.value()  
    }
    else {
        switch(params.mapping) { 
            case "graphmap": 
            //GRAPHMAP cannot align RNA
            dna_bc_fastq = RNA2DNA(bc_fastq)
            aln_reads = GRAPHMAP(dna_bc_fastq, reference)
            break
            case "graphmap2": 
            aln_reads = GRAPHMAP2(bc_fastq, reference)
            break
            case "minimap2": 
            aln_reads = MINIMAP2(bc_fastq, reference)
            break
            case "bwa": 
            aln_reads = BWA(reference, bc_fastq)
            break
            default: 
            println "ERROR ################################################################"
            println "${params.mapping} is not a supported alignment"
            println "ERROR ################################################################"
            println "Exiting ..."
            System.exit(0)
            break

        }    

        // Concatenate bamfiles
        if (params.demultiplexing == "NO" ) reshaped_aln_reads = reshapeSamples(aln_reads)
        else reshaped_aln_reads = reshapeDemuxSamples(aln_reads)

        jaln_reads = SAMTOOLS_CAT(reshaped_aln_reads.groupTuple())

        // Perform SORTING and INDEXING on bam files
        sorted_alns = SAMTOOLS_SORT(jaln_reads)
        aln_indexes = SAMTOOLS_INDEX(sorted_alns)

        // Converting BAM to CRAM and 
        if (params.cram_conv == "YES") {
            good_ref = checkRef(reference)
            bam2Cram(good_ref, params.subsampling_cram, sorted_alns.join(aln_indexes))
        }
        // OPTIONAL CLEANING BAM FILES
        if (params.saveSpace == "YES") {
            bamCleanFile(reshaped_aln_reads.groupTuple(), jaln_reads.map{it[1]}.collect(), ".bam")
        }
        // Perform bam2stats on sorted bams
        aln_stats = bam2stats(sorted_alns)
        stats_aln = joinAlnStats(aln_stats.map{ it[1]}.collect())
    
        // Perform NanoPlot on sorted bams
        nanoplot_qcs = NANOPLOT_QC(sorted_alns)
        multiqc_data = multiqc_data.mix(stats_aln)
    }   

    // Concatenate fastq files
    if (params.demultiplexing == "NO" ) reshaped_bc_fastq = reshapeSamples(bc_fastq)
    else reshaped_bc_fastq = reshapeDemuxSamples(bc_fastq)

    fastq_files = concatenateFastQFiles(reshaped_bc_fastq.groupTuple())

    // Perform fastqc QC on fastq
    fastqc_files = FASTQC(fastq_files)
    multiqc_data = multiqc_data.mix(fastqc_files.map{it[1]})

    // OPTIONAL CLEANING FASTQC FILES
    if (params.saveSpace == "YES") {
        fastqCleanFile(reshaped_bc_fastq.groupTuple(), fastq_files.map{it[1]}.collect().mix(fastqc_files.map{it[1]}.collect(), jaln_reads.map{it[1]}.collect()).collect(), ".gz")
    }
    // OPTIONAL Perform COUNTING / ASSIGNMENT
    if (params.counting == "nanocount" && params.ref_type == "transcriptome") {
        read_counts = NANOCOUNT(sorted_alns.join(aln_indexes))
        assignments = AssignReads(sorted_alns, "nanocount")
        stat_counts = countStats(assignments)
        stats_counts = joinCountStats(stat_counts.map{ it[1]}.collect())
        multiqc_data = multiqc_data.mix(stats_counts)
    }
    else if (params.counting == "htseq" && params.ref_type == "genome") {
        htseq_out = HTSEQ_COUNT(annotation, sorted_alns.join(aln_indexes))
        read_counts = htseq_out.counts
        assignments = AssignReads(htseq_out.bam, "htseq")
        stat_counts = countStats(assignments)
        stats_counts = joinCountStats(stat_counts.map{ it[1]}.collect())
        multiqc_data = multiqc_data.mix(stats_counts)
    } else if (params.counting == "NO") {
    } else {
        println "ERROR ################################################################"
        println "${params.counting} is not compatible with ${params.ref_type}"
        println "htseq requires a genome as reference and an annotation in GTF"
        println "nanocount requires a transcriptome as a reference"     
        println "ERROR ################################################################"
        println "Exiting ..."
        System.exit(0)
    } 
    if (params.discovery == "bambu" && params.ref_type == "genome"){
        sorted_alns.map{
            [it[1]]
        }.collect().map{
            ["assembly", it]
        }.set{data_to_bambu}    
        bambu_out = BAMBU_ASSEMBLE(reference, annotation, data_to_bambu)
    } else if (params.discovery == "isoquant" && params.ref_type == "genome"){
        aln_indexes.map{
            [it[1]]
        }.collect().map{
            ["assembly", it]
        }.set{ixd_4_bambu}
        
        sorted_alns.map{
            [it[1]]
        }.collect().map{
            ["assembly", it]
        }.join(ixd_4_bambu).set{data_to_isoquant}
        //data_to_isoquant.view()
    
        bambu_out = ISOQUANT_ASSEMBLE(reference, annotation, data_to_isoquant)
    } else if (params.discovery == "NO") {
    } else {
        println "ERROR ################################################################"
        println "${params.discovery} is not compatible with ${params.ref_type}"
        println "bambu requires a genome as reference and an annotation in GTF"
        println "ERROR ################################################################"
        println "Exiting ..."
        System.exit(0)
    }
    
    // Perform MULTIQC report
    MULTIQC(multiqc_data.collect())
    
}


workflow preprocess_simple {
    take:
        bc_fastq
        
    main:   

    // Perform Fastqc QC on fastq
    fastqc_files = FASTQC(bc_fastq)

    // Perform mapping on fastq files
    if (params.mapping == "NO") {
        stats_aln = Channel.value() 
        sorted_alns = Channel.value()   
        nanoplot_qcs = Channel.value()  
    }
    else {
        switch(params.mapping) { 
            case "graphmap": 
            dna_bc_fastq = RNA2DNA(bc_fastq)
                        aln_reads = GRAPHMAP(dna_bc_fastq, reference)
            break
            case "graphmap2": 
            aln_reads = GRAPHMAP2(bc_fastq, reference)
            break
                        case "minimap2": 
            aln_reads = MINIMAP2(bc_fastq, reference)
            break
            case "bwa": 
            aln_reads = BWA(reference, bc_fastq)
            break
            default: 
            println "ERROR ################################################################"
            println "${params.mapping} is not a supported alignment"
            println "ERROR ################################################################"
            println "Exiting ..."
            System.exit(0)
            break
        }    

        // Perform SORTING and INDEXING on bam files
        sorted_alns = SAMTOOLS_SORT(aln_reads)
        aln_indexes = SAMTOOLS_INDEX(sorted_alns)

        // Converting BAM to CRAM and 
        if (params.cram_conv == "YES") {
            good_ref = checkRef(reference)
            bam2Cram(good_ref, params.subsampling_cram, sorted_alns.join(aln_indexes))
        }
        
        // Perform bam2stats on sorted bams
        aln_stats = bam2stats(sorted_alns)
        stats_aln = joinAlnStats(aln_stats.map{ it[1]}.collect())
    
        // Perform NanoPlot on sorted bams
        nanoplot_qcs = NANOPLOT_QC(sorted_alns)
    }   


    // OPTIONAL Perform COUNTING / ASSIGNMENT
    if (params.counting == "nanocount" && params.ref_type == "transcriptome") {
        read_counts = NANOCOUNT(sorted_alns.join(aln_indexes))
                //read_counts = NANOCOUNT(sorted_alns)
        assignments = AssignReads(sorted_alns, "nanocount")
        stat_counts = countStats(assignments)
        stats_counts = joinCountStats(stat_counts.map{ it[1]}.collect())
    }
    else if (params.counting == "htseq" && params.ref_type == "genome") {
        htseq_out = HTSEQ_COUNT(params.annotation, sorted_alns.join(aln_indexes))
        read_counts = htseq_out.counts
        assignments = AssignReads(htseq_out.bam, "htseq")
        stat_counts = countStats(assignments)
        stats_counts = joinCountStats(stat_counts.map{ it[1]}.collect())
    } 
    else if (params.counting == "NO") {
        // Default empty channels for reporting
        stats_counts = Channel.value()
    } else {
        println "ERROR ################################################################"
        println "${params.counting} is not compatible with ${params.ref_type}"
        println "htseq requires a genome as reference and an annotation in GTF"
        println "nanocount requires a transcriptome as a reference"     
        println "ERROR ################################################################"
        println "Exiting ..."
        System.exit(0)
    } 
    

    // Perform MULTIQC report
    fastqc_files.map{it[1]}.set{qcs}
    all_res = qcs.mix(multiqc_info,stats_counts, stats_aln)
    MULTIQC(all_res.collect())
    
}


 workflow {
    if (params.fast5 != "" && params.fastq == "") {

        Channel
            .fromPath( params.fast5)                                             
            .ifEmpty { error "Cannot find any file matching: ${params.fast5}" }
            .set {fast5_files}

        fast5_files.map { 
            def filepath = file(it)
            def file_parts = "${filepath}".tokenize("/")
            def folder_name  = filepath[-2]
            [folder_name, it]
        }.groupTuple().set{ fast5_per_folder}
    
        // Check tools
        checkTools(tools, progPars)

        def num = 0
        fast5_per_folder.map{
            def folder_name = it[0]
            def buffer_files = it[1].flatten().collate(params.granularity)
            [folder_name, buffer_files]
        }.transpose().map{
            num++ 
            [ "${it[0]}---${num}", it[1] ]
        }.set{ fast5_4_analysis }

        //GET_WORKFLOWS(params.flowcell, params.kit).view()
        if (params.basecalling == "guppy" && params.demultiplexing == "NO" ) outf = flow1(fast5_4_analysis)
        else outf = flow2(fast5_4_analysis, barcodes_to_include)
        
        def bc_fast5 = outf.basecalled_fast5
        def bc_fastq = outf.basecalled_fastq
        def basecalled_stats = outf.basecalled_stats
        preprocess_flow(bc_fast5, bc_fastq, basecalled_stats)
        
    } else if(params.fast5 == "" && params.fastq != "") {

        // Check tools
        tools["basecalling"] = "NO"
        tools["demultiplexing"] = "NO"
        checkTools(tools, progPars)
        Channel.fromFilePairs( params.fastq , size: 1)
            .ifEmpty { error "Cannot find any file matching: ${params.fastq}" }
            .set {fastq_files}
        
        preprocess_simple(fastq_files)
    
        
    } else {
            println "ERROR ################################################################"
            println "Please choose one between fast5 and fastq as input!!!" 
            println "ERROR ################################################################"
            println "Exiting ..."
            System.exit(0)
        
    }

    //all_ver = BAMBU_VER().mix(DEMULTIPLEX_VER()).mix(NANOQ_VER()).mix(NANOFILT_VER())
    //.mix(GRAPHMAP_VER()).mix(GRAPHMAP2_VER())
    //.mix(MINIMAP2_VER()).mix(BWA_VER()).mix(FASTQC_VER())
    //.mix(SAMTOOLS_VERSION()).mix(NANOPLOT_VER()).mix(NANOCOUNT_VER()).mix(HTSEQ_VER()).mix(MULTIQC_VER())
    //.collectFile(name: 'tool_version.txt', newLine: false, storeDir:outputMultiQC)
    
 
 }

workflow.onComplete {
    println "Pipeline BIOCORE@CRG Master of Pore - preprocess completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

/*
* Mail notification
*/

if (params.email == "yourmail@yourdomain" || params.email == "") { 
    log.info 'Skipping the email\n'
}
else {
    log.info "Sending the email to ${params.email}\n"

    workflow.onComplete {

    def msg = """\
        Pipeline BIOCORE@CRG Master of Pore 2 preprocess execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
        .stripIndent()

        sendMail(to: params.email, subject: "Master of Pore 3 execution", body: msg, attach: "${outputMultiQC}/multiqc_report.html")
    }
}


