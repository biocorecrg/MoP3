#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* 
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '3.0'

params.help            = false
params.resume          = false

log.info """

╔╦╗╔═╗╔═╗  ╔═╗┬─┐┌─┐┌─┐┬─┐┌─┐┌─┐┌─┐┌─┐┌─┐
║║║║ ║╠═╝  ╠═╝├┬┘├┤ ├─┘├┬┘│ ││  ├┤ └─┐└─┐
╩ ╩╚═╝╩    ╩  ┴└─└─┘┴  ┴└─└─┘└─┘└─┘└─┘└─┘
                                                                                       
====================================================
BIOCORE@CRG Master of Pores 3. Preprocessing - N F  ~  version ${version}
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

email                     : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// include functions, outdirs from other files
evaluate(new File("../outdirs.nf"))
def local_modules = file("${projectDir}/../local_modules.nf")
def subworkflowsDir = "${projectDir}/../BioNextflow/subworkflows"
joinScript = file("${projectDir}/bin/join.r")

// get and check input files
if (params.mapping != "NO") {
    reference = file(params.reference)
    if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
} else {
    reference = ""
}

// INIZIALIZE MULTIQC REPORT
config_report = file("${projectDir}/config.yaml")
if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("${projectDir}/../img/logo_small.png")

Channel.from( config_report, logo )
    .collect().set{multiqc_info}

outputReport   = file("${outputMultiQC}/multiqc_report.html")

if( outputReport.exists() ) {
  log.info "Moving old report to multiqc_report.html multiqc_report.html.old"
  outputReport.moveTo("${outputMultiQC}/multiqc_report.html.old")
}

// Get models
Channel.fromPath( "${projectDir}/deeplexicon/*.h5").set{deepmodels}
Channel.fromPath( "${projectDir}/dorado_models/*", type: 'dir').collect().set{dorado_models}
Channel.fromPath( "${projectDir}/seqtagger_models").set{seqtagger_models}

dorado_models.view()

// check GPU usage. 
if (params.GPU != "cuda11" && params.GPU != "cuda10" && params.GPU != "OFF" && params.GPU != "ON") exit 1, "Please specify cuda11, cuda10, ON or OFF if GPU processors are available. ON is legacy for cuda10"
def gpu = (params.GPU != 'OFF' ? 'ON' : 'OFF')
def cuda_cont = (params.GPU == 'cuda11' ? 'biocorecrg/mopbasecallc11:0.3' : 'biocorecrg/mopbasecall:0.3')



// CHECK INCOMPATIBILITIES AMONG PARAMETERS

if (params.ref_type == "genome") {
    if (params.annotation != "") {
        annotation = file(params.annotation)
        if( !annotation.exists() ) exit 1, "Missing annotation file: ${params.annotation}!"
    }
}

outmode = "copy"

include { get_barcode_list; RNA2DNA; preparing_demultiplexing_fast5_seqtagger; preparing_demultiplexing_fast5_deeplexicon; extract_seqtagger_fastq; extract_deeplexicon_fastq; parseFinalSummary; checkTools; reshapeSamples; reshapeDemuxSamples; checkRef; getParameters; homogenizeVals } from "${local_modules}" 
include { extract_demultiplexed_fast5 as extract_demultiplexed_fast5_seqtagger;  extract_demultiplexed_fast5 as extract_demultiplexed_fast5_deeplexicon } from "${local_modules}" addParams(OUTPUTF5: outputFast5, OUTPUTST: outputQual, LABEL: 'big_cpus')
include { extract_demultiplexed_fast5_guppy } from "${local_modules}" addParams(OUTPUT: outputFast5, LABEL: 'big_cpus')

def demulti_fast5_opt = homogenizeVals(params.demulti_fast5)

def guppy_basecall_label = (params.GPU != 'OFF' ? 'basecall_gpus' : 'big_cpus')
def deeplexi_basecall_label = (params.GPU != 'OFF' ? 'demulti_gpus' : '')
def output_bc = (demulti_fast5_opt == 'ON' ? '' : outputFast5)
def outputMinionQC = (demulti_fast5_opt == 'ON' ? '': outputQual)


def guppypars = ""
// GET PROGRAM PARS AND VERIFY 
def tools = [:]
tools["basecalling"] = homogenizeVals(params.basecalling)
tools["demultiplexing"] = homogenizeVals(params.demultiplexing)
tools["mapping"] = homogenizeVals(params.mapping)
tools["filtering"] = homogenizeVals(params.filtering)
tools["counting"] = homogenizeVals(params.counting)
tools["discovery"] = homogenizeVals(params.discovery)

// Remove basecalling and demultiplexing in case of fastq input
if(params.fast5 == "" && params.fastq != "") {
    tools["basecalling"] = "NO"
    tools["demultiplexing"] = "NO"
} else {
    guppypars = parseFinalSummary(params.conffile)
    // Create a channel for tool options
    if (workflow.profile == "awsbatch") guppypars = guppypars + " --data_path /nextflow-bin/ont-guppy/data"
}

progPars = getParameters(params.pars_tools)
checkTools(tools, progPars)

// Create a channel for excluded ids
barcodes_to_include = get_barcode_list(params.barcodes)

// CHECK GUPPY VERSION
include { GET_VERSION as GUPPY_VERSION } from "${subworkflowsDir}/basecalling/guppy" 

def guppy_basecall_pars = guppypars + " " + progPars["basecalling--guppy"]
def guppy6_basecall_pars = guppy_basecall_pars + "--disable_qscore_filtering"


// INCLUDE MODULES / SUBWORKFLOWS
include { final_message; notify_slack } from "${subworkflowsDir}/global_functions.nf"
include { GET_WORKFLOWS; BASECALL as GUPPY_BASECALL; BASECALL_DEMULTI as GUPPY_BASECALL_DEMULTI } from "${subworkflowsDir}/basecalling/guppy" addParams(EXTRAPARS_BC: guppy_basecall_pars, EXTRAPARS_DEM: progPars["demultiplexing--guppy"], LABEL: guppy_basecall_label, GPU: gpu, MOP: "YES", OUTPUT: output_bc, CONTAINER: cuda_cont, OUTPUTMODE: outmode)
include { BASECALL as GUPPY6_BASECALL; BASECALL_DEMULTI as GUPPY6_BASECALL_DEMULTI } from "${subworkflowsDir}/basecalling/guppy" addParams(EXTRAPARS_BC: guppy6_basecall_pars, EXTRAPARS_DEM: progPars["demultiplexing--guppy"], LABEL: guppy_basecall_label, GPU: gpu, MOP: "YES", OUTPUT: output_bc, CONTAINER: cuda_cont, OUTPUTMODE: outmode)
include { BASECALL as DORADO_BASECALL } from "${subworkflowsDir}/basecalling/dorado" addParams(EXTRAPARS: progPars["basecalling--dorado"], LABEL: guppy_basecall_label, GPU: gpu, MOP: "YES", OUTPUT: output_bc, OUTPUTMODE: outmode)
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
include { getFast5; filterBarcodes } from "${local_modules}" 


// ADD A CHECK FOR GUPPY FOR DISABLING SCORE

def separateGuppy (fast5) {

	data_and_ver = GUPPY_VERSION().map{
        def vals = it.split("\\.")
    	vals[0] 
    }.toInteger().combine(fast5)
		
    newer = data_and_ver.map{ if (it[0] >= 6) [ it[1], it[2] ]}
    older = data_and_ver.map{ if (it[0] < 6 ) [ it[1], it[2] ]}
 
    return([newer, older])   
}
    

/*
* Wrapper for basecalling
*/
workflow BASECALL {
    take: 
        fast5_4_analysis
        
    main:
        switch(params.basecalling) {                      
           case "guppy":
               (newer, older) = separateGuppy(fast5_4_analysis)    
               outbc = GUPPY6_BASECALL(newer)
               outbc6 = GUPPY_BASECALL(older)  
               
               basecalled_fastq = outbc.basecalled_fastq.concat(outbc6.basecalled_fastq)    
               basecalled_fast5 = outbc.basecalled_fast5.concat(outbc6.basecalled_fast5)    
               basecalling_stats = outbc.basecalling_stats.concat(outbc6.basecalling_stats)    
 
               bc_fast5 = reshapeSamples(basecalled_fast5)
               bc_stats = reshapeSamples(basecalling_stats)
               break; 
           case "dorado": 
       	       outbc = DORADO_BASECALL (fast5_4_analysis, dorado_models)
               basecalled_fastq = outbc.basecalled_fastq
               bc_fast5 = Channel.empty()
               bc_stats = Channel.empty()
               break; 
        }        
		
    emit:
       basecalled_fast5 = bc_fast5
       basecalled_fastq = basecalled_fastq
       basecalling_stats = bc_stats
}

/*
* Wrapper for demultiplexing
*/
workflow DEMULTIPLEX {
    take: 
        fast5_4_analysis
        
    main:
        switch(params.demultiplexing) {                      
           case "deeplexicon":
               outbc = BASECALL(fast5_4_analysis)
               demux = DEMULTIPLEX_DEEPLEXICON(deepmodels, fast5_4_analysis)
               demufq = extract_deeplexicon_fastq(demux.join(outbc.basecalled_fastq))
               basecalling_stats = outbc.basecalling_stats    
               basecalled_fast5 = outbc.basecalled_fast5 
               break;
           case "seqtagger": 
               outbc = BASECALL(fast5_4_analysis)
               demux = SEQTAGGER_DEMULTIPLEX(fast5_4_analysis, seqtagger_models)
               demufq = extract_seqtagger_fastq(demux.join(outbc.basecalled_fastq))
               basecalling_stats = outbc.basecalling_stats    
               basecalled_fast5 = outbc.basecalled_fast5 
               break;
           default:
               (newer, older) = separateGuppy(fast5_4_analysis)   
               outbc = GUPPY_BASECALL_DEMULTI(older)
               outbc6 = GUPPY6_BASECALL_DEMULTI(newer)

               demufq = outbc.basecalled_fastq.concat(outbc6.basecalled_fastq)    
               basecalled_fast5 = outbc.basecalled_fast5.concat(outbc6.basecalled_fast5)    
               basecalling_stats = outbc.basecalling_stats.concat(outbc6.basecalling_stats)    
               if (params.demultiplexing == "readucks") {
                    demufq = READUCKS_DEMULTIPLEX(demufq)
               }
               break;
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

    emit:
        demultiplexed_fastqs =  demufq
        basecalled_fastq = reshapedDemufq
        basecalling_stats = basecalling_stats
        basecalled_fast5 = basecalled_fast5
}

workflow preprocess_flow {
    take:
        bc_fast5
        bc_fastq
        basecalling_stats
        
    main:   
    // Perform MinIONQC on basecalling stats
    basecall_qc = MinIONQC(basecalling_stats.groupTuple())   
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

        // Concatenate bamfiles differently depending on if demultiplexed or not
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
        // Perform bam2stats on sorted bams
        aln_stats = bam2stats(sorted_alns)
        stats_aln = joinAlnStats(aln_stats.map{ it[1]}.collect())
    
        // Perform NanoPlot on sorted bams
        nanoplot_qcs = NANOPLOT_QC(sorted_alns)
        multiqc_data = multiqc_data.mix(stats_aln)
    }   

    // Concatenate fastq files differently depending on if demultiplexed or not
    if (params.demultiplexing == "NO" ) reshaped_bc_fastq = reshapeSamples(bc_fastq)
    else reshaped_bc_fastq = reshapeDemuxSamples(bc_fastq)
    fastq_files = concatenateFastQFiles(reshaped_bc_fastq.groupTuple())

    // Perform fastqc QC on fastq
    fastqc_files = FASTQC(fastq_files)
    multiqc_data = multiqc_data.mix(fastqc_files.map{it[1]})

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
 
    // INPUT IS FAST5
    if (params.fast5 != "" && params.fastq == "") {
        fast5_4_analysis = getFast5(params.fast5)
        
        if (params.demultiplexing == "NO" ) outf = BASECALL(fast5_4_analysis)
        else outf = DEMULTIPLEX(fast5_4_analysis)

       // Optional fastq filtering 
        switch(params.filtering) {                      
          case "nanofilt": 
            basecalled_fastq = NANOFILT_FILTER(outf.basecalled_fastq)
            break; 
          case "nanoq": 
            basecalled_fastq = NANOQ_FILTER(outf.basecalled_fastq)
            break;  
          default:
            basecalled_fastq = outf.basecalled_fastq
            break;         
        }   
        
        preprocess_flow(outf.basecalled_fast5, basecalled_fastq, outf.basecalling_stats)
        
    } else if(params.fast5 == "" && params.fastq != "") {
    // INPUT IS FASTQ
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

    def text = final_message("MoP3")
    println text
    if (params.hook != "") {
       notify_slack(text, params.hook)
    }
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
 	   def msg = final_message("MoP3")	
        sendMail(to: params.email, subject: "MoP3 - preprocess execution", body: msg, attach: "${outputMultiQC}/multiqc_report.html")
    }
}
