#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastq_dir = "${projectDir}/../fastq"
params.ref_dir   = "${projectDir}/../reference"
params.outdir    = "${projectDir}/../results"
params.threads   = 4

process ALIGN_AND_QC {
    tag "$sample ($type)"
    cpus params.threads
    publishDir "${params.outdir}/bam", mode: 'copy'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.txt'

    input:
    tuple val(sample), val(type), path(reads)
    path genome

    output:
    path "${sample}.sorted.bam", emit: bam
    path "${sample}.sorted.bam.bai", emit: bai
    path "${sample}.flagstat.txt", emit: flagstat
    path "${sample}.stats.txt", emit: stats

    script:
    def preset = type == 'cdna' ? '-ax splice:hq' : '-ax splice -uf -k14'
    """
    minimap2 -t ${task.cpus} ${preset} "$genome" "$reads" |
        samtools sort -@ ${task.cpus} -o "${sample}.sorted.bam" -
    samtools index "${sample}.sorted.bam"
    samtools flagstat "${sample}.sorted.bam" > "${sample}.flagstat.txt"
    samtools stats "${sample}.sorted.bam" > "${sample}.stats.txt"
    """
}

process BAMBU {
    tag 'Bambu discovery and quantification with reference annotations'
    cpus 4
    publishDir "${params.outdir}/bambu", mode: 'copy'

    input:
    path bams
    path genome
    path gtf

    output:
    path "counts_transcript.txt"
    path "counts_gene.txt"
    path "extended_annotations.gtf"

    script:
    """
    cat > run_bambu.R <<'RSCRIPT'
    suppressPackageStartupMessages(library(bambu))

    # Compatibility patch for Bambu 3.4.1 with dplyr >= 1.1.
    ns <- asNamespace("bambu")
    f <- get("makeUnsplicedTibble", envir = ns)
    body_text <- paste(deparse(body(f)), collapse = "\\n")
    old <- "group_by\\\\(chr, strand, start,\\\\s*end\\\\) %>% summarise\\\\("
    if (grepl(old, body_text, perl = TRUE)) {
        body_text <- sub(old, "group_by(chr, strand, start, end) %>% reframe(", body_text, perl = TRUE)
        body(f) <- parse(text = body_text)[[1L]]
        unlockBinding("makeUnsplicedTibble", ns)
        assign("makeUnsplicedTibble", f, envir = ns)
        lockBinding("makeUnsplicedTibble", ns)
    }

    staged_files <- list.files(".", full.names = TRUE)
    reads <- staged_files[grepl("[.]bam", staged_files)]
    genome_file <- staged_files[grepl("[.]fa\$", staged_files)][1]
    annotation_file <- staged_files[grepl("[.]gtf\$", staged_files)][1]
    annotations <- prepareAnnotations(annotation_file)
    se <- bambu(
        reads = reads,
        annotations = annotations,
        genome = genome_file,
        discovery = TRUE,
        quant = TRUE,
        ncore = 4,
        lowMemory = TRUE,
        verbose = TRUE
    )
    writeBambuOutput(se, path = "./")
    RSCRIPT
    Rscript run_bambu.R
    test -s counts_transcript.txt
    test -s counts_gene.txt
    test -s extended_annotations.gtf
    """
}

workflow {
    fastqs = channel.fromPath("${params.fastq_dir}/*")
        .filter { it.name ==~ /(?i).*\.(fastq|fq)(\.gz)?$/ }
        .map { fq ->
            def lower = fq.name.toLowerCase()
            def type = lower.contains('directrna') || lower.contains('direct_rna') ? 'direct_rna' :
                       (lower.contains('cdna') ? 'cdna' : null)
            if (type == null) error "FASTQ filename must contain cDNA or direct_RNA: ${fq.name}"
            def sample = fq.name.replaceFirst(/(?i)\.(fastq|fq)(\.gz)?$/, '')
            tuple(sample, type, fq)
        }

    genome = channel.value(file("${params.ref_dir}/grch38.fa"))
    gtf = channel.value(file("${params.ref_dir}/Homo_sapiens.GRCh38.91.gtf"))

    aligned = ALIGN_AND_QC(fastqs, genome)
    BAMBU(aligned.bam.collect(), genome, gtf)
}
