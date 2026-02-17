#!/usr/bin/env nextflow
/*
══════════════════════════════════════════════════════════════════════════════════
    SPLIT — SNP-Level Inspection of Parental Transcripts
══════════════════════════════════════════════════════════════════════════════════
    Allele-specific RNA-seq pipeline  •  Nextflow DSL2
    Ultra-minimalist design for solo bioinformaticians.

    Steps:
      0a. [Optional] GEO → SRR resolution  (NCBI E-utilities, Python stdlib)
      0b. [Optional] SRA download           (fasterq-dump + bgzip)
      1.  Download references               (genome FASTA, GTF, VCF from Ensembl/MGP)
      2.  SNPsplit genome preparation        (N-masked genome for dual-hybrid)
      3.  STAR index                         (N-masked + reference genomes)
      4.  STAR alignment                     (both N-masked and reference tracks)
      5.  Sort + optional deduplication      (samtools sort + markdup -r)
      6.  SNPsplit allele separation         (on N-masked BAMs only)
      7.  featureCounts                      (genome1, genome2, reference — gene_name)

    Input modes:
      --sra_ids     SRR/ERR/DRR accessions or GSE/GSM (comma-separated or file)
      --fastq_dir   Directory of *.fastq.gz (single-end or paired-end auto-detect)
      --input       CSV samplesheet: sample,fastq_1,fastq_2

    Dual-hybrid allele-specific analysis:
      --strain1     First strain in VCF  (default: CAST_EiJ → genome1)
      --strain2     Second strain in VCF (default: C57BL_6NJ → genome2)
══════════════════════════════════════════════════════════════════════════════════
*/

nextflow.enable.dsl = 2

// ── Print banner ────────────────────────────────────────────────────────────

log.info """
╔══════════════════════════════════════════════════════════════════╗
║   S P L I T                                                      ║
║   SNP-Level Inspection of Parental Transcripts                   ║
║   v1.0.0                                                       ║
╚══════════════════════════════════════════════════════════════════╝
Input        : ${params.input ?: params.fastq_dir ?: params.sra_ids ?: 'NONE'}
Output       : ${params.outdir}
Strain 1     : ${params.strain1} (→ genome1)
Strain 2     : ${params.strain2} (→ genome2)
Dedup        : ${params.dedup}
Strandedness : ${params.strandedness}
CPUs         : ${params.max_cpus}
──────────────────────────────────────────────────────────────────
""".stripIndent()

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 0a : GEO → SRR resolution (NCBI E-utilities, Python stdlib)
//  Accepts GSE or GSM accessions and returns a list of SRR IDs.
// ─────────────────────────────────────────────────────────────────────────────

process RESOLVE_GEO {
    tag "${geo_id}"
    label 'process_low'
    errorStrategy 'retry'
    maxRetries 2

    input:
    val geo_id

    output:
    stdout

    script:
    """
    #!/usr/bin/env python3
    import urllib.request, xml.etree.ElementTree as ET, sys, time

    def esearch(db, term):
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={db}&term={term}&retmax=10000&usehistory=y"
        for attempt in range(3):
            try:
                return ET.parse(urllib.request.urlopen(url, timeout=30))
            except:
                time.sleep(5 * (attempt + 1))
        sys.exit(f"Failed to query NCBI for {term}")

    def elink(dbfrom, db, ids):
        id_str = ",".join(ids)
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom={dbfrom}&db={db}&id={id_str}"
        for attempt in range(3):
            try:
                return ET.parse(urllib.request.urlopen(url, timeout=30))
            except:
                time.sleep(5 * (attempt + 1))
        sys.exit(f"Failed elink for {ids}")

    geo_id = "${geo_id}".strip()
    srr_ids = []

    # GSE → GDS → SRA
    if geo_id.startswith("GSE"):
        tree = esearch("gds", geo_id + "[ACCN]")
        gds_ids = [e.text for e in tree.findall(".//Id")]
        if not gds_ids:
            sys.exit(f"No GDS entries for {geo_id}")
        tree = elink("gds", "sra", gds_ids)
        sra_ids = [e.text for e in tree.findall(".//Link/Id")]
        if sra_ids:
            from urllib.request import urlopen
            id_str = ",".join(sra_ids)
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={id_str}&rettype=runinfo&retmode=text"
            data = urlopen(url, timeout=60).read().decode()
            for line in data.strip().split("\\n")[1:]:
                fields = line.split(",")
                if len(fields) > 0 and fields[0].startswith(("SRR","ERR","DRR")):
                    srr_ids.append(fields[0])

    # GSM → SRA
    elif geo_id.startswith("GSM"):
        tree = esearch("gds", geo_id + "[ACCN]")
        gds_ids = [e.text for e in tree.findall(".//Id")]
        if gds_ids:
            tree = elink("gds", "sra", gds_ids)
            sra_ids = [e.text for e in tree.findall(".//Link/Id")]
            if sra_ids:
                from urllib.request import urlopen
                id_str = ",".join(sra_ids)
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={id_str}&rettype=runinfo&retmode=text"
                data = urlopen(url, timeout=60).read().decode()
                for line in data.strip().split("\\n")[1:]:
                    fields = line.split(",")
                    if len(fields) > 0 and fields[0].startswith(("SRR","ERR","DRR")):
                        srr_ids.append(fields[0])

    if not srr_ids:
        sys.exit(f"No SRR accessions found for {geo_id}")

    print(",".join(sorted(set(srr_ids))), end="")
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 0b : SRA download (fasterq-dump + bgzip)
//  Downloads a single SRR accession and compresses to .fastq.gz
// ─────────────────────────────────────────────────────────────────────────────

process SRA_DOWNLOAD {
    tag "${srr_id}"
    label 'process_medium'
    publishDir "${params.outdir}/00_sra_fastq", mode: 'copy'

    input:
    val srr_id

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads

    script:
    meta = [id: srr_id]
    """
    # Download SRA and split into separate read files
    fasterq-dump ${srr_id} --threads ${task.cpus} --split-files --progress

    # Compress all resulting FASTQ files with bgzip
    for fq in *.fastq; do
        bgzip -@ ${task.cpus} -f "\${fq}"
    done

    # Auto-detect SE vs PE based on output file count
    COUNT=\$(ls -1 *.fastq.gz | wc -l)
    if [ "\${COUNT}" -eq 1 ]; then
        mv *.fastq.gz ${srr_id}.fastq.gz
    elif [ "\${COUNT}" -ge 2 ]; then
        # Rename to standard R1/R2 naming if needed
        if [ -f "${srr_id}_1.fastq.gz" ]; then
            mv ${srr_id}_1.fastq.gz ${srr_id}_R1_001.fastq.gz
            mv ${srr_id}_2.fastq.gz ${srr_id}_R2_001.fastq.gz
        fi
    fi
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 1 : Download references (genome FASTA, GTF, VCF)
//  Downloads from Ensembl (genome + GTF) and MGP (VCF) for the specified
//  species/assembly. Uses storeDir for persistent caching.
// ─────────────────────────────────────────────────────────────────────────────

process DOWNLOAD_REFERENCES {
    label 'process_low'
    storeDir "${params.outdir}/reference"

    output:
    path "genome.fa",              emit: genome_fa
    path "genome.fa.fai",          emit: genome_fai
    path "genes.gtf",              emit: gtf
    path "snps.vcf.gz",            emit: vcf
    path "snps.vcf.gz.csi",        emit: vcf_index

    script:
    """
    # ── Genome FASTA ──
    echo "[DL] Genome FASTA: ${params.genome_url}"
    wget -q -O genome.fa.gz "${params.genome_url}"
    gunzip -f genome.fa.gz
    samtools faidx genome.fa

    # ── GTF annotation ──
    echo "[DL] GTF: ${params.gtf_url}"
    wget -q -O genes.gtf.gz "${params.gtf_url}"
    gunzip -f genes.gtf.gz

    # ── VCF (MGP SNPs) ──
    echo "[DL] VCF: ${params.vcf_url}"
    wget -q -O snps.vcf.gz "${params.vcf_url}"

    echo "[DL] VCF index: ${params.vcf_index_url}"
    wget -q -O snps.vcf.gz.csi "${params.vcf_index_url}"
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 2 : SNPsplit genome preparation (dual-hybrid N-masking)
//  Creates an N-masked genome and a universal SNP file for allele assignment.
//  --dual_hybrid mode: both strains are processed, genome positions with
//  strain-discriminating SNPs are replaced by N in the output FASTA.
// ─────────────────────────────────────────────────────────────────────────────

process SNPSPLIT_GENOME_PREP {
    label 'process_high'
    storeDir "${params.outdir}/reference/snpsplit_prep"

    input:
    path genome_fa
    path genome_fai
    path vcf
    path vcf_index

    output:
    path "genome.N-masked.fa",     emit: nmask_fa
    path "genome.N-masked.fa.fai", emit: nmask_fai
    path "all_SNPs_*.txt",         emit: snp_file

    script:
    """
    # Create a directory structure that SNPsplit expects
    mkdir -p ref_genome
    cp ${genome_fa} ref_genome/
    cp ${genome_fai} ref_genome/

    # Run SNPsplit genome preparation in dual-hybrid mode
    # This generates per-chromosome N-masked FASTA files and a SNP annotation file
    SNPsplit_genome_preparation \\
        --reference_genome ref_genome \\
        --vcf_file ${vcf} \\
        --strain "${params.strain1}" \\
        --strain2 "${params.strain2}" \\
        --dual_hybrid \\
        --nmasking

    # Find the dual hybrid output directory
    DUAL_DIR=\$(find . -maxdepth 1 -type d -name "*dual_hybrid*" | head -1)

    if [ -z "\${DUAL_DIR}" ]; then
        echo "ERROR: dual_hybrid directory not found after SNPsplit_genome_preparation"
        exit 1
    fi

    # Concatenate per-chromosome N-masked FASTAs into a single genome file
    cat \${DUAL_DIR}/*.fa > genome.N-masked.fa
    samtools faidx genome.N-masked.fa

    # Locate the SNP annotation file (used by SNPsplit for allele assignment)
    SNP_FILE=\$(find . -maxdepth 1 -name "all_*SNPs_*reference.based_on_*.txt" | head -1)
    if [ -z "\${SNP_FILE}" ]; then
        echo "ERROR: SNP file not found after SNPsplit_genome_preparation"
        exit 1
    fi
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 3 : STAR genomeGenerate (index building)
//  Builds a STAR genome index. Called twice: once for N-masked, once for ref.
//  sjdbOverhang is auto-estimated from the first FASTQ file.
// ─────────────────────────────────────────────────────────────────────────────

process STAR_INDEX {
    tag "${index_name}"
    label 'process_high'
    storeDir "${params.outdir}/reference/star_${index_name}"

    input:
    val index_name
    path genome_fa
    path gtf
    val sjdb_overhang

    output:
    path "star_idx", emit: index

    script:
    """
    mkdir -p star_idx

    STAR \\
        --runThreadN ${task.cpus} \\
        --runMode genomeGenerate \\
        --genomeDir star_idx \\
        --genomeFastaFiles ${genome_fa} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang ${sjdb_overhang} \\
        --genomeSAsparseD 2 \\
        --limitGenomeGenerateRAM ${params.star_limit_genome_ram}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 3b : Estimate read length from a FASTQ file
//  Reads the first 50,000 lines to determine the most common read length.
//  Returns sjdbOverhang = read_length - 1.
// ─────────────────────────────────────────────────────────────────────────────

process ESTIMATE_READ_LENGTH {
    label 'process_low'

    input:
    path fastq

    output:
    stdout

    script:
    """
    # Sample the first 50000 lines (12500 reads) and find the most common length
    gzip -cd ${fastq} | head -n 50000 | \\
        awk 'NR%4==2 { l[length(\$0)]++ } END { max=0; best=0; for(k in l) { if(l[k]>max) { max=l[k]; best=k } } print best - 1 }'
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 4 : STAR alignment
//  Aligns reads to a STAR index. Works for both SE and PE.
//  --alignEndsType EndToEnd: no soft-clipping (recommended for SNPsplit)
//  --outFilterMultimapNmax 1: unique mappers only
// ─────────────────────────────────────────────────────────────────────────────

process STAR_ALIGN {
    tag "${meta.id}_${track}"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path star_index
    val track   // "nmask" or "ref"

    output:
    tuple val(meta), val(track), path("*.Aligned.out.bam"), emit: bam
    path "*Log.final.out",                                  emit: log

    script:
    def prefix = "${meta.id}_${track}_"
    def read_files = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${read_files} \\
        --readFilesCommand gzip -cd \\
        --alignEndsType EndToEnd \\
        --outSAMattributes NH HI NM MD \\
        --outFilterMultimapNmax 1 \\
        --outSAMtype BAM Unsorted \\
        --outFileNamePrefix ${prefix}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 5 : Sort + optional deduplication (samtools sort + markdup -r)
//  If --dedup is true, removes PCR/optical duplicates after coordinate sort.
//  Produces an indexed BAM ready for downstream analysis.
// ─────────────────────────────────────────────────────────────────────────────

process SORT_DEDUP {
    tag "${meta.id}_${track}"
    label 'process_medium'
    publishDir "${params.outdir}/${track == 'nmask' ? '04_aln_nmask' : '05_aln_ref'}/${meta.id}", mode: 'copy', pattern: "*.{bam,bai}"

    input:
    tuple val(meta), val(track), path(bam)

    output:
    tuple val(meta), val(track), path("*.sorted*.bam"), path("*.sorted*.bam.bai"), emit: bam

    script:
    def prefix = "${meta.id}_${track}"
    if (params.dedup)
        """
        # Coordinate sort
        samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${bam}

        # Remove PCR/optical duplicates (-r removes them, not just marks)
        samtools markdup -@ ${task.cpus} -r ${prefix}.sorted.bam ${prefix}.sorted.dedup.bam
        samtools index ${prefix}.sorted.dedup.bam

        # Clean up intermediate
        rm -f ${prefix}.sorted.bam
        """
    else
        """
        samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${bam}
        samtools index ${prefix}.sorted.bam
        """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 6 : SNPsplit (allele separation)
//  Splits N-mask–aligned BAMs into genome1 (strain1) and genome2 (strain2).
//  Operates on the deduplicated BAM from the N-mask track.
// ─────────────────────────────────────────────────────────────────────────────

process SNPSPLIT {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/06_snpsplit/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), val(track), path(bam), path(bai)
    path snp_file

    output:
    tuple val(meta), path("*.genome1.bam"),     emit: genome1_bam
    tuple val(meta), path("*.genome2.bam"),     emit: genome2_bam
    tuple val(meta), path("*.unassigned.bam"),  emit: unassigned_bam
    path "*.SNPsplit_report.yaml",              emit: report, optional: true
    path "*SNPsplit_report.txt",                emit: report_txt, optional: true
    path "*SNPsplit_sort.txt",                  emit: sort_report, optional: true

    script:
    def se_flag = meta.single_end ? "--single_end" : "--paired"
    """
    SNPsplit \\
        ${se_flag} \\
        --snp_file ${snp_file} \\
        ${bam}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 7 : featureCounts (read quantification per gene)
//  Counts reads overlapping gene features. Uses gene_name for human-readable
//  output (e.g., Gapdh instead of ENSMUSG00000057666).
//  -s strandedness: 0 = unstranded, 1 = forward, 2 = reverse
// ─────────────────────────────────────────────────────────────────────────────

process FEATURECOUNTS {
    tag "${count_label}"
    label 'process_medium'
    publishDir "${params.outdir}/07_counts", mode: 'copy'

    input:
    path bam_files
    path gtf
    val count_label   // e.g., "genome1_CAST_EiJ", "genome2_C57BL_6NJ", "reference"

    output:
    path "counts_${count_label}.txt",         emit: counts
    path "counts_${count_label}.txt.summary", emit: summary

    script:
    def pe_flag = params.force_se ? "" : "-p --countReadPairs"
    """
    featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -s ${params.strandedness} \\
        -g gene_name \\
        ${pe_flag} \\
        -o counts_${count_label}.txt \\
        ${bam_files}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS 8 : MultiQC (aggregate all QC reports)
// ─────────────────────────────────────────────────────────────────────────────

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/08_multiqc", mode: 'copy'

    input:
    path ('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    """
    multiqc . --force --no-data-dir 2>/dev/null || multiqc . --force
    """
}


// ═══════════════════════════════════════════════════════════════════════════════
//  WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════════

workflow {

    // ── Build input channel ───────────────────────────────────────────────
    // Three mutually exclusive input modes (same as STREAM)

    if (params.sra_ids) {

        // Parse SRA IDs: can be a file (one per line) or comma-separated string
        def raw_ids
        if (file(params.sra_ids).exists()) {
            raw_ids = Channel.fromPath(params.sra_ids)
                .splitText()
                .map { it.trim() }
                .filter { it }
        } else {
            raw_ids = Channel.of(params.sra_ids.toString().split(','))
                .flatten()
                .map { it.trim() }
                .filter { it }
        }

        // Separate GEO accessions (need resolution) from direct SRR/ERR/DRR
        def geo_ids   = raw_ids.filter { it.startsWith('GSE') || it.startsWith('GSM') }
        def direct_ids = raw_ids.filter { it.startsWith('SRR') || it.startsWith('ERR') || it.startsWith('DRR') }

        // Resolve GEO → SRR, then merge with direct IDs
        resolved = RESOLVE_GEO(geo_ids)
            .splitCsv()
            .flatten()
            .map { it.trim() }
            .filter { it }

        all_srr = direct_ids.mix(resolved)

        // Download each SRR
        ch_reads = SRA_DOWNLOAD(all_srr).reads

    } else if (params.fastq_dir) {

        // Auto-detect PE (matching _R1/_R2) vs SE (everything else)
        ch_pe = Channel.fromFilePairs("${params.fastq_dir}/*_R{1,2}_001.fastq.gz", flat: false)
            .map { id, files -> [[id: id, single_end: false], files] }

        ch_se = Channel.fromPath("${params.fastq_dir}/*.fastq.gz")
            .filter { !it.getName().matches(/.*_R[12]_001\.fastq\.gz/) }
            .map { fq -> [[id: fq.getBaseName().replaceAll(/\.fastq$/, ''), single_end: true], [fq]] }

        ch_reads = ch_pe.mix(ch_se)

    } else if (params.input) {

        ch_reads = Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample, single_end: !row.fastq_2]
                def files = row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)]
                [meta, files]
            }

    } else {
        error "Please specify --sra_ids, --fastq_dir, or --input"
    }

    // ── 1. Download references ──────────────────────────────────────────
    DOWNLOAD_REFERENCES()

    // ── 2. SNPsplit genome preparation (N-masked genome) ────────────────
    SNPSPLIT_GENOME_PREP(
        DOWNLOAD_REFERENCES.out.genome_fa,
        DOWNLOAD_REFERENCES.out.genome_fai,
        DOWNLOAD_REFERENCES.out.vcf,
        DOWNLOAD_REFERENCES.out.vcf_index
    )

    // ── 3. Estimate read length from first FASTQ ────────────────────────
    first_fq = ch_reads.first().map { meta, reads -> reads[0] }
    sjdb_overhang = ESTIMATE_READ_LENGTH(first_fq).map { it.trim() as Integer }

    // ── 3b. Build STAR indexes (N-masked + reference) ───────────────────
    star_idx_nmask = STAR_INDEX(
        'nmask',
        SNPSPLIT_GENOME_PREP.out.nmask_fa,
        DOWNLOAD_REFERENCES.out.gtf,
        sjdb_overhang
    ).index

    star_idx_ref = STAR_INDEX(
        'ref',
        DOWNLOAD_REFERENCES.out.genome_fa,
        DOWNLOAD_REFERENCES.out.gtf,
        sjdb_overhang
    ).index

    // ── 4. STAR alignment — two parallel tracks ─────────────────────────
    // Track 1: N-masked (for SNPsplit allele separation)
    bam_nmask = STAR_ALIGN(ch_reads, star_idx_nmask, 'nmask')

    // Track 2: Reference (for standard gene counting)
    bam_ref   = STAR_ALIGN(ch_reads, star_idx_ref,   'ref')

    // ── 5. Sort + optional deduplication ────────────────────────────────
    sorted_nmask = SORT_DEDUP(bam_nmask.bam)
    sorted_ref   = SORT_DEDUP(bam_ref.bam)

    // ── 6. SNPsplit (allele separation on N-mask BAMs) ──────────────────
    SNPSPLIT(
        sorted_nmask.bam,
        SNPSPLIT_GENOME_PREP.out.snp_file
    )

    // ── 7. featureCounts ────────────────────────────────────────────────
    // Collect all genome1 BAMs → count table for strain1
    genome1_bams = SNPSPLIT.out.genome1_bam.map { meta, bam -> bam }.collect()
    genome2_bams = SNPSPLIT.out.genome2_bam.map { meta, bam -> bam }.collect()
    ref_bams     = sorted_ref.bam.map { meta, track, bam, bai -> bam }.collect()

    FEATURECOUNTS(
        genome1_bams,
        DOWNLOAD_REFERENCES.out.gtf,
        "genome1_${params.strain1}"
    )

    FEATURECOUNTS(
        genome2_bams,
        DOWNLOAD_REFERENCES.out.gtf,
        "genome2_${params.strain2}"
    )

    FEATURECOUNTS(
        ref_bams,
        DOWNLOAD_REFERENCES.out.gtf,
        "reference"
    )

    // ── 8. MultiQC ─────────────────────────────────────────────────────
    mqc_files = Channel.empty()
        .mix(
            bam_nmask.log.collect(),
            bam_ref.log.collect(),
            FEATURECOUNTS.out.summary.collect()
        )
        .collect()

    MULTIQC(mqc_files)
}

// ── Completion message ──────────────────────────────────────────────────────

workflow.onComplete {
    log.info """
══════════════════════════════════════════════════════════════════
  Duration : ${workflow.duration}
  Output   : ${params.outdir}
  Counts   : ${params.outdir}/07_counts/
══════════════════════════════════════════════════════════════════
""".stripIndent()
}
