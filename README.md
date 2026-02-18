# SPLIT 🧬

**S**NP-Level **I**nspection of **P**arental **T**ranscripts

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04-23aa62.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Nextflow DSL2 pipeline for **allele-specific RNA-seq** analysis using STAR, SNPsplit, and featureCounts.
Ultra-minimalist — 2 files only (`main.nf` + `nextflow.config`). Designed for solo bioinformaticians.

---

## Pipeline Overview

```mermaid
flowchart TD
    subgraph INPUT ["Input"]
        direction TB
        SRA["SRR / ERR / DRR"]
        GEO["GSE / GSM"]
        FQ_DIR["FASTQ directory"]
        CSV["CSV samplesheet"]
    end

    GEO -->|"NCBI E-utils"| RESOLVE_GEO["RESOLVE_GEO"]
    SRA --> SRA_DOWNLOAD["SRA_DOWNLOAD"]
    RESOLVE_GEO -->|"SRR IDs"| SRA_DOWNLOAD
    FQ_DIR --> FASTQS(("FASTQs"))
    CSV --> FASTQS
    SRA_DOWNLOAD --> FASTQS

    subgraph REFS ["References"]
        direction TB
        DOWNLOAD["DOWNLOAD_REFERENCES"]
    end

    DOWNLOAD -->|"FASTA + VCF"| SNPSPLIT_PREP["SNPSPLIT_GENOME_PREP"]

    FASTQS -->|"first FASTQ"| READ_LEN["ESTIMATE_READ_LENGTH"]

    SNPSPLIT_PREP -->|"N-masked FASTA"| IDX_NMASK["STAR_INDEX\n(N-masked)"]
    DOWNLOAD -->|"ref FASTA"| IDX_REF["STAR_INDEX\n(reference)"]
    DOWNLOAD -->|"GTF"| IDX_NMASK
    DOWNLOAD -->|"GTF"| IDX_REF
    READ_LEN -->|"sjdbOverhang"| IDX_NMASK
    READ_LEN -->|"sjdbOverhang"| IDX_REF

    subgraph NMASK_TRACK ["N-masked track"]
        direction TB
        ALIGN_NM["STAR_ALIGN\n(N-masked)"]
        SORT_NM["SORT_DEDUP"]
        SPLIT_ALLELE["SNPSPLIT"]
    end

    FASTQS --> ALIGN_NM
    IDX_NMASK --> ALIGN_NM
    ALIGN_NM --> SORT_NM
    SORT_NM --> SPLIT_ALLELE
    SNPSPLIT_PREP -->|"SNP file"| SPLIT_ALLELE

    subgraph REF_TRACK ["Reference track"]
        direction TB
        ALIGN_REF["STAR_ALIGN\n(reference)"]
        SORT_REF["SORT_DEDUP"]
    end

    FASTQS --> ALIGN_REF
    IDX_REF --> ALIGN_REF
    ALIGN_REF --> SORT_REF

    SPLIT_ALLELE -->|"genome1 BAMs"| FC_G1["FEATURECOUNTS\n(genome1)"]
    SPLIT_ALLELE -->|"genome2 BAMs"| FC_G2["FEATURECOUNTS\n(genome2)"]
    SORT_REF --> FC_REF["FEATURECOUNTS\n(reference)"]

    DOWNLOAD -->|"GTF"| FC_G1
    DOWNLOAD -->|"GTF"| FC_G2
    DOWNLOAD -->|"GTF"| FC_REF

    ALIGN_NM -->|"STAR logs"| MULTIQC["MULTIQC"]
    ALIGN_REF -->|"STAR logs"| MULTIQC
    FC_G1 -->|"summary"| MULTIQC
    FC_G2 -->|"summary"| MULTIQC
    FC_REF -->|"summary"| MULTIQC

    FC_G1 --> OUT_G1["genome1 counts\n(strain1)"]
    FC_G2 --> OUT_G2["genome2 counts\n(strain2)"]
    FC_REF --> OUT_REF["reference counts"]
    MULTIQC --> OUT_QC["MultiQC report"]

    classDef inputStyle fill:#3498db,stroke:#2c6fbb,color:#fff,stroke-width:2px
    classDef processStyle fill:#2c3e50,stroke:#1a252f,color:#fff,stroke-width:2px
    classDef keyProcess fill:#c0392b,stroke:#922b21,color:#fff,stroke-width:3px
    classDef outputStyle fill:#27ae60,stroke:#1e8449,color:#fff,stroke-width:2px
    classDef dataNode fill:#f39c12,stroke:#d68910,color:#fff,stroke-width:2px

    class SRA,GEO,FQ_DIR,CSV inputStyle
    class RESOLVE_GEO,SRA_DOWNLOAD,DOWNLOAD,SNPSPLIT_PREP,IDX_NMASK,IDX_REF,READ_LEN,ALIGN_NM,SORT_NM,ALIGN_REF,SORT_REF,FC_G1,FC_G2,FC_REF,MULTIQC processStyle
    class SPLIT_ALLELE keyProcess
    class OUT_G1,OUT_G2,OUT_REF,OUT_QC outputStyle
    class FASTQS dataNode
```

## Quick Start

```bash
# From a FASTQ directory (auto-detects PE/SE)
nextflow run IPNP-BIPN/SPLIT --fastq_dir /path/to/fastqs --outdir results -resume

# From SRA accessions
nextflow run IPNP-BIPN/SPLIT --sra_ids "SRR1234567,SRR1234568" --outdir results -resume

# From a GEO dataset (auto-resolves GSE → SRR)
nextflow run IPNP-BIPN/SPLIT --sra_ids GSE80810 --outdir results -resume

# From a samplesheet CSV
nextflow run IPNP-BIPN/SPLIT --input samplesheet.csv --outdir results -resume

# Custom strains
nextflow run IPNP-BIPN/SPLIT \
    --fastq_dir /path/to/fastqs \
    --strain1 CAST_EiJ \
    --strain2 C57BL_6NJ \
    --dedup true \
    --outdir results \
    -resume
```

### Samplesheet format (CSV)

```csv
sample,fastq_1,fastq_2
sampleA,/path/to/sampleA_R1_001.fastq.gz,/path/to/sampleA_R2_001.fastq.gz
sampleB,/path/to/sampleB.fastq.gz,
```

> Leave `fastq_2` empty for single-end reads.

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | `null` | Samplesheet CSV (`sample,fastq_1,fastq_2`) |
| `--fastq_dir` | `null` | Directory of FASTQs (`*.fastq.gz`) |
| `--sra_ids` | `null` | SRA/GEO accessions (comma-separated or file) |
| `--outdir` | `results` | Output directory |
| `--strain1` | `CAST_EiJ` | First strain in VCF → genome1 |
| `--strain2` | `C57BL_6NJ` | Second strain in VCF → genome2 |
| `--dedup` | `true` | Remove PCR duplicates (samtools markdup -r) |
| `--strandedness` | `0` | featureCounts strandedness (0/1/2) |
| `--force_se` | `false` | Force single-end counting in featureCounts |
| `--genome_url` | Ensembl GRCm39 | Genome FASTA URL |
| `--gtf_url` | Ensembl 2023_04 | GTF annotation URL |
| `--vcf_url` | MGP REL2021 v8 | VCF SNPs URL |
| `--star_limit_genome_ram` | `60000000000` | STAR --limitGenomeGenerateRAM |
| `--max_cpus` | auto | Maximum number of CPUs |
| `--max_memory` | `64 GB` | Maximum memory (scales all process labels) |

### Pre-built references (skip downloads)

| Parameter | Description |
|-----------|-------------|
| `--genome_fa` | Pre-downloaded genome FASTA |
| `--gtf` | Pre-downloaded GTF |
| `--vcf` | Pre-downloaded VCF |
| `--star_index_nmask` | Pre-built STAR index (N-masked) |
| `--star_index_ref` | Pre-built STAR index (reference) |
| `--snp_file` | Pre-existing SNPsplit SNP annotation file |

---

## Output Structure

```
results/
├── 00_sra_fastq/           # Downloaded FASTQs (if SRA input)
├── 04_aln_nmask/           # N-mask aligned BAMs (sorted, deduped)
├── 05_aln_ref/             # Reference aligned BAMs (sorted, deduped)
├── 06_snpsplit/            # SNPsplit output (genome1, genome2, unassigned)
├── 07_counts/
│   ├── counts_genome1_CAST_EiJ.txt      # Allele-specific counts (strain1)
│   ├── counts_genome2_C57BL_6NJ.txt     # Allele-specific counts (strain2)
│   └── counts_reference.txt              # Standard reference counts
├── 08_multiqc/             # Aggregated QC report
├── reference/              # Downloaded + cached references
│   ├── genome.fa           # GRCm39 soft-masked
│   ├── genes.gtf           # Ensembl annotation
│   ├── snps.vcf.gz         # MGP REL2021 SNPs
│   ├── snpsplit_prep/      # N-masked genome + SNP file
│   ├── star_nmask/         # STAR index (N-masked)
│   └── star_ref/           # STAR index (reference)
└── pipeline_info/          # Nextflow timeline, trace, DAG, report
```

---

## Requirements

**Core** (always required):
`STAR` `samtools` `SNPsplit` `featureCounts` (subread) `multiqc` `wget`

**Optional**:
`sra-tools` `bgzip` (htslib/tabix) — for SRA download

**Nextflow** ≥ 23.04

---

## How it works

1. **References** are automatically downloaded from Ensembl (GRCm39 genome + GTF) and MGP (SNP VCF). Cached via `storeDir` — only downloaded once.

2. **SNPsplit genome preparation** creates an N-masked genome where strain-discriminating SNP positions are replaced by N. This prevents alignment bias toward the reference allele.

3. **Two parallel alignment tracks**:
   - **N-masked track**: alignments used for allele-specific analysis (SNPsplit)
   - **Reference track**: standard alignments for total gene expression

4. **SNPsplit** assigns each read from the N-mask track to genome1 (strain1), genome2 (strain2), or unassigned based on informative SNP positions.

5. **featureCounts** produces three count tables using `gene_name` attribute for human-readable gene symbols (e.g., *Gapdh* instead of ENSMUSG00000057666).

---

## Resume & Cache

The pipeline natively leverages Nextflow's cache (`-resume`). Already completed steps are automatically skipped. References (genome, GTF, VCF, STAR indexes, N-masked genome) are persisted via `storeDir` and reused across runs.

```bash
# Re-run after a crash — picks up exactly where it left off
nextflow run main.nf --fastq_dir fastqs --outdir results -resume
```

---

## License

MIT
