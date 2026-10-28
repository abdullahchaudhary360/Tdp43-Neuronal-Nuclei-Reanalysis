#!/usr/bin/env bash
# Reference FASTQ-to-counts workflow for paired-end RNA-seq.
# IMPORTANT: this configurable script was not executed for the present reanalysis.
# The manuscript's statistical results start from the supplied counts_matrix.txt.
# Verify library strandedness and record all tool versions before a new run.

set -euo pipefail

THREADS="${THREADS:-8}"
STRAND="${STRAND:-0}"  # featureCounts: 0=unstranded, 1=stranded, 2=reverse stranded
OUTDIR="${OUTDIR:-$PWD/rnaseq_run}"
SAMPLES_FILE="${SAMPLES_FILE:-samples.txt}"  # one SRR accession per line
GENOME_FASTA="${GENOME_FASTA:-}"
GTF="${GTF:-}"
GENOME_DIR="${GENOME_DIR:-$OUTDIR/star_index}"
ADAPTERS="${ADAPTERS:-}"  # Trimmomatic adapter FASTA

require_var() {
  local name="$1"
  if [[ -z "${!name}" ]]; then
    echo "Set $name before running." >&2
    exit 2
  fi
}
require_var GENOME_FASTA
require_var GTF
require_var ADAPTERS
[[ -s "$SAMPLES_FILE" ]] || { echo "Missing $SAMPLES_FILE" >&2; exit 2; }

for cmd in fasterq-dump fastqc trimmomatic STAR featureCounts samtools; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "Missing command: $cmd" >&2; exit 2; }
done

mkdir -p "$OUTDIR"/{fastq,qc_raw,qc_trimmed,trimmed,aligned,counts,logs}

if [[ ! -s "$GENOME_DIR/genomeParameters.txt" ]]; then
  STAR --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang 119
fi

while read -r sample; do
  [[ -z "$sample" || "$sample" =~ ^# ]] && continue
  echo "Processing $sample"

  if [[ ! -s "$OUTDIR/fastq/${sample}_1.fastq" && ! -s "$OUTDIR/fastq/${sample}_1.fastq.gz" ]]; then
    fasterq-dump "$sample" --split-files --threads "$THREADS" \
      --outdir "$OUTDIR/fastq" >"$OUTDIR/logs/${sample}.download.log" 2>&1
  fi
  [[ -s "$OUTDIR/fastq/${sample}_1.fastq" ]] && gzip -f "$OUTDIR/fastq/${sample}_1.fastq"
  [[ -s "$OUTDIR/fastq/${sample}_2.fastq" ]] && gzip -f "$OUTDIR/fastq/${sample}_2.fastq"

  fastqc -t "$THREADS" -o "$OUTDIR/qc_raw" \
    "$OUTDIR/fastq/${sample}_1.fastq.gz" "$OUTDIR/fastq/${sample}_2.fastq.gz"

  trimmomatic PE -threads "$THREADS" -phred33 \
    "$OUTDIR/fastq/${sample}_1.fastq.gz" "$OUTDIR/fastq/${sample}_2.fastq.gz" \
    "$OUTDIR/trimmed/${sample}_1P.fastq.gz" "$OUTDIR/trimmed/${sample}_1U.fastq.gz" \
    "$OUTDIR/trimmed/${sample}_2P.fastq.gz" "$OUTDIR/trimmed/${sample}_2U.fastq.gz" \
    "ILLUMINACLIP:${ADAPTERS}:2:30:10" LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 MINLEN:36 >"$OUTDIR/logs/${sample}.trimmomatic.log" 2>&1

  fastqc -t "$THREADS" -o "$OUTDIR/qc_trimmed" \
    "$OUTDIR/trimmed/${sample}_1P.fastq.gz" "$OUTDIR/trimmed/${sample}_2P.fastq.gz"

  STAR --runThreadN "$THREADS" \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$OUTDIR/trimmed/${sample}_1P.fastq.gz" "$OUTDIR/trimmed/${sample}_2P.fastq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$OUTDIR/aligned/${sample}_" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS NM MD \
    --sjdbGTFfile "$GTF"
  samtools index "$OUTDIR/aligned/${sample}_Aligned.sortedByCoord.out.bam"
done < "$SAMPLES_FILE"

mapfile -t bams < <(find "$OUTDIR/aligned" -name '*_Aligned.sortedByCoord.out.bam' | sort)
featureCounts -T "$THREADS" -p --countReadPairs -s "$STRAND" \
  -a "$GTF" -t exon -g gene_id \
  -o "$OUTDIR/counts/raw_counts.txt" "${bams[@]}"

tail -n +2 "$OUTDIR/counts/raw_counts.txt" | cut -f1,7- > "$OUTDIR/counts/counts_matrix_with_header.txt"
echo "Completed. Record software versions and verify $OUTDIR/counts/raw_counts.txt.summary."
