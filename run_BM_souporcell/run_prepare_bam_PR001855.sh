#!/usr/bin/env bash
#SBATCH -J sort_bam_PR001855  # Job name
#SBATCH -o sort_%j.out         # Standard output (%j = job ID)
#SBATCH -e sort_%j.err         # Standard error
#SBATCH -c 16                  # Number of CPU cores
#SBATCH --mem=150G             # Memory limit
#SBATCH -t 1:00:00            # Time limit (hh:mm:ss)
#SBATCH --mail-type=END,FAIL   # (optional) email notifications
#SBATCH --mail-user=wenjie.sun@vai.org
set -euo pipefail

threads=10

in_dir="../data/Janssens_Lab_Data/BMMC_Beds/PR001855_sciCoCnT_JAND_DNA_ALIGNMENT_WELLID_BARCODES/BED"
out_dir="./tmp/PR001855_bam/"
mkdir -p "$out_dir"

## Remove duplicates
for in_bam in "$in_dir"/*.bam; do
    echo "Processing $in_bam"

    in_bam_name=$(basename "$in_bam")
    in_bam_name_prefix=${in_bam_name%.bam}

    samtools sort -n -m 5G -@ $threads "$in_bam" \
      | samtools fixmate -m -@ $threads - - \
      | samtools sort -m 5G -@ $threads - \
      | samtools markdup -r -s \
          -f "${out_dir}/${in_bam_name}.dup.out" \
          --barcode-tag CB \
          -@ $threads \
          - "${out_dir}/${in_bam_name_prefix}.dedup.bam"
done

## Update the read groups
for dedup_bam in "$out_dir"/*.dedup.bam; do
    echo "Updating read groups for $dedup_bam"
    python ../../lib/bin/bin_add_tags.py "$dedup_bam" | samtools view -b -@ $threads -o "${dedup_bam%.bam}.rg.bam" -
done

## Merge all BAM files into one
samtools merge -@ $threads -f "${out_dir}/merged.bam" "$out_dir"/*.dedup.rg.bam

## Sort and index the merged BAM file
samtools sort -m 10G -@ $threads -o "${out_dir}/merged.sorted.bam" "${out_dir}/merged.bam"
samtools index -@ $threads "${out_dir}/merged.sorted.bam"

