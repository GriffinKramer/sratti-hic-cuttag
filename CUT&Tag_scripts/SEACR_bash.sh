#!/bin/bash
set -e
CUTTAG_DIR="CUTTag"
GENOME="genome/chrom.sizes"
SEACR_SCRIPT="SEACR_1.3.sh"
SEACR_ENV_BASH="/opt/anaconda3/envs/samtools_env/bin/bash"
OUT_DIR="SEACR_non"

mkdir -p "$OUT_DIR"

for bam in $CUTTAG_DIR/**/*.dedup.bam; do
    sample=$(basename "$bam" .dedup.bam)
    dir=$(dirname "$bam")
    echo "Processing $sample..."

    samtools sort -n -o "$dir/${sample}.namesorted.bam" "$bam"
    samtools view -b -f 0x2 -F 0x4 "$dir/${sample}.namesorted.bam" > "$dir/${sample}.mapped.bam"
    bedtools bamtobed -bedpe -i "$dir/${sample}.mapped.bam" > "$dir/${sample}.bed"
    awk '$1==$4 && $6-$2 < 1000 {print $0}' "$dir/${sample}.bed" > "$dir/${sample}.clean.bed"
    cut -f 1,2,6 "$dir/${sample}.clean.bed" | sort -k1,1 -k2,2n -k3,3n > "$dir/${sample}.fragments.bed"
    bedtools genomecov -bg -i "$dir/${sample}.fragments.bed" -g "$GENOME" > "$dir/${sample}.fragments.bedgraph"
    rm "$dir/${sample}.bed" "$dir/${sample}.clean.bed" "$dir/${sample}.fragments.bed" "$dir/${sample}.mapped.bam"
    echo "Done: $dir/${sample}.fragments.bedgraph"
done
echo "All samples preprocessed."
echo

for bedgraph in $CUTTAG_DIR/**/*.fragments.bedgraph; do
    sample=$(basename "$bedgraph" .fragments.bedgraph)
    dir=$(dirname "$bedgraph")

    if [[ "$sample" == *IgG* ]]; then
        continue
    fi

    matches=($(find "$dir" -maxdepth 1 -iname "*IgG*.fragments.bedgraph"))

    if [[ ${#matches[@]} -eq 0 ]]; then
        echo "WARNING: No IgG control found for $sample in $dir, skipping."
        continue
    elif [[ ${#matches[@]} -gt 1 ]]; then
        echo "WARNING: Multiple IgG controls found for $sample in $dir, using ${matches[0]}."
    fi
    control="${matches[0]}"

    echo "Calling SEACR peaks for $sample vs $(basename "$control")..."
    "$SEACR_ENV_BASH" "$SEACR_SCRIPT" \
        "$bedgraph" \
        "$control" \
        non stringent \
        "$OUT_DIR/${sample}__non_stringent"
done
echo "All SEACR peak calls complete."
