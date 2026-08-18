#!/bin/bash
set -e

CUTTAG_DIR="CUTTag"
OUT_DIR="MACS3"
GENOME_SIZE=44000000

mkdir -p "$OUT_DIR"

BROAD_MARKS=("K27" "K36" "K92" "K9")

for bam in $CUTTAG_DIR/**/*.dedup.bam; do
    sample=$(basename "$bam" .dedup.bam)
    dir=$(dirname "$bam")

    if [[ "$sample" == *IgG* ]]; then
        continue
    fi

    matches=($(find "$dir" -maxdepth 1 -iname "*IgG*.dedup.bam"))
    if [[ ${#matches[@]} -eq 0 ]]; then
        echo "WARNING: No IgG control found for $sample in $dir, skipping."
        continue
    elif [[ ${#matches[@]} -gt 1 ]]; then
        echo "WARNING: Multiple IgG controls found for $sample in $dir, using ${matches[0]}."
    fi
    control="${matches[0]}"

    broad_flag=""
    for mark in "${BROAD_MARKS[@]}"; do
        if [[ "$sample" == *"$mark" ]]; then
            broad_flag="--broad"
            break
        fi
    done

    echo "Calling MACS3 peaks for $sample vs $(basename "$control") (mode: ${broad_flag:-narrow})..."
    macs3 callpeak \
        -t "$bam" \
        -c "$control" \
        -f BAMPE --nomodel \
        -g "$GENOME_SIZE" \
        --slocal 500 --llocal 5000 \
        $broad_flag \
        -n "$sample" --outdir "$OUT_DIR"
done

echo "All MACS3 peak calls complete."
