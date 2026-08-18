import re
from pathlib import Path

CUTTAG_DIR = Path("CUTTag")
IN_DIR = Path("SEACR")
OUT_DIR = Path("SEACR_filtered")
FE_THRESH = 2.0
SUFFIX_PATTERN = re.compile(r"(_stringent\.stringent\.bed|\.stringent\.bed|\.relaxed\.bed)$")

def find_bedgraph(sample):
    matches = list(CUTTAG_DIR.glob(f"*/{sample}.fragments.bedgraph"))
    return matches[0] if matches else None


def genome_mean_depth(bedgraph_path):
    total_signal = 0.0
    total_width = 0.0
    with open(bedgraph_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            start, end, depth = int(parts[1]), int(parts[2]), float(parts[3])
            width = end - start
            total_signal += width * depth
            total_width += width
    return total_signal / total_width if total_width > 0 else 0.0


def main():
    OUT_DIR.mkdir(exist_ok=True)

    for peak_file in sorted(IN_DIR.glob("*.bed")):
        sample = SUFFIX_PATTERN.sub("", peak_file.name)

        bedgraph = find_bedgraph(sample)
        if bedgraph is None:
            print(f"WARNING: no matching bedgraph for {peak_file.name} "
                  f"(looked for {sample}.fragments.bedgraph), skipping.")
            continue

        gm = genome_mean_depth(bedgraph)

        total = 0
        kept = 0
        outfile = OUT_DIR / peak_file.name.replace("_", "_filtered_", 1)
        fe_outfile = outfile.with_name(outfile.name + ".fe_values.tsv")
        with open(peak_file) as fin, open(outfile, "w") as fout, open(fe_outfile, "w") as ffe:
            ffe.write("chrom\tstart\tend\tfold_enrichment\n")
            for line in fin:
                total += 1
                parts = line.rstrip("\n").split("\t")
                chrom, start, end, signal = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
                width = end - start
                avg_depth = signal / width if width > 0 else 0.0
                fe = avg_depth / gm if gm > 0 else 0.0
                if fe >= FE_THRESH:
                    fout.write(line)
                    ffe.write(f"{chrom}\t{start}\t{end}\t{fe:.4f}\n")
                    kept += 1

        print(f"{peak_file.name}: total={total}  kept={kept} (FE>={FE_THRESH}, "
              f"genome_mean={gm:.4f})  -> {outfile}  (FE values: {fe_outfile})")

    print(f"Done. Filtered peak files written to {OUT_DIR}/")


if __name__ == "__main__":
    main()
