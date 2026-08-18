import os

import cooler
import pandas as pd

RUNS = [
    {
        "sample": "DovetailFLF",
        "resolution": 5000,
        "mcool": "HiC/Dovetail_FLF_Merged/Dovetail_FLF_Merged_balanced.mcool",
    },
    {
        "sample": "DovetailPF",
        "resolution": 5000,
        "mcool": "HiC/Dovetail_PF_Merged/Dovetail_PF_Merged_balanced.mcool",
    },
    {
        "sample": "DovetailFLF",
        "resolution": 10000,
        "mcool": "HiC/Dovetail_FLF_Merged/Dovetail_FLF_Merged_balanced.mcool",
    },
    {
        "sample": "DovetailPF",
        "resolution": 10000,
        "mcool": "HiC/Dovetail_PF_Merged/Dovetail_PF_Merged_balanced.mcool",
    },
    {
        "sample": "ArimaFLF",
        "resolution": 10000,
        "mcool": "HiC/Arima_FLF/Arima_FLF_unbalanced.mcool",
    },
    {
        "sample": "ArimaPF",
        "resolution": 10000,
        "mcool": "HiC/Arima_PF/Arima_PF_unbalanced.mcool",
    },
]
OUT_DIR = "fithic"

def generate_fithic_inputs(sample, resolution, mcool_path, out_dir):
    label = f"{sample}_{resolution // 1000}kb"
    print(f"--- {label} ---")

    clr = cooler.Cooler(f"{mcool_path}::resolutions/{resolution}")
    bins = clr.bins()[:]
    pixels = clr.pixels()[:]

    pixels = pixels[pixels["bin1_id"] != pixels["bin2_id"]].reset_index(drop=True)

    counts1 = pixels.groupby("bin1_id")["count"].sum()
    counts2 = pixels.groupby("bin2_id")["count"].sum()
    bin_counts = counts1.add(counts2, fill_value=0)
    bins["mid"] = (bins["start"] + bins["end"]) // 2
    bins["count"] = bins.index.map(bin_counts).fillna(0).astype(int)
    bins["mappable"] = (bins["count"] > 0).astype(int)
    fragments = pd.DataFrame(
        {
            "chr": bins["chrom"],
            "extra": 0,
            "mid": bins["mid"],
            "count": bins["count"],
            "mappable": bins["mappable"],
        }
    )
    fragments.to_csv(
        f"{out_dir}/{label}_fragments.txt.gz",
        sep="\t",
        header=False,
        index=False,
        compression="gzip",
    )
    print("Fragments done")

    mid1 = bins["mid"].iloc[pixels["bin1_id"]].values
    mid2 = bins["mid"].iloc[pixels["bin2_id"]].values
    chr1 = bins["chrom"].iloc[pixels["bin1_id"]].values
    chr2 = bins["chrom"].iloc[pixels["bin2_id"]].values
    interactions = pd.DataFrame(
        {
            "chr1": chr1,
            "mid1": mid1,
            "chr2": chr2,
            "mid2": mid2,
            "count": pixels["count"].values,
        }
    )
    interactions.to_csv(
        f"{out_dir}/{label}_interactions.txt.gz",
        sep="\t",
        header=False,
        index=False,
        compression="gzip",
    )
    print("Interactions done")


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    for run in RUNS:
        generate_fithic_inputs(run["sample"], run["resolution"], run["mcool"], OUT_DIR)
        print()
