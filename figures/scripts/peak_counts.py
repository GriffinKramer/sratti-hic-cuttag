
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

MACS3_DIR = Path("MACS3")
SEACR_DIR = Path("SEACR")
SEACR_NON_DIR = Path("SEACR_non")
SEACR_FILTERED_DIR = Path("SEACR_filtered")
OUT_DIR = Path("peak_compare")
OUT_DIR.mkdir(exist_ok=True)

def count_lines(path):
    with open(path) as f:
        return sum(1 for _ in f)

def clean_sample(filename):
    return filename.split("_")[0]

CONDITIONS = ["FL2", "FL", "PF"]
MARK_ORDER = ["K4", "K27", "K36", "K92", "K9"]


def split_condition_mark(sample):
    for cond in CONDITIONS:
        if sample.startswith(cond):
            return cond, sample[len(cond):]
    return "", sample


def sample_sort_key(sample):
    condition, mark = split_condition_mark(sample)
    mark_rank = MARK_ORDER.index(mark) if mark in MARK_ORDER else len(MARK_ORDER)
    return (mark_rank, condition)


def collect_macs3():
    counts = {}
    files = sorted(MACS3_DIR.glob("*_peaks.narrowPeak")) + sorted(MACS3_DIR.glob("*_peaks.broadPeak"))
    for path in files:
        counts[clean_sample(path.name)] = count_lines(path)
    return counts


def collect_seacr_non():
    counts = {}
    for path in sorted(SEACR_NON_DIR.glob("*__non_stringent.stringent.bed")):
        counts[clean_sample(path.name)] = count_lines(path)
    return counts


def collect_seacr_norm():
    counts = {}
    for path in sorted(SEACR_DIR.glob("*_stringent.stringent.bed")):
        counts[clean_sample(path.name)] = count_lines(path)
    return counts


def collect_seacr_filtered():
    counts = {}
    for path in sorted(SEACR_FILTERED_DIR.glob("*_filtered_*stringent.stringent.bed")):
        counts[clean_sample(path.name)] = count_lines(path)
    return counts


def make_plot(df, out_path):
    ax = df.plot(kind="bar", figsize=(max(12, 0.6 * len(df.index)), 6), width=0.8)
    ax.set_ylabel("Peak count")
    ax.set_xlabel("Sample")
    ax.set_title("Peak counts by peak-calling method")
    ax.legend(title="Method", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved figure: {out_path}")


def main():
    methods = {
        "MACS3": collect_macs3(),
        "SEACR (non)": collect_seacr_non(),
        "SEACR (norm)": collect_seacr_norm(),
        "SEACR (norm, FE>=2)": collect_seacr_filtered(),
    }

    for name, counts in methods.items():
        print(f"{name}: {len(counts)} samples found")

    all_samples = sorted(set().union(*[set(c) for c in methods.values()]), key=sample_sort_key)
    if not all_samples:
        raise SystemExit(
            "No peak files found in MACS3/, SEACR/, or SEACR_filtered/ -- "
            "check that the paths above match where you actually ran the scripts."
        )

    df = pd.DataFrame(
        {name: [counts.get(s, 0) for s in all_samples] for name, counts in methods.items()},
        index=all_samples,
    )
    df.index.name = "sample"
    df.to_csv(OUT_DIR / "peak_counts_by_method.csv")
    print(f"\nSaved counts table: {OUT_DIR / 'peak_counts_by_method.csv'}")

    make_plot(df, out_path=OUT_DIR / "peak_counts_by_method.png")


if __name__ == "__main__":
    main()
