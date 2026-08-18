import glob
import os
from collections import defaultdict
import pandas as pd

# config
ISLANDS_PATH = "genome/sratti_parasitism_islands.bed"
WINDOWS = {"exact": 0, "500bp": 500, "1000bp": 1000}
Q_THRESHOLD = 0.01

CONTACT_SOURCES = {
    "ArimaFLF": {
        "path": "fithic/ArimaFLF/ArimaFLF_10kb.spline_pass2.res10000.significances.txt.gz",
        "half_window": 5000,
    },
    "ArimaPF": {
        "path": "fithic/ArimaPF/ArimaPF_10kb.spline_pass2.res10000.significances.txt.gz",
        "half_window": 5000,
    },
    "DovetailFLF": {
        "path": "fithic/DovetailFLF/DovetailFLF_10kb.spline_pass2.res10000.significances.txt.gz",
        "half_window": 5000,
    },
    "DovetailPF": {
        "path": "fithic/DovetailPF/DovetailPF_10kb.spline_pass2.res10000.significances.txt.gz",
        "half_window": 5000,
    },
}

LOOP_SOURCES = {
    "DovetailFLF": "cooltools/DovetailFLF_dots.bedpe",
    "DovetailPF": "cooltools/DovetailPF_dots.bedpe",
}

SEACR_DIR = "SEACR_filtered"
PEAK_SAMPLES = {
    ("FL", "FL1", "K4"): "FLK4",
    ("FL", "FL1", "K9"): "FLK9",
    ("FL", "FL1", "K27"): "FLK27",
    ("FL", "FL1", "K36"): "FLK36",
    ("FL", "FL1", "K92"): "FLK92",
    ("FL", "FL2", "K4"): "FL2K4",
    ("FL", "FL2", "K9"): "FL2K9",
    ("FL", "FL2", "K27"): "FL2K27",
    ("FL", "FL2", "K36"): "FL2K36",
    ("FL", "FL2", "K92"): "FL2K92",
    ("PF", None, "K4"): "PFK4",
    ("PF", None, "K9"): "PFK9",
    ("PF", None, "K27"): "PFK27",
    ("PF", None, "K36"): "PFK36",
    ("PF", None, "K92"): "PFK92",
}

OUT_PATH = "overlap_analysis_filtered.xlsx"

# load files
def find_seacr_file(prefix):
    """SEACR filenames carry long EKDL sequencer codes, so glob by prefix."""
    matches = glob.glob(f"{SEACR_DIR}/{prefix}_*_stringent.stringent.bed")
    if len(matches) == 0:
        raise FileNotFoundError(
            f"No SEACR file found for prefix '{prefix}' in {SEACR_DIR}/"
        )
    if len(matches) > 1:
        raise ValueError(f"Multiple SEACR files match prefix '{prefix}': {matches}")
    return matches[0]


def load_peaks(prefix):
    path = find_seacr_file(prefix)
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "auc", "max", "region"],
    )


def load_islands():
    islands = pd.read_csv(
        ISLANDS_PATH,
        sep=r"\s+",
        header=None,
        names=["chr", "start", "end", "genes"],
        dtype={"start": int, "end": int},
    )
    print(f"Loaded {len(islands)} parasitism islands")
    return islands


def extend(df, window):
    ext = df.copy()
    ext["start"] = (ext["start"] - window).clip(lower=0)
    ext["end"] = ext["end"] + window
    return ext


def load_contacts(source_name):
    cfg = CONTACT_SOURCES[source_name]
    df = pd.read_csv(cfg["path"], sep="\t")
    sig = df[df["q-value"] < Q_THRESHOLD].copy()
    hw = cfg["half_window"]
    sig["start1"] = sig["fragmentMid1"] - hw
    sig["end1"] = sig["fragmentMid1"] + hw
    sig["start2"] = sig["fragmentMid2"] - hw
    sig["end2"] = sig["fragmentMid2"] + hw
    print(f"{source_name}: {len(sig)} significant contacts (q<{Q_THRESHOLD})")
    return sig


def load_loops(source_name):
    df = pd.read_csv(LOOP_SOURCES[source_name], sep="\t")
    df = df.rename(columns={"chrom1": "chr1", "chrom2": "chr2"})
    print(f"{source_name}: {len(df)} loops")
    return df

def build_single_index(df, chr_col="chr", start_col="start", end_col="end"):
    idx = defaultdict(list)
    for c, s, e in zip(df[chr_col], df[start_col], df[end_col]):
        idx[c].append((s, e))
    return idx


def any_overlap(chrom, start, end, idx):
    for s, e in idx.get(chrom, ()):
        if s < end and e > start:
            return True
    return False


def build_pair_index(
    df, chr1="chr1", s1="start1", e1="end1", chr2="chr2", s2="start2", e2="end2"
):
    idx = defaultdict(list)
    for c1, a, b, c, d in zip(df[chr1], df[s1], df[e1], df[s2], df[e2]):
        idx[c1].append((a, b, c, d))
    return idx

# overlap analysis
def twoanchor_vs_single(twoanchor_df, single_idx):
    a1 = [
        any_overlap(c, s, e, single_idx)
        for c, s, e in zip(
            twoanchor_df["chr1"], twoanchor_df["start1"], twoanchor_df["end1"]
        )
    ]
    a2 = [
        any_overlap(c, s, e, single_idx)
        for c, s, e in zip(
            twoanchor_df["chr2"], twoanchor_df["start2"], twoanchor_df["end2"]
        )
    ]
    total = len(twoanchor_df)
    both = sum(x and y for x, y in zip(a1, a2))
    either = sum(x or y for x, y in zip(a1, a2))
    return both, either, total - either, total


def single_vs_single(single_df, ref_idx):
    n = sum(
        any_overlap(c, s, e, ref_idx)
        for c, s, e in zip(single_df["chr"], single_df["start"], single_df["end"])
    )
    total = len(single_df)
    return n, total - n, total


def twoanchor_vs_twoanchor(a_df, b_idx):
    matches = 0
    for c, s1, e1, s2, e2 in zip(
        a_df["chr1"], a_df["start1"], a_df["end1"], a_df["start2"], a_df["end2"]
    ):
        for bs1, be1, bs2, be2 in b_idx.get(c, ()):
            if (s1 < be1 and e1 > bs1 and s2 < be2 and e2 > bs2) or (
                s1 < be2 and e1 > bs2 and s2 < be1 and e2 > bs1
            ):
                matches += 1
                break
    return matches, len(a_df)


def pct(n, total):
    return round(n / total * 100, 1) if total else float("nan")


# main
os.makedirs("pi_overlap", exist_ok=True)

islands = load_islands()
islands_windows = {w: extend(islands, size) for w, size in WINDOWS.items()}
islands_idx = {w: build_single_index(df) for w, df in islands_windows.items()}

peaks = {key: load_peaks(prefix) for key, prefix in PEAK_SAMPLES.items()}
peaks_idx = {key: build_single_index(df) for key, df in peaks.items()}

contacts = {name: load_contacts(name) for name in CONTACT_SOURCES}
loops = {name: load_loops(name) for name in LOOP_SOURCES}
loops_pair_idx = {name: build_pair_index(df) for name, df in loops.items()}

summary_rows = []
sheets = {}


def add_rows(comparison, rows):
    for r in rows:
        r["comparison"] = comparison
    summary_rows.extend(rows)
    sheets[comparison] = pd.DataFrame(rows)


for stage, contact_key, sheet_name in [
    ("FL", "ArimaFLF", "ArimaFLF_Contacts_vs_Peaks"),
    ("PF", "ArimaPF", "ArimaPF_Contacts_vs_Peaks"),
]:
    rows = []
    for (pstage, rep, mark), idx in peaks_idx.items():
        if pstage != stage:
            continue
        both, either, neither, total = twoanchor_vs_single(contacts[contact_key], idx)
        rows.append(
            {
                "replicate": rep,
                "mark": mark,
                "both": both,
                "either": either,
                "neither": neither,
                "total": total,
                "pct_both": pct(both, total),
                "pct_either": pct(either, total),
            }
        )
    add_rows(sheet_name, rows)

for contact_key, sheet_name in [
    ("ArimaFLF", "ArimaFLF_Contacts_vs_PI"),
    ("ArimaPF", "ArimaPF_Contacts_vs_PI"),
]:
    rows = []
    for w, idx in islands_idx.items():
        both, either, neither, total = twoanchor_vs_single(contacts[contact_key], idx)
        rows.append(
            {
                "window": w,
                "both": both,
                "either": either,
                "neither": neither,
                "total": total,
                "pct_both": pct(both, total),
                "pct_either": pct(either, total),
            }
        )
    add_rows(sheet_name, rows)

for stage, contact_key, loop_key, sheet_name in [
    ("FL", "DovetailFLF", "DovetailFLF", "DoveFLF_Contacts_vs_Loops"),
    ("PF", "DovetailPF", "DovetailPF", "DovePF_Contacts_vs_Loops"),
]:
    c_matches, c_total = twoanchor_vs_twoanchor(
        contacts[contact_key], loops_pair_idx[loop_key]
    )
    l_matches, l_total = twoanchor_vs_twoanchor(
        loops[loop_key], build_pair_index(contacts[contact_key])
    )
    rows = [
        {
            "contacts_with_loop_match": c_matches,
            "total_contacts": c_total,
            "pct_contacts_matched": pct(c_matches, c_total),
            "loops_with_contact_match": l_matches,
            "total_loops": l_total,
            "pct_loops_matched": pct(l_matches, l_total),
        }
    ]
    add_rows(sheet_name, rows)

for stage, contact_key, sheet_name in [
    ("FL", "DovetailFLF", "DoveFLF_Contacts_vs_Peaks"),
    ("PF", "DovetailPF", "DovePF_Contacts_vs_Peaks"),
]:
    rows = []
    for (pstage, rep, mark), idx in peaks_idx.items():
        if pstage != stage:
            continue
        both, either, neither, total = twoanchor_vs_single(contacts[contact_key], idx)
        rows.append(
            {
                "replicate": rep,
                "mark": mark,
                "both": both,
                "either": either,
                "neither": neither,
                "total": total,
                "pct_both": pct(both, total),
                "pct_either": pct(either, total),
            }
        )
    add_rows(sheet_name, rows)

for contact_key, sheet_name in [
    ("DovetailFLF", "DoveFLF_Contacts_vs_PI"),
    ("DovetailPF", "DovePF_Contacts_vs_PI"),
]:
    rows = []
    for w, idx in islands_idx.items():
        both, either, neither, total = twoanchor_vs_single(contacts[contact_key], idx)
        rows.append(
            {
                "window": w,
                "both": both,
                "either": either,
                "neither": neither,
                "total": total,
                "pct_both": pct(both, total),
                "pct_either": pct(either, total),
            }
        )
    add_rows(sheet_name, rows)

for stage, loop_key, sheet_name in [
    ("FL", "DovetailFLF", "DoveFLF_Loops_vs_Peaks"),
    ("PF", "DovetailPF", "DovePF_Loops_vs_Peaks"),
]:
    rows = []
    for (pstage, rep, mark), idx in peaks_idx.items():
        if pstage != stage:
            continue
        both, either, neither, total = twoanchor_vs_single(loops[loop_key], idx)
        rows.append(
            {
                "replicate": rep,
                "mark": mark,
                "both": both,
                "either": either,
                "neither": neither,
                "total": total,
                "pct_both": pct(both, total),
                "pct_either": pct(either, total),
            }
        )
    add_rows(sheet_name, rows)

for loop_key, sheet_name in [
    ("DovetailFLF", "DoveFLF_Loops_vs_PI"),
    ("DovetailPF", "DovePF_Loops_vs_PI"),
]:
    rows = []
    for w, idx in islands_idx.items():
        both, either, neither, total = twoanchor_vs_single(loops[loop_key], idx)
        rows.append(
            {
                "window": w,
                "both": both,
                "either": either,
                "neither": neither,
                "total": total,
                "pct_both": pct(both, total),
                "pct_either": pct(either, total),
            }
        )
    add_rows(sheet_name, rows)

for stage, sheet_name in [("FL", "FLF_Peaks_vs_PI"), ("PF", "PF_Peaks_vs_PI")]:
    rows = []
    for (pstage, rep, mark), pdf in peaks.items():
        if pstage != stage:
            continue
        for w, idx in islands_idx.items():
            n, not_n, total = single_vs_single(pdf, idx)
            rows.append(
                {
                    "replicate": rep,
                    "mark": mark,
                    "window": w,
                    "overlapping": n,
                    "not_overlapping": not_n,
                    "total": total,
                    "pct_overlapping": pct(n, total),
                }
            )
    add_rows(sheet_name, rows)

summary_df = pd.DataFrame(summary_rows)
cols = ["comparison"] + [c for c in summary_df.columns if c != "comparison"]
summary_df = summary_df[cols]

with pd.ExcelWriter(OUT_PATH, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="Summary", index=False)
    for sheet_name, df in sheets.items():
        df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

print(f"\nSaved {OUT_PATH} with {len(sheets) + 1} sheets.")
