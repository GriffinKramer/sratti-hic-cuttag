# Integrating Hi-C and CUT&Tag to Profile Chromatin Organization and Parasitism Regulation in *Strongyloides ratti*

## Overview

This repository contains the analysis pipeline for a study integrating Hi-C and CUT&Tag to characterize 3D genome organization and histone modification landscapes in the parasitic nematode *Strongyloides ratti*, and their relationship to parasitism gene regulation.

Two Hi-C library preparation kits (**Arima-HiC+** and **Dovetail LinkPrep**) were compared across free-living and parasitic female life stages, and CUT&Tag was used to profile five histone marks (H3K4me3, H3K36me3, H3K9me2, H3K9me3, H3K27me3) across the same conditions. Contact matrices, A/B compartments, and chromatin loops were compared against CUT&Tag peaks and against known parasitism gene islands to test whether parasitism genes occupy a distinct chromatin environment.

Developed as part of the MSc Bioinformatics program at the University of Bath (Hunt Lab).

**Key findings:**
- Activating histone marks (H3K4me3, H3K36me3) show significant Hi-C contact enrichment consistently across kits, life-cycle stages, and replicates
- Repressive marks (H3K9me2, H3K9me3, H3K27me3) show minimal contact overlap
- Hi-C contacts are significantly depleted, rather than enriched, at parasitism gene islands
- Parasitism islands occupy a distinct, repressive chromatin environment relative to the rest of the genome

## Repository Structure
.
├── HiC_scripts/
│   ├── contacts_per_crom.py          # Per-chromosome Hi-C contact counting
│   └── generate_FitHiC_inputs.py     # Builds input files for FitHiC significant contact/loop calling
│
├── CUT&Tag_scripts/
│   ├── SEACR_bash.sh                 # SEACR peak calling
│   ├── filter_SEACR.py               # Post-processing/filtering of SEACR peak calls
│   └── macs3_bash.sh                 # MACS3 peak calling (benchmark comparison)
│
├── overlap_analysis_scripts/
│   ├── circular_genomic_permutation.py  # Monte Carlo circular permutation testing
│   └── overlap_analysis.py              # Overlap of Hi-C contacts/loops, CUT&Tag peaks,
│                                         # and parasitism gene islands
│
├── functional_annotation_scripts/
│   ├── HiC_functional_analysis.py    # PFAM enrichment for genes at significant Hi-C contacts
│   ├── peak_functional_analysis.py   # PFAM enrichment for genes at CUT&Tag peaks
│   └── x2_functional_analysis.py     # PFAM enrichment for genes in A/B compartments
│
├── figures/                          # Figures and the scripts used to generate them
└── README.md

## Methods Summary

**Hi-C processing** — Reads trimmed with Trim Galore, aligned and filtered with HiCUP, and processed into contact matrices with HiFive (BWA-backed alignment) and pairtools. A/B compartments assigned with FanC using GC content as the phasing track. Significant contacts and chromatin loops identified with FitHiC.

**CUT&Tag processing** — Reads aligned with Bowtie2 and processed with SAMtools. Peaks called with SEACR (normalized mode), benchmarked against MACS3.

**Overlap analysis** — Significance of overlap between Hi-C contacts/loops, CUT&Tag peaks, and parasitism gene islands assessed via Monte Carlo circular genomic permutation testing (fold enrichment, z-score, two-tailed p-value).

**Functional annotation** — PFAM domain enrichment among genes at significant contacts/peaks tested with Fisher's exact test and Benjamini-Hochberg correction; domain sharing across conditions visualized with UpSet plots.

## Dependencies

**Hi-C:** `Trim Galore` `HiCUP` `HiFive` `BWA` `pairtools` `cooler` `cooltools` `FanC` `FitHiC`
**CUT&Tag:** `Bowtie2` `SAMtools` `SEACR` `MACS3`
**Analysis:** `Python` `pandas` `scipy` `matplotlib` `UpSetPlot`

## Author

Griffin Kramer — [linkedin.com/in/griffin-d-kramer](https://www.linkedin.com/in/griffin-d-kramer)
