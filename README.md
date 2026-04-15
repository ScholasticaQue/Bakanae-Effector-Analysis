# Bakanae-Effector-Analysis
This project presents an integrated genomic and evolutionary framework for identifying, characterizing, and prioritizing high-risk effector genes in Fusarium fujikuroi across 14 global isolates.

By combining secretome prediction, effector filtering, orthogroup analysis, and molecular evolution (dN/dS), this work reveals lineage-specific adaptation and uncovers candidate virulence factors driving pathogen evolution.
## Background
Fusarium fujikuroi is a major seed-borne pathogen causing Bakanae disease in rice. Understanding its effector repertoire is critical for:

- Seed health surveillance
- Pathogen risk assessment
- Biosecurity in rice production systems

## Objectives
- Identify high-confidence effector candidates using a consensus secretome pipeline
- Characterize lineage-specific vs core effectors
- Quantify evolutionary pressure using dN/dS (PAML CodeML)
- Perform structural and functional annotation
- Prioritize high-risk virulence candidates for surveillance

## Analysis Overview
A multi-stage bioinformatics pipeline was developed and implemented to identify high-confidence effectors.

**1. Secretome Prediction:** Protein filtering using SignalP, TargetP, and DeepTMHMM.

**2. Effector Identification:** Using EffectorP 3.0 and physical property filters (>=4 cysteines, length 50-300aa)

**3. Evolutionary Analysis:** OrthoFinder → orthogroups, MAFFT → alignment, PAL2NAL → codon alignment, PAML (CodeML) → dN/dS

**4. Structural & Functional Annotation:** AlphaFold → structure prediction, Foldseek → structural similarity, PHI-base → virulence annotation.

**5. Genomic Context Analysis:** EDTA → transposable elements, Bedtools → TE proximity, Statistical testing: Mann-Whitney U, Fisher’s Exact Test



## Significant Results Summary

This repository highlights the unique evolutionary and pathogenic adaptation strategy of the **Kenyan lineage of *Fusarium fujikuroi***.

---

#### 1. Effector Expansion in Specific Lineages
5B, 5C, Baka → **~250 effectors each**

Other isolates → **5–10 effectors**

Indicates lineage-specific pathogenic expansion

### 2. No Core Effector Repertoire
- 0 core effectors across all isolates

Suggests:

No universal infection strategy — highly adaptive pathogen

### 3. Strong Positive Selection
- 24 orthogroups under selection

- Some with **ω ≈ 999 (extreme evolution)**

### 4. Structural Virulence Factors Identified

Examples:

- Necrosis-inducing proteins (NLPs)

- Ecp2 effector

- Pectinesterases

- LPMOs

### 5. High Novelty

- 62.5% of high-risk candidates = no PHI-base match

Indicates:

Potential novel virulence mechanisms

### 6. No Genome Compartmentalization
- No TE proximity bias

- No gene sparsity differences

Evolution occurs:

Across the whole genome, not just hotspots.


## Repository Structure
The project is organized into modular directories corresponding to each major analytical step.

- `/scripts/secretome/` (SignalP, TargetP, cysteine filtering)
- `/scripts/pangenome/` (OrthoFinder scripts and core genome extraction)
- `/scripts/beast/` (BEAST2 XML configuration and log processing)
- `/scripts/viz/` (R code for the **Chord Diagram** and **Presence/Absence Heatmap**)

