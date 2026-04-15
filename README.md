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
5B, 5C, Baka → ~250 effectors each
Other isolates → 5–10 effectors

Indicates lineage-specific pathogenic expansion

### 🇰🇪 The Kenyan Effector Expansion

#### Isolate Redundancy
- A **25× expansion** of classical **SCRP (Small Cysteine-Rich Protein) effectors** was identified in Kenyan isolates **5B** and **5C** compared to global isolates.
- This suggests a lineage-specific adaptation strategy potentially enhancing host interaction and virulence potential.

#### Paralog Discovery
- **CD-HIT clustering (90% sequence identity)** reduced:
  - **1,377 raw effector hits → 288 unique effector clusters**
- Results indicate that the observed expansion is primarily driven by:
  - High paralog redundancy  
  - Local gene duplication events  

#### The "Borrowed Arsenal"
- Several effector clusters (**Cluster 102, 100, 104**) show strong homology to effectors found in:
  - *Fusarium equiseti*
  - *Fusarium oxysporum*
- These findings suggest possible:
  - Horizontal gene transfer
  - Shared ancestral effector repertoires
  - Adaptive gene acquisition

---

### Structural Homology Alignment

#### Stable Evolutionary Hub — Scaffold 10
- Enriched in **But2 domain proteins**
- Functional relevance:
  - Protein turnover
  - Cellular homeostasis
  - Potential regulatory stability within the genome

#### Plastic Adaptive Hub — Scaffold 7
- Enriched in **Hydrophobin proteins**
- Functional relevance:
  - Spore dispersal and environmental survival
  - Surface interaction and host colonization
  - Rapid environmental adaptability

---

### Biological Interpretation

Together, these results suggest that the Kenyan lineage combines:

- **Effector gene amplification**
- **Cross-species effector acquisition**
- **Genome compartmentalization into stable and adaptive hubs**

This dual evolutionary strategy may contribute to enhanced pathogenic success and ecological flexibility.


## Repository Structure
The project is organized into modular directories corresponding to each major analytical step.

- `/scripts/secretome/` (SignalP, TargetP, cysteine filtering)
- `/scripts/pangenome/` (OrthoFinder scripts and core genome extraction)
- `/scripts/beast/` (BEAST2 XML configuration and log processing)
- `/scripts/viz/` (R code for the **Chord Diagram** and **Presence/Absence Heatmap**)

