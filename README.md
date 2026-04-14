# Bakanae-Effector-Analysis
This project presents an integrated genomic and evolutionary framework for identifying, characterizing, and prioritizing high-risk effector genes in Fusarium fujikuroi across 114 global isolates.

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

## Computational Pipeline
A multi-stage bioinformatics pipeline was developed and implemented to identify high-confidence effectors.

**1. Secretome Prediction:** Protein filtering using SignalP, TargetP, and DeepTMHMM.

**2. Effector Identification:** Using EffectorP 3.0 and physical property filters (>=4 cysteines, length 50-300aa)

**3. Functional Annotation:** Cross-referencing against PHI-base and PFAM databases.

**4. Comparative Genomics:** Clustered presence/absence analysis across 13 global isolates.

## Evolutionary & Phylogeographic Diffusion Modeling

To investigate the global movement and evolutionary history of *Fusarium fujikuroi*, a phylogeographic reconstruction framework was implemented as outlined below.

### 1. Pangenome Clustering
Protein-coding sequences (translated CDS) from all 13 *F. fujikuroi* isolates were clustered into orthologous groups using **OrthoFinder**.  
This step enabled the identification of shared and accessory gene content across isolates.

- Input: Translated CDS FASTA files  
- Tool: OrthoFinder  
- Output: Orthogroups and inferred gene families

### 2. Core Genome Alignment
Orthogroups present in all isolates were defined as the **core genome**.  
These conserved genes were extracted and concatenated to construct a high-confidence alignment representing the evolutionary backbone of the species.

- Criteria: Single-copy orthologs present in all 13 isolates  
- Output: Concatenated core genome alignment

### 3. Bayesian Diffusion Modeling
A Bayesian phylogeographic analysis was conducted using **BEAST2** to reconstruct spatial and temporal diffusion dynamics.

- Input:
  - Core genome alignment
  - Sampling dates (temporal data)
  - Geographic origin metadata (10 countries)
- Model:
  - Discrete trait diffusion model
  - Bayesian MCMC framework
- Output:
  - Time-calibrated phylogeny
  - Ancestral location state reconstructions
  - BEAST log and tree files

### 4. Pathogen Connectivity Visualization
To visualize inferred migration patterns, **diffusion events** between geographical regions were extracted from BEAST log files.  
A **Chord Diagram** was generated to represent the intensity and directionality of transitions between countries.

- Input: BEAST diffusion logs
- Visualization: Chord diagram
- Interpretation: Strength and frequency of pathogen movement between regions

## Evolutionary and Selection Pressure Analysis
Evolutionary pressure ($dN/dS$) was analyzed across predicted effector proteins in the 13 Fusarium fujikuroi isolates. 
The goal was to identify signatures of host-adaptation and positive selection within the fungal secretome. By comparing the rates of non-synonymous ($dN$) and synonymous ($dS$) mutations, we characterized the selective forces acting on fungal "weapons" (effectors) used during plant infection.

**1. Clustering:** Homologous sequences were grouped using **CD-HIT** (90% identity threshold) to define effector "clusters" across isolates.

**2. Codon-Aware Alignment:** Nucleotide sequences (CDS) were aligned using **MAFFT** and **PAL2NAL**, ensuring the preservation of the triplet codon reading frame.

**3. Phylogeny:** Maximum Likelihood trees were generated for each cluster using FastTree (GTR model).

**4. Selection Testing:** Evolutionary rates were calculated using PAML (CodeML), comparing:
- **M0 (One-ratio):** To determine the global $\omega$ per cluster.
- **M7 vs M8 (Site models):** To detect specific amino acids under positive selection.

### Functional Prioritization

To refine candidate effector targets and predict their biological relevance in host infection, structural bioinformatics approaches were applied.

#### Structural Validation
- Predicted three-dimensional protein structures using **AlphaFold3** to obtain high-confidence structural models of prioritized effector candidates.
- Structural predictions enabled:
  - Identification of conserved functional folds
  - Detection of structural motifs not apparent from sequence analysis alone
  - Improved confidence in effector annotation

#### Structural Homology Screening
- Conducted structural similarity searches using **FoldSeek** to compare predicted effector structures against known protein databases.
- This analysis facilitated:
  - Identification of potential functional analogs
  - Discovery of conserved virulence-associated structural domains
  - Detection of distant evolutionary relationships beyond sequence similarity

#### Host Target Inference
- Combined structural predictions and homology results to prioritize effectors likely to interact with **rice host biological pathways**.
- This approach supports:
  - Identification of candidate host manipulation strategies
  - Selection of high-value targets for downstream functional validation
  - Improved understanding of pathogen–host molecular interactions

---

### Biological Significance

Structural prioritization provides a higher-resolution framework for identifying candidate virulence factors and helps bridge genomic discovery with functional pathogenicity mechanisms.


## Significant Results Summary

This repository highlights the unique evolutionary and pathogenic adaptation strategy of the **Kenyan lineage of *Fusarium fujikuroi***.

---

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

