# Bakanae-Effector-Analysis
A comparative genomics and phylogeographic pipeline for rice bakanae disease. Features secretome filtering (SignalP/TargetP/EffectorP), PHI-base virulence prioritization, and OrthoFinder-based pangenome clustering coupled with BEAST2 diffusion analysis to model global pathogen spread, by analysing 13 isolates from 10 countries, including the Ivory Coast, India, and Japan.

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

**Clustering:** Homologous sequences were grouped using **CD-HIT** (90% identity threshold) to define effector "clusters" across isolates.

**Codon-Aware Alignment:** Nucleotide sequences (CDS) were aligned using **MAFFT** and **PAL2NAL**, ensuring the preservation of the triplet codon reading frame.

**Phylogeny:** Maximum Likelihood trees were generated for each cluster using FastTree (GTR model).

**Selection Testing:** Evolutionary rates were calculated using PAML (CodeML), comparing:
- M0 (One-ratio): To determine the global $\omega$ per cluster.
- M7 vs M8 (Site models): To detect specific amino acids under positive selection.

## Repository Structure
The project is organized into modular directories corresponding to each major analytical step.

- `/scripts/secretome/` (SignalP, TargetP, cysteine filtering)
- `/scripts/pangenome/` (OrthoFinder scripts and core genome extraction)
- `/scripts/beast/` (BEAST2 XML configuration and log processing)
- `/scripts/viz/` (R code for the **Chord Diagram** and **Presence/Absence Heatmap**)

