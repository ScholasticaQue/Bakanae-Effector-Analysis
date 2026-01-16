# Bakanae-Effector-Analysis
A comparative genomics and phylogeographic pipeline for rice bakanae disease. Features secretome filtering (SignalP/TargetP/EffectorP), PHI-base virulence prioritization, and OrthoFinder-based pangenome clustering coupled with BEAST2 diffusion analysis to model global pathogen spread. By analysing 13 isolates from 10 countries, including Ivory Coast, India, and Japan.

# Computational Pipeline
A multi-stage bioinformatics pipeline was developed and implemented to identify high-confidence effectors.
1. Secretome Prediction: Protein filtering using SignalP, TargetP, and DeepTMHMM.

2. Effector Identification: Using EffectorP 3.0 and physical property filters (>=4 cysteines, length 50-300aa)

3. Functional Annotation: Cross-referencing against PHI-base and PFAM databases.

4. Comparative Genomics: Clustered presence/absence analysis across 13 global isolates.
