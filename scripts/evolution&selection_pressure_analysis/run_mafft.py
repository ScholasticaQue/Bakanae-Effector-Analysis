import os
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path

# ==============================
# PATHS & CONFIGURATION
# ==============================
TSV_FILE = "/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results/mapped_effector_orthogroups.tsv"
TXT_FILE = "/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results/refined_effector_ids.txt"
PROT_FASTA = "/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results/refined_classical_effectors.fasta"

# --- NEW: Path to the directory containing individual isolate CDS files ---
# Change this to the actual folder path where 5B.fasta, 5C.fasta, etc., are stored
CDS_INPUT_DIR = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input")

# Output Directories
DIR_190_PROT = Path("orthogroups_5B_5C_Baka_aligned")
DIR_16_PROT = Path("orthogroups_presence_5plus_aligned")
DIR_190_CDS = Path("orthogroups_5B_5C_Baka_CDS")
DIR_16_CDS = Path("orthogroups_presence_5plus_CDS")

for folder in [DIR_190_PROT, DIR_16_PROT, DIR_190_CDS, DIR_16_CDS]:
    folder.mkdir(parents=True, exist_ok=True)

ISOLATES = ["5B", "5C", "Baka", "B20", "C1995", "CFF1", "CFF2",
            "E282", "FSU48", "IMI58289", "KSU3368", "M567",
            "MRC2276", "NCIM1100"]

# ==============================
# HELPER FUNCTIONS
# ==============================
def normalize_id(gene_id):
    gene_id = str(gene_id).strip().replace('"', '')
    return gene_id.split("|")[-1] if "|" in gene_id else gene_id

def align_with_mafft(input_fasta, output_fasta):
    try:
        cmd = f"mafft --auto {input_fasta} > {output_fasta}"
        subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error aligning {input_fasta}: {e}")

# ==============================
# 1. LOAD DATA
# ==============================
print("Loading sequence data...")

# Load valid effector list
with open(TXT_FILE) as f:
    valid_effectors = {normalize_id(line) for line in f if line.strip()}

# Load Protein Sequences
prot_dict = {normalize_id(rec.id): rec for rec in SeqIO.parse(PROT_FASTA, "fasta")}

# --- NEW: Load CDS sequences from individual Isolate files ---
cds_dict = {}
for iso in ISOLATES:
    # Adjust file extension if necessary (e.g., .fna, .cds, .fasta)
    iso_cds_file = CDS_INPUT_DIR / f"{iso}.cds.fasta" 
    
    if iso_cds_file.exists():
        for rec in SeqIO.parse(iso_cds_file, "fasta"):
            clean_id = normalize_id(rec.id)
            if clean_id in valid_effectors:
                cds_dict[clean_id] = rec
    else:
        print(f"⚠️ Warning: Could not find CDS file for isolate {iso}: {iso_cds_file}")

print(f"Loaded {len(cds_dict)} total CDS sequences across all isolates.")

# Load Orthogroup TSV
df = pd.read_csv(TSV_FILE, sep="\t", dtype=str).fillna("")

# ==============================
# 2. PROCESS, ALIGN & EXTRACT CDS
# ==============================
print("Processing orthogroups...")

for _, row in df.iterrows():
    og_id = row["Orthogroup"]
    og_prot_records = []
    og_cds_records = []
    present_isolates = []

    for iso in ISOLATES:
        cell = row.get(iso, "")
        if cell:
            genes = [normalize_id(g) for g in str(cell).split(",") if g.strip()]
            for g in genes:
                if g in valid_effectors:
                    if g in prot_dict:
                        og_prot_records.append(prot_dict[g])
                    if g in cds_dict:
                        og_cds_records.append(cds_dict[g])
                    
                    if iso not in present_isolates:
                        present_isolates.append(iso)

    if not og_prot_records:
        continue

    # Determine selection criteria
    targets = []
    if all(iso in present_isolates for iso in ["5B", "5C", "Baka"]):
        targets.append((DIR_190_PROT, DIR_190_CDS))
    if len(present_isolates) >= 5:
        targets.append((DIR_16_PROT, DIR_16_CDS))

    for prot_folder, cds_folder in targets:
        # Protein Alignment
        temp_prot = prot_folder / f"{og_id}_temp.fasta"
        aln_prot = prot_folder / f"{og_id}.fasta"
        SeqIO.write(og_prot_records, temp_prot, "fasta")
        
        if len(og_prot_records) > 1:
            align_with_mafft(temp_prot, aln_prot)
            os.remove(temp_prot)
        else:
            os.rename(temp_prot, aln_prot)

        # CDS Extraction
        if og_cds_records:
            cds_out = cds_folder / f"{og_id}_cds.fasta"
            SeqIO.write(og_cds_records, cds_out, "fasta")

print("-" * 30)
print(f"✅ Analysis Complete.")
