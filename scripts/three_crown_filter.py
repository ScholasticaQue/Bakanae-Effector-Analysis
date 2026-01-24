import os
from pathlib import Path
from Bio import SeqIO

# --- PATH CONFIGURATION ---
# Base directory where all your analysis folders sit
BASE_DIR = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input")

# Paths to the tool-specific result folders
SIGNALP_DIR = BASE_DIR / "signalp_results"
TARGETP_DIR = BASE_DIR / "targetp_results"
DEEPTMHMM_DIR = BASE_DIR / "deeptmhmm_final"

# Where to save the final EffectorP input files
OUTPUT_DIR = BASE_DIR / "effectorp_input"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Your 13 Fusarium isolates
ISOLATES = ["5B", "5C", "B20", "C1995", "CFF1", "CFF2", "E282", "FSU48", 
            "IMI58289", "KSU3368", "M567", "MRC2276", "NCIM1100"]

def get_signalp_passers(isolate):
    """Path: signalp_results/{isolate}/prediction_results.txt"""
    passers = set()
    file_path = SIGNALP_DIR / isolate / "prediction_results.txt"
    if file_path.exists():
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip(): continue
                parts = line.split()
                # Col 0: ID, Col 1: Prediction (Sec/SPI)
                if len(parts) > 1 and "Sec/SPI" in parts[1]:
                    passers.add(parts[0])
    return passers

def get_targetp_passers(isolate):
    """Path: targetp_results/{isolate}/{isolate}_summary.targetp2"""
    passers = set()
    file_path = TARGETP_DIR / isolate / f"{isolate}_summary.targetp2"
    if file_path.exists():
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip(): continue
                parts = line.split()
                # Col 0: ID, Col 1: Prediction (SP)
                if len(parts) > 1 and parts[1] == "SP":
                    passers.add(parts[0])
    return passers

def get_deeptmhmm_passers(isolate):
    """Criteria: Header contains 'GLOB' (Globular/No TM helices)"""
    passers = set()
    file_path = DEEPTMHMM_DIR / isolate / f"{isolate}_merged.3line"
    if file_path.exists():
        with open(file_path, 'r') as f:
            for line in f:
                # We only need to check the header lines
                if line.startswith(">"):
                    # Example line: >5B_000001 | GLOB
                    if "GLOB" in line:
                        # Extract the ID before the first space
                        protein_id = line.replace(">", "").split()[0]
                        passers.add(protein_id)
    return passers

# --- MAIN EXECUTION ---
for iso in ISOLATES:
    print(f"[+] Processing Isolate: {iso}")
    
    # 1. Collect IDs that passed each tool
    sp_ids = get_signalp_passers(iso)
    tp_ids = get_targetp_passers(iso)
    tm_ids = get_deeptmhmm_passers(iso)
    
    # 2. Find the Intersection (Triple-Pass)
    secretome_ids = sp_ids.intersection(tp_ids).intersection(tm_ids)
    
    # 3. Extract sequences from original sanitized FASTA
    # Assuming sanitized files are in the BASE_DIR
    original_fasta = BASE_DIR / f"{iso}.faa"
    output_fasta = OUTPUT_DIR / f"{iso}.secretome.faa"
    
    secretome_records = []
    if original_fasta.exists():
        for record in SeqIO.parse(original_fasta, "fasta"):
            if record.id in secretome_ids:
                secretome_records.append(record)
        
        SeqIO.write(secretome_records, output_fasta, "fasta")
        
        print(f"    - Found in SignalP: {len(sp_ids)}")
        print(f"    - Found in TargetP: {len(tp_ids)}")
        print(f"    - Found in DeepTMHMM: {len(tm_ids)}")
        print(f"    - Final Secretome: {len(secretome_records)} sequences -> {output_fasta.name}")
        print("-" * 40)
    else:
        print(f"    [!] Error: Original FASTA {iso}.sanitized.faa not found in {BASE_DIR}")

print("\nDone! All isolate secretomes are ready in 'effectorp_input' folder.")

