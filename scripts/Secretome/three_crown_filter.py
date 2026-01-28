import os
import re
from pathlib import Path
from Bio import SeqIO

# --- PATH CONFIGURATION ---
BASE_DIR = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input")
SIGNALP_DIR = BASE_DIR / "signalp_results"
TARGETP_DIR = BASE_DIR / "targetp_results"
DEEPTMHMM_DIR = BASE_DIR / "deeptmhmm_final"
OUTPUT_DIR = BASE_DIR / "effectorp_input"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

ISOLATES = ["5B", "5C", "B20", "C1995", "CFF1", "CFF2", "E282", "FSU48", 
            "IMI58289", "KSU3368", "M567", "MRC2276", "NCIM1100"]

def get_standardized_id(raw_id, iso_prefix):
    """
    Standardizes: 'B20_num=7', 'B20_7', 'B20_000007' -> 'B20_7'
    """
    raw_id = str(raw_id).strip().replace(">", "")
    # Find all numbers in the string
    numbers = re.findall(r'\d+', raw_id)
    if numbers:
        # Use the isolate prefix + the numeric value (as an integer to remove leading zeros)
        return f"{iso_prefix}_{int(numbers[-1])}"
    return raw_id

def get_signalp_passers(isolate):
    passers = set()
    file_path = SIGNALP_DIR / isolate / "prediction_results.txt"
    if file_path.exists():
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip(): continue
                parts = line.split()
                if len(parts) > 1 and parts[1] == "SP":
                    passers.add(get_standardized_id(parts[0], isolate))
    return passers

def get_targetp_passers(isolate):
    passers = set()
    # Check for both {isolate}_summary.targetp2 and {isolate}.targetp2 (Fixes KSU3368)
    folder = TARGETP_DIR / isolate
    files = list(folder.glob("*.targetp2")) if folder.exists() else []
    
    for file_path in files:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip(): continue
                parts = line.split()
                if len(parts) > 1 and parts[1] == "SP":
                    passers.add(get_standardized_id(parts[0], isolate))
    return passers

def get_deeptmhmm_passers(isolate):
    passers = set()
    folder = DEEPTMHMM_DIR / isolate
    files = list(folder.glob("*.3line")) if folder.exists() else []
    
    for file_path in files:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith(">") and ("GLOB" in line or "SP" in line):
                    passers.add(get_standardized_id(line, isolate))
    return passers

# --- MAIN EXECUTION ---
print(f"{'Isolate':<12} | {'SignalP':<8} | {'TargetP':<8} | {'DeepTMHMM':<9} | {'Final'}")
print("-" * 70)

for iso in ISOLATES:
    sp_ids = get_signalp_passers(iso)
    tp_ids = get_targetp_passers(iso)
    tm_ids = get_deeptmhmm_passers(iso)
    
    secretome_ids = sp_ids.intersection(tp_ids).intersection(tm_ids)
    
    # Locate FASTA (Check .sanitized.faa then .faa)
    fasta_path = BASE_DIR / f"{iso}.sanitized.faa"
    if not fasta_path.exists():
        fasta_path = BASE_DIR / f"{iso}.faa"
        
    output_fasta = OUTPUT_DIR / f"{iso}.secretome.faa"
    final_records = []
    
    if fasta_path.exists():
        for record in SeqIO.parse(fasta_path, "fasta"):
            if get_standardized_id(record.id, iso) in secretome_ids:
                final_records.append(record)
        
        if final_records:
            SeqIO.write(final_records, output_fasta, "fasta")
            
        print(f"{iso:<12} | {len(sp_ids):<8} | {len(tp_ids):<8} | {len(tm_ids):<9} | {len(final_records)}")
    else:
        print(f"{iso:<12} | [!] FASTA MISSING")

print("\nDone! Files generated in 'effectorp_input'.")
