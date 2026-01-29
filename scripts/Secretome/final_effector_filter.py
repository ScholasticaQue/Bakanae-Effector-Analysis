import os
from Bio import SeqIO
from pathlib import Path

# --- PATHS ---
input_dir = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_input")
results_dir = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results")
results_dir.mkdir(parents=True, exist_ok=True)

# New Output Name to reflect filtering
master_fasta = results_dir / "refined_classical_effectors.fasta"

ISOLATES = ["5B", "5C", "B20", "C1995", "CFF1", "CFF2", "E282", "FSU48", 
            "IMI58289", "KSU3368", "M567", "MRC2276", "NCIM1100"]

# --- FILTERING FUNCTION ---
def is_classical_effector(seq_record):
    sequence = str(seq_record.seq)
    length = len(sequence)
    cys_count = sequence.upper().count("C")
    
    # Criteria: 50-300 amino acids AND at least 4 Cysteines
    if 50 <= length <= 300 and cys_count >= 4:
        return True
    return False

# Load valid effector IDs per isolate (EffectorP results)
effector_map = {}
for iso in ISOLATES:
    res_file = results_dir / f"{iso}_effectorp_results.txt"
    if res_file.exists():
        with open(res_file) as f:
            effector_map[iso] = {line.split()[0] for line in f 
                                if "effector" in line.lower() and "non-effector" not in line.lower()}

# Extract and Filter sequences
filtered_records = []
total_effectorp_hits = 0

for iso in ISOLATES:
    fasta_file = input_dir / f"{iso}.secretome.faa"
    if fasta_file.exists():
        iso_matches = 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in effector_map.get(iso, set()):
                total_effectorp_hits += 1
                # --- NEW FILTERING STEP ---
                if is_classical_effector(record):
                    record.id = f"{iso}|{record.id}"
                    record.description = ""
                    filtered_records.append(record)
                    iso_matches += 1
        print(f"Isolate {iso}: {iso_matches} classical effectors passed filters.")
    else:
        print(f"⚠️ Warning: Could not find fasta for {iso}")

if filtered_records:
    SeqIO.write(filtered_records, master_fasta, "fasta")
    print("-" * 30)
    print(f"📊 Summary Statistics:")
    print(f"Total EffectorP raw hits: {total_effectorp_hits}")
    print(f"Refined Classical Effectors: {len(filtered_records)}")
    print(f"📍 Location: {master_fasta}")
else:
    print("❌ Error: No sequences passed the filters. Adjust criteria if necessary.")

