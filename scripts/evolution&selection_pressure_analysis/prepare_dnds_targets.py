import os
from Bio import SeqIO
from pathlib import Path

# --- CONFIGURATION (Adjusted to your paths) ---
base_path = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results")
clstr_file = base_path / "effectors_clustered.clstr"
fasta_file = base_path / "refined_classical_effectors.fasta"
output_dir = base_path / "dnds_input_groups"

os.makedirs(output_dir, exist_ok=True)

# --- 1. Load all sequences into memory ---
print("Loading sequences...")
seq_records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# --- 2. Parse CD-HIT .clstr file to map Cluster -> SeqIDs ---
print("Parsing clusters...")
clusters = {}
current_cluster = ""

with open(clstr_file) as f:
    for line in f:
        if line.startswith(">Cluster"):
            current_cluster = line.strip().replace(">", "").replace(" ", "_")
            clusters[current_cluster] = []
        else:
            # Extract the ID between '>' and '...'
            # Format is usually: >Isolate|ID...
            full_id = line.split(">")[1].split("...")[0]
            clusters[current_cluster].append(full_id)

# --- 3. Identify "Mixed" Clusters and Extract Sequences ---
kenya_isolates = {"5B", "5C"}
extracted_count = 0

print(f"Checking {len(clusters)} clusters for dN/dS criteria...")

for cluster_id, member_ids in clusters.items():
    # Identify which isolates are in this cluster
    isolates_in_cluster = {mid.split("|")[0] for mid in member_ids}
    
    has_kenya = any(iso in kenya_isolates for iso in isolates_in_cluster)
    has_global = any(iso not in kenya_isolates for iso in isolates_in_cluster)
    
    # Selection analysis needs at least 2 sequences to compare
    if has_kenya and has_global and len(member_ids) >= 2:
        output_file = output_dir / f"{cluster_id}.fasta"
        
        cluster_records = []
        for mid in member_ids:
            if mid in seq_records:
                cluster_records.append(seq_records[mid])
            else:
                print(f"Warning: {mid} not found in FASTA!")

        if cluster_records:
            SeqIO.write(cluster_records, output_file, "fasta")
            extracted_count += 1

print("-" * 30)
print(f"Success! Extracted {extracted_count} mixed clusters to:")
print(f" {output_dir}")

