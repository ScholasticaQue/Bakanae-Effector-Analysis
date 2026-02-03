"""
Pathogen Effector Analysis: Evolutionary Selection Pressure (dN/dS)
Author: Scholastica Quaicoe
Institution: Africa Rice Center

Description:
This script automates PAML/CodeML analysis for predicted effector clusters.
It compares site models M7 (null) vs M8 (selection) to identify specific 
amino acids under positive selection (ω > 1). This is critical for 
identifying signatures of host-adaptation in global isolates.
"""
import os
import subprocess
from pathlib import Path

# Paths
# Aligned codon sequences and phylogenetic trees should be in this directory
WORK_DIR = Path("analysis/evolutionary_selection/alignments")
os.chdir(work_dir)

def create_ctl(name, aln_file, tree_file):
    """Generates the .ctl configuration file required for PAML Codeml."""
    ctl_content = f"""
      seqfile = {aln_file}
      treefile = {tree_file}
      outfile = {name}.sites.txt
      noisy = 0
      verbose = 0
      runmode = 0
      seqtype = 1
      CodonFreq = 2
      model = 0
      NSsites = 7 8    * Run both M7 (null) and M8 (selection)
      icode = 0
      fix_kappa = 0
      kappa = 2
      fix_omega = 0
      omega = 0.5
      ncatG = 10
    """
    with open(f"{name}.ctl", "w") as f:
        f.write(ctl_content)

print("🧬 Starting dN/dS Analysis...")

for aln in work_dir.glob("*.codon.aln"):
    # Fix: Correctly extract just 'Cluster_XX' by removing '.codon.aln'
    base_name = aln.name.replace(".codon.aln", "")
    tree = work_dir / f"{base_name}.tree"
    
    if tree.exists():
        print(f" Running CodeML for {base_name}...")
        create_ctl(base_name, aln.name, tree.name)
        
        # Run codeml
        subprocess.run(["codeml", f"{base_name}.ctl"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # Clean up PAML junk
        for temp in ["2NG.dn", "2NG.ds", "2NG.t", "4fold.nuc", "lnf", "rst", "rst1", "rub"]:
            if os.path.exists(temp): os.remove(temp)
    else:
        # Debug print to see exactly what it was looking for
        print(f" Missing tree: expected {tree.name}, but not found.")

print(" Analysis finished.")
