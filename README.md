CONFADS2C (v2.1)Universal Interactive Conformational Search & Quantum Chemistry Input Preparation WorkflowCONFADS2C is an automated, interactive Python workflow designed to bridge fast semiempirical conformational searching with high-level Quantum Mechanics (QM) calculations. It streamlines the transition from raw SMILES strings to cluster-ready ORCA input files and SLURM job submission scripts.🌟 Key FeaturesUniversal Molecular Support: Handles neutral organic molecules, cations, anions, radicals, open-shell systems, and metal-containing complexes.Automatic Electronic State Detection: Computes formal charge ($q$) and estimates default spin multiplicity ($2S + 1$) directly from SMILES, with interactive confirmation/override options.In-Memory 3D Embedding: Employs RDKit's ETKDGv3 algorithm coupled with MMFF94/UFF force field pre-minimization.xTB & CREST Integration: Fast geometry pre-optimization with GFN2-xTB followed by extensive conformational ensemble sampling using CREST.Machine-Learned Prioritization: Integrated with CONFPASS (Goodman Lab, Cambridge) for structural clustering and priority ranking of conformers.Interactive CLI & Smart Presets: Prompts for ORCA theory levels (e.g., B3LYP-D3BJ, $\omega$B97X-D3, r2SCAN-3c, or custom lines), memory allocation (%maxcore), and SLURM parameters.Cluster-Ready Output Generation: Automatically formats individual .xyz geometries, .inp ORCA files, and executable .sh SLURM scripts for batch execution (sbatch).Detailed Execution Summary: Generates performance reports (*_summary.txt) tracking time elapsed per step and total runtime percentages.🔄 Workflow Pipeline Architecture                       ┌──────────────────────────────┐
                       │     SMILES Input / Parsing   │
                       └──────────────┬───────────────┘
                                      │
                                      ▼
                       ┌──────────────────────────────┐
                       │  RDKit 3D Embedding (ETKDGv3)│
                       └──────────────┬───────────────┘
                                      │
                                      ▼
                       ┌──────────────────────────────┐
                       │   xTB Pre-Optimization (GFN2)│
                       └──────────────┬───────────────┘
                                      │
                                      ▼
                       ┌──────────────────────────────┐
                       │ CREST Conformational Search  │
                       └──────────────┬───────────────┘
                                      │
                                      ▼
                       ┌──────────────────────────────┐
                       │  OpenBabel Format Conversion │
                       │          (XYZ ──► SDF)       │
                       └──────────────┬───────────────┘
                                      │
                                      ▼
                       ┌──────────────────────────────┐
                       │  CONFPASS Ensemble Ranking   │
                       │     & Priority Extraction    │
                       └──────────────┬───────────────┘
                                      │
                                      ▼
                       ┌──────────────────────────────┐
                       │ Top Conformer Selection &    │
                       │   ORCA / SLURM Generation    │
                       └──────────────────────────────┘
🛠️ Prerequisites & Installation1. Debian / Ubuntu Native SetupTo install all standard operating system dependencies via standard APT packaging:sudo apt update && sudo apt install -y \
  python3-rdkit \
  openbabel \
  python3-openbabel \
  python3-pandas \
  python3-numpy \
  python3-sklearn \
  python3-natsort \
  python3-networkx \
  python3-scipy \
  git wget tar
2. Install Grimme Lab Binaries (xTB & CREST)Download the official static binaries compiled with OpenMP support:# Install xtb (v6.7.1)
wget https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz
tar -xf xtb-6.7.1-linux-x86_64.tar.xz
sudo cp xtb-dist/bin/xtb /usr/local/bin/
rm -rf xtb-dist xtb-6.7.1-linux-x86_64.tar.xz

# Install CREST (v3.0.2)
wget https://github.com/crest-lab/crest/releases/download/v3.0.2/crest-gnu-12-ubuntu-latest.tar.xz
tar -xf crest-gnu-12-ubuntu-latest.tar.xz
sudo cp crest/crest /usr/local/bin/
sudo chmod +x /usr/local/bin/crest
rm -rf crest crest-gnu-12-ubuntu-latest.tar.xz
Verify binary installations:xtb --version && crest --version
3. Install CONFPASSClone and install the CONFPASS library into your system Python environment:git clone https://github.com/Goodman-lab/CONFPASS.git /tmp/CONFPASS
SITE_DIR=$(python3 -c "import site; print(site.getsitepackages()[0])")
sudo mkdir -p $SITE_DIR/confpass
sudo cp -r /tmp/CONFPASS/* $SITE_DIR/confpass/
rm -rf /tmp/CONFPASS
Verify CONFPASS module functionality:python3 -m confpass -h
🚀 UsageRun the workflow script from your terminal:python3 confads2c.v2.1.py
Interactive Example Session======================================================================
          CONFADS2C v2.1 — UNIVERSAL CONFORMATIONAL SEARCH            
======================================================================

Enter SMILES string: CC(=O)Oc1ccccc1C(=O)O

  ℹ Molecule Analysis:
    • Formula           : C9H8O4
    • Molecular Weight  : 180.16 g/mol
    • Heavy Atom Count  : 13
    • Detected Charge   : 0
    • Est. Multiplicity : 1 (Singlet/Closed-Shell)

Enter identifier name [Default: C9H8O4]: Aspirin

Physical / Electronic State Settings:
 Confirm/Override Formal Charge [0]: 0
 Confirm/Override Spin Multiplicity (1=Singlet, 2=Doublet, 3=Triplet) [1]: 1

Conformational Search Settings:
 Number of top conformers to extract for ORCA [Default: 10]: 5
 CREST Energy Window (kcal/mol) [Default: 10.0]: 6.0

ORCA Quantum Chemistry Settings:
 Presets for Level of Theory:
   [1] B3LYP D3BJ def2-SVP OPT TightOPT FREQ (Fast & Standard)
   [2] ωB97X-D3 def2-TZVP OPT TightOPT FREQ (High Accuracy)
   [3] r2SCAN-3c OPT TightOPT FREQ (Modern Composite Method)
   [4] Custom input line
 Select preset [1-4, Default: 1]: 1

SLURM Cluster & Execution Settings:
 Walltime limit (HH:MM:SS) [Default: 24:00:00]: 12:00:00
 CPU Cores (--ntasks) [Default: 16]: 16
 SLURM Partition (leave blank if default): standard
 ORCA %maxcore memory per core in MB [Default: Auto/Dynamic]: 4000
📁 Output Directory StructureFor a run named Aspirin, CONFADS2C creates an isolated workspace directory with the following structure:Aspirin_workflow_20260917_160000/
├── Aspirin_initial.xyz             # Initial RDKit ETKDG 3D geometry
├── xtbopt.xyz                      # xTB pre-optimized geometry
├── crest_conformers.xyz            # Full CREST conformer ensemble
├── Aspirin_conformers.sdf          # OpenBabel converted ensemble file
├── Aspirin_confpass_output.txt     # CONFPASS ranking & priority logs
├── Aspirin_summary.txt             # Comprehensive timing & execution report
└── Aspirin_conformers/             # Output QM calculation folder
    ├── Aspirin_1.xyz
    ├── Aspirin_1.inp               # ORCA Input file
    ├── Aspirin_1.sh                # Executable SLURM submit script
    ├── Aspirin_2.xyz
    ├── Aspirin_2.inp
    ├── Aspirin_2.sh
    └── ...
Submitting Jobs to a ClusterNavigate to the generated conformers folder and submit all generated SLURM jobs at once:cd Aspirin_workflow_20260917_160000/Aspirin_conformers/
for f in *.sh; do sbatch "$f"; done
📊 Summary Report Example======================================================================
CONFADS2C v2.1 — WORKFLOW SUMMARY REPORT
======================================================================

Molecule Identifier : Aspirin
SMILES String       : CC(=O)Oc1ccccc1C(=O)O
Execution Date      : 2026-09-17 16:02:15

CHEMICAL & PHYSICAL PROPERTIES:
  • Formula           : C9H8O4
  • Formal Charge     : 0
  • Spin Multiplicity : 1

CONFORMATIONAL SEARCH RESULTS:
  • CREST Conformers  : 14 total ensemble structures
  • Extracted Top     : 5 conformers for QM re-optimization
  • Energy Window     : 6.0 kcal/mol

ORCA & SLURM CONFIGURATION:
  • Level of Theory   : ! B3LYP D3BJ def2-SVP OPT TightOPT FREQ
  • Parallel Cores    : 16 cores (%pal nprocs 16)
  • Memory per Core   : 4000 MB
  • SLURM Walltime    : 12:00:00

TIMING BREAKDOWN:
  • 3D Embedding (RDKit)       :    0.18 s (  0.5%)
  • xTB Optimization           :    4.12 s ( 11.2%)
  • CREST Search               :   31.25 s ( 84.8%)
  • OpenBabel Conversion       :    0.32 s (  0.9%)
  • CONFPASS Evaluation        :    0.85 s (  2.3%)
  • Extraction & QM Inputs     :    0.12 s (  0.3%)
  --------------------------------------------------
  • TOTAL WORKFLOW TIME       :   36.84 s

======================================================================
📜 Citations & ReferencesIf you use this workflow script in academic publications, please cite the underlying tools:RDKit: RDKit: Open-source cheminformatics toolkit. https://www.rdkit.orgxTB: Grimme, S.; Bannwarth, C.; Shushkov, P. J. Chem. Theory Comput. 2017, 13 (5), 1989–2009.CREST: Pracht, P.; Bohle, F.; Grimme, S. Phys. Chem. Chem. Phys. 2020, 22 (14), 7169–7192.CONFPASS: Goodman Lab, University of Cambridge. https://github.com/Goodman-lab/CONFPASSOpenBabel: O'Boyle, N. M. et al. J. Cheminf. 2011, 3 (1), 33.ORCA: Neese, F. Wiley Interdiscip. Rev. Comput. Mol. Sci. 2012, 2 (1), 73–78.📄 LicenseDistributed under the MIT License. See LICENSE for details.
