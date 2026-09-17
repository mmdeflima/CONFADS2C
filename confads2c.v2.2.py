#!/usr/bin/env python3
"""
=============================================================================
CONFADS2C: Interactive Conformational Search & QM Preparation Workflow (v2.2)
=============================================================================
Supported Tools: RDKit, OpenBabel, xTB, CREST, CONFPASS, ORCA, SLURM
Features:
  - Universal molecule support (neutral, cations, anions, radicals, metal complexes)
  - Automatic charge & spin multiplicity detection with user confirmation
  - High-performance native 3D embedding via RDKit ETKDGv3
  - Thread safety controls preventing OpenBLAS/MKL/OpenMP log inflation bug (v2.2)
  - Interactive smart prompts with preset choices and safety checks
  - Process safety via safe subprocess execution (shell=False)
  - Detailed terminal timing, progress spinners, and output summary report
=============================================================================
"""

import os
import sys

# ------------------------------------------------------------
# Environment Safety & Thread Oversaturation Control (v2.2 Fix)
# Prevent OpenBLAS / MKL / OMP thread explosion and infinite log growth
# ------------------------------------------------------------
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_STACKSIZE"] = "1G"
os.environ["OPENBLAS_VERBOSE"] = "0"

import re
import time
import threading
import subprocess
import multiprocessing
from pathlib import Path
from datetime import datetime

# ------------------------------------------------------------
# RDKit Import Verification
# ------------------------------------------------------------
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
except ImportError:
    print("\033[0;31mERROR: RDKit is not installed in the current Python environment.\033[0m")
    print("Install via APT: sudo apt install python3-rdkit")
    print("Or activate your computational chemistry environment.")
    sys.exit(1)

# ANSI Color Formatting & Terminal Styling
GREEN = '\033[0;32m'
YELLOW = '\033[1;33m'
BLUE = '\033[0;34m'
CYAN = '\033[0;36m'
RED = '\033[0;31m'
BOLD = '\033[1m'
DIM = '\033[2m'
NC = '\033[0m'

ICON_OK = f"{GREEN}✓{NC}"
ICON_FAIL = f"{RED}✗{NC}"
ICON_WARN = f"{YELLOW}⚠{NC}"
ICON_INFO = f"{CYAN}ℹ{NC}"

class Spinner:
    """Animated progress spinner for long-running subprocess tasks."""
    def __init__(self, message=" Executing..."):
        self.message = message
        self.spinning = False
        self.thread = None

    def spin(self):
        chars = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
        idx = 0
        while self.spinning:
            sys.stdout.write(f"\r {CYAN}{chars[idx % len(chars)]}{NC} {self.message}")
            sys.stdout.flush()
            time.sleep(0.12)
            idx += 1
        sys.stdout.write("\r" + " " * (len(self.message) + 10) + "\r")
        sys.stdout.flush()

    def start(self):
        self.spinning = True
        self.thread = threading.Thread(target=self.spin, daemon=True)
        self.thread.start()

    def stop(self):
        self.spinning = False
        if self.thread:
            self.thread.join(timeout=0.5)

def run_command_safe(cmd_list, step_name, stdout_log=None, stderr_log=None, use_spinner=True):
    """Executes external commands safely as argument lists without shell injection risks."""
    print(f" {ICON_INFO} {BOLD}{step_name}{NC} starting...")
    spinner = Spinner(f"{step_name} in progress...")
    if use_spinner:
        spinner.start()
        
    start_time = time.perf_counter()
    
    stdout_target = open(stdout_log, 'w') if stdout_log else subprocess.PIPE
    stderr_target = open(stderr_log, 'w') if stderr_log else subprocess.PIPE

    try:
        process = subprocess.run(
            cmd_list,
            stdout=stdout_target,
            stderr=stderr_target,
            text=True,
            env=os.environ,
            check=False
        )
    finally:
        if stdout_log:
            stdout_target.close()
        if stderr_log:
            stderr_target.close()
            
    end_time = time.perf_counter()
    if use_spinner:
        spinner.stop()
        
    elapsed = end_time - start_time

    # Log size safety check
    if stdout_log and Path(stdout_log).exists():
        size_mb = Path(stdout_log).stat().st_size / (1024 * 1024)
        if size_mb > 500:
            print(f" {ICON_WARN} {YELLOW}Warning: Log file {stdout_log} is large ({size_mb:.1f} MB).{NC}")

    if process.returncode != 0:
        if stderr_log and Path(stderr_log).exists():
            with open(stderr_log, 'r') as err_f:
                err_content = err_f.read()
                if "Missing comma between descriptors" in err_content:
                    print(f" {ICON_WARN} {YELLOW}Notice: Minor xTB output format warning ignored.{NC}")
                    return elapsed

        print(f" {ICON_FAIL} {RED}Failed during step: {step_name}{NC}")
        if stderr_log and Path(stderr_log).exists():
            print(f"{DIM}--- Error Log ({stderr_log}) ---{NC}")
            with open(stderr_log, 'r') as err_f:
                print(err_f.read()[-1000:])
        sys.exit(1)

    print(f" {ICON_OK} {GREEN}{step_name} completed in {elapsed:.2f}s{NC}")
    return elapsed

def smiles_to_3d_xyz(mol, output_xyz_path: Path) -> tuple[int, int]:
    """Generates initial 3D coordinates using RDKit ETKDGv3 and MMFF/UFF force fields."""
    mol_with_h = Chem.AddHs(mol)
    
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    embed_res = AllChem.EmbedMolecule(mol_with_h, params)
    
    if embed_res != 0:
        AllChem.EmbedMolecule(mol_with_h, useRandomCoords=True)

    if AllChem.MMFFHasAllMoleculeParams(mol_with_h):
        AllChem.MMFFOptimizeMolecule(mol_with_h, maxIters=500)
    else:
        AllChem.UFFOptimizeMolecule(mol_with_h, maxIters=500)

    num_atoms = mol_with_h.GetNumAtoms()
    conf = mol_with_h.GetConformer()
    
    with open(output_xyz_path, 'w') as f:
        f.write(f"{num_atoms}\nInitial RDKit ETKDGv3 3D Geometry\n")
        for atom_idx in range(num_atoms):
            atom = mol_with_h.GetAtomWithIdx(atom_idx)
            pos = conf.GetAtomPosition(atom_idx)
            f.write(f"{atom.GetSymbol():<2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n")

    return num_atoms

def analyze_smiles(smiles_input: str):
    """Parses SMILES and calculates charge, formula, weight, and default spin multiplicity."""
    mol = Chem.MolFromSmiles(smiles_input)
    if mol is None:
        return None
    
    charge = Chem.GetFormalCharge(mol)
    total_electrons = sum(atom.GetAtomicNum() for atom in Chem.AddHs(mol).GetAtoms()) - charge
    spin_multiplicity = 1 if (total_electrons % 2 == 0) else 2
    
    formula = rdMolDescriptors.CalcMolFormula(Chem.AddHs(mol))
    mw = Descriptors.MolWt(mol)
    
    return {
        "mol": mol,
        "formula": formula,
        "mw": mw,
        "charge": charge,
        "spin": spin_multiplicity,
        "num_heavy_atoms": mol.GetNumHeavyAtoms()
    }

def get_user_inputs():
    """Interactive prompt collector with instant structure verification and smart defaults."""
    print(f"\n{BOLD}{CYAN}======================================================================{NC}")
    print(f"{BOLD}{CYAN}          CONFADS2C v2.2 — UNIVERSAL CONFORMATIONAL SEARCH            {NC}")
    print(f"{BOLD}{CYAN}======================================================================{NC}\n")

    while True:
        smiles = input(f"{BOLD}Enter SMILES string:{NC} ").strip()
        if not smiles:
            print(f"{RED}SMILES string cannot be empty.{NC}")
            continue
            
        analysis = analyze_smiles(smiles)
        if analysis is None:
            print(f"{RED}Invalid SMILES format. Please re-enter.{NC}")
            continue
            
        print(f"\n  {ICON_INFO} {BOLD}Molecule Analysis:{NC}")
        print(f"    • Formula           : {GREEN}{analysis['formula']}{NC}")
        print(f"    • Molecular Weight  : {analysis['mw']:.2f} g/mol")
        print(f"    • Heavy Atom Count  : {analysis['num_heavy_atoms']}")
        print(f"    • Detected Charge   : {CYAN}{analysis['charge']}{NC}")
        print(f"    • Est. Multiplicity : {CYAN}{analysis['spin']}{NC} ({'Singlet/Closed-Shell' if analysis['spin']==1 else 'Doublet/Unpaired'})")
        break

    while True:
        default_name = re.sub(r'[^a-zA-Z0-9_-]', '', analysis['formula'])
        name = input(f"\n{BOLD}Enter identifier name [Default: {default_name}]:{NC} ").strip()
        if not name:
            name = default_name
        name = re.sub(r'[^a-zA-Z0-9_-]', '_', name)
        break

    print(f"\n{BOLD}Physical / Electronic State Settings:{NC}")
    charge_in = input(f" Confirm/Override Formal Charge [{analysis['charge']}]: ").strip()
    charge = int(charge_in) if charge_in.lstrip('-').isdigit() else analysis['charge']

    spin_in = input(f" Confirm/Override Spin Multiplicity (1=Singlet, 2=Doublet, 3=Triplet) [{analysis['spin']}]: ").strip()
    spin = int(spin_in) if spin_in.isdigit() and int(spin_in) > 0 else analysis['spin']

    print(f"\n{BOLD}Conformational Search Settings:{NC}")
    n_str = input(f" Number of top conformers to extract for ORCA [Default: 10]: ").strip()
    n_conformers = int(n_str) if n_str.isdigit() else 10

    ewin_str = input(f" CREST Energy Window (kcal/mol) [Default: 10.0]: ").strip()
    ewin = float(ewin_str) if ewin_str.replace('.', '', 1).isdigit() else 10.0

    print(f"\n{BOLD}ORCA Quantum Chemistry Settings:{NC}")
    print(" Presets for Level of Theory:")
    print("   [1] B3LYP D3BJ def2-SVP OPT TightOPT FREQ (Fast & Standard)")
    print("   [2] ωB97X-D3 def2-TZVP OPT TightOPT FREQ (High Accuracy)")
    print("   [3] r2SCAN-3c OPT TightOPT FREQ (Modern Composite Method)")
    print("   [4] Custom input line")
    
    method_choice = input(" Select preset [1-4, Default: 1]: ").strip()
    if method_choice == "2":
        method = "wB97X-D3 def2-TZVP OPT TightOPT FREQ"
    elif method_choice == "3":
        method = "r2SCAN-3c OPT TightOPT FREQ"
    elif method_choice == "4":
        method = input(" Enter custom ORCA keywords line: ").strip()
        if not method:
            method = "B3LYP D3BJ def2-SVP OPT TightOPT FREQ"
    else:
        method = "B3LYP D3BJ def2-SVP OPT TightOPT FREQ"

    print(f"\n{BOLD}SLURM Cluster & Execution Settings:{NC}")
    walltime = input(" Walltime limit (HH:MM:SS) [Default: 24:00:00]: ").strip() or "24:00:00"
    ntasks = input(" CPU Cores (--ntasks) [Default: 16]: ").strip() or "16"
    partition_str = input(" SLURM Partition (leave blank if default): ").strip()
    partition = f"#SBATCH -p {partition_str}" if partition_str else "# No specific partition"
    
    maxcore_str = input(" ORCA %maxcore memory per core in MB [Default: Auto/Dynamic]: ").strip()
    maxcore = maxcore_str if maxcore_str.isdigit() else None

    return {
        "smiles": smiles,
        "mol": analysis["mol"],
        "name": name,
        "charge": charge,
        "spin": spin,
        "n_conformers": n_conformers,
        "ewin": ewin,
        "method": method,
        "walltime": walltime,
        "ntasks": ntasks,
        "partition": partition,
        "maxcore": maxcore
    }

def extract_and_generate(sdf_file: Path, pass_output: Path, cfg: dict) -> tuple[int, Path]:
    """Extracts prioritized conformers and writes XYZ, ORCA input, and SLURM submission scripts."""
    with open(pass_output, 'r') as f:
        content = f.read()
        
    match = re.search(r'\[([\d,\s]+)\]', content)
    if not match:
        print(f"{RED}ERROR: CONFPASS priority array not found in output.{NC}")
        sys.exit(1)
        
    plist = [int(x.strip()) for x in match.group(1).split(',')]

    suppl = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
    mols = [m for m in suppl if m is not None]
    n_ext = min(cfg["n_conformers"], len(plist), len(mols))

    output_dir = Path(f"{cfg['name']}_conformers")
    output_dir.mkdir(exist_ok=True)

    print(f"\n{CYAN} Writing top {n_ext} conformers to {output_dir}/{NC}")

    for i, idx in enumerate(plist[:n_ext], 1):
        if idx <= len(mols):
            mol = mols[idx-1]
            conf = mol.GetConformer()

            base_name = f"{cfg['name']}_{i}"
            xyz_path = output_dir / f"{base_name}.xyz"
            inp_path = output_dir / f"{base_name}.inp"
            sh_path = output_dir / f"{base_name}.sh"

            num_atoms = mol.GetNumAtoms()
            with open(xyz_path, 'w') as f:
                f.write(f"{num_atoms}\n{cfg['name']} - Priority Conformer {i} (CONFPASS ID {idx})\n")
                for ai in range(num_atoms):
                    pos = conf.GetAtomPosition(ai)
                    sym = mol.GetAtomWithIdx(ai).GetSymbol()
                    f.write(f"{sym:<2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n")

            with open(inp_path, 'w') as f:
                f.write(f"# {cfg['name']} - Conformer {i} (CONFPASS priority {i})\n")
                f.write(f"# Generated by CONFADS2C v2.2 on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                f.write(f"! {cfg['method']}\n\n")
                f.write(f"%pal nprocs {cfg['ntasks']} end\n")
                if cfg["maxcore"]:
                    f.write(f"%maxcore {cfg['maxcore']}\n")
                f.write(f"\n* xyzfile {cfg['charge']} {cfg['spin']} ./{base_name}.xyz\n")

            with open(sh_path, 'w') as f:
                f.write(f"""#!/bin/bash
#SBATCH -J {base_name}
#SBATCH -t {cfg['walltime']}
#SBATCH -n {cfg['ntasks']}
{cfg['partition'] if not cfg['partition'].startswith('#') else ''}

export INPUT="{base_name}.inp"
export OUTPUT="{base_name}.out"

module load orca/5.0.4 2>/dev/null || module load orca/4.2.1 2>/dev/null || true

if command -v orca &> /dev/null; then
    orca $INPUT > $OUTPUT
elif [ -n "$ORCA_BIN" ]; then
    $ORCA_BIN $INPUT > $OUTPUT
else
    echo "ERROR: ORCA executable not found in PATH." >&2
    exit 1
fi
""")
            sh_path.chmod(0o755)

            if i % 5 == 0 or i == n_ext:
                print(f"   [{i:2d}/{n_ext}] Generated {base_name}.xyz / .inp / .sh")

    return n_ext, output_dir

def main():
    total_start = time.perf_counter()

    cfg = get_user_inputs()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    work_dir = Path(f"{cfg['name']}_workflow_{timestamp}")
    work_dir.mkdir(exist_ok=True)
    os.chdir(work_dir)

    print(f"\n{BOLD}{CYAN}📁 Working Directory Created:{NC} {work_dir.resolve()}\n")

    timings = {}

    # STEP 1: SMILES -> Initial 3D Geometry
    print(f"{BOLD}{GREEN}[1/6] Native SMILES ⟶ 3D Embedding (RDKit ETKDGv3){NC}")
    step1_start = time.perf_counter()
    initial_xyz = Path(f"{cfg['name']}_initial.xyz")
    num_atoms = smiles_to_3d_xyz(cfg["mol"], initial_xyz)
    elapsed = time.perf_counter() - step1_start
    timings["3D Embedding (RDKit)"] = elapsed
    print(f" {ICON_OK} {GREEN}Generated initial 3D structure ({num_atoms} atoms) in {elapsed:.2f}s{NC}\n")

    # STEP 2: Fast Pre-Optimization via xTB
    print(f"{BOLD}{GREEN}[2/6] Geometry Optimization (xTB / GFN2){NC}")
    cmd_xtb = [
        "xtb", str(initial_xyz),
        "--opt", "--gfn2",
        "--chrg", str(cfg["charge"]),
        "--uhf", str(cfg["spin"] - 1)
    ]
    elapsed = run_command_safe(cmd_xtb, "xTB Optimization", stdout_log="xtb_opt.log", stderr_log="xtb_err.log")
    timings["xTB Optimization"] = elapsed
    
    xtb_opt_xyz = Path("xtbopt.xyz")
    if xtb_opt_xyz.exists():
        print(f" {ICON_OK} Optimized geometry saved to xtbopt.xyz\n")
    else:
        print(f" {ICON_FAIL} {RED}Error: xtbopt.xyz was not generated.{NC}")
        sys.exit(1)

    # STEP 3: Conformational Ensemble Sampling via CREST
    print(f"{BOLD}{GREEN}[3/6] Conformational Search (CREST){NC}")
    n_cores = multiprocessing.cpu_count()
    cmd_crest = [
        "crest", "xtbopt.xyz",
        "--gfn2",
        "--ewin", str(cfg["ewin"]),
        "--chrg", str(cfg["charge"]),
        "--uhf", str(cfg["spin"] - 1),
        "-T", str(n_cores)
    ]
    elapsed = run_command_safe(cmd_crest, "CREST Sampling", stdout_log="crest.log", stderr_log="crest_err.log")
    timings["CREST Search"] = elapsed

    crest_xyz = Path("crest_conformers.xyz")
    if crest_xyz.exists():
        with open(crest_xyz, 'r') as f:
            n_found = sum(1 for line in f if re.match(r'^\s*\d+\s*$', line))
        print(f" {ICON_OK} {GREEN}Discovered {n_found} unique conformers.{NC}\n")
    else:
        print(f" {ICON_FAIL} {RED}Error: crest_conformers.xyz was not generated.{NC}")
        sys.exit(1)

    # STEP 4: Format Conversion
    print(f"{BOLD}{GREEN}[4/6] Format Conversion (OpenBabel XYZ ⟶ SDF){NC}")
    sdf_out = Path(f"{cfg['name']}_conformers.sdf")
    cmd_obabel = ["obabel", "crest_conformers.xyz", "-O", str(sdf_out)]
    elapsed = run_command_safe(cmd_obabel, "XYZ to SDF Conversion")
    timings["OpenBabel Conversion"] = elapsed
    print(f" {ICON_OK} Conformers converted to {sdf_out}\n")

    # STEP 5: Prioritization & Clustering
    print(f"{BOLD}{GREEN}[5/6] Prioritization & Clustering (CONFPASS){NC}")
    confpass_log = Path(f"{cfg['name']}_confpass_output.txt")
    cmd_confpass = [sys.executable, "-m", "confpass", str(sdf_out)]
    elapsed = run_command_safe(cmd_confpass, "CONFPASS Prioritization", stdout_log=str(confpass_log))
    timings["CONFPASS Evaluation"] = elapsed
    print(f" {ICON_OK} CONFPASS output logged to {confpass_log}\n")

    # STEP 6: Top Conformer Extraction
    print(f"{BOLD}{GREEN}[6/6] Generating ORCA Inputs & SLURM Scripts{NC}")
    step6_start = time.perf_counter()
    n_ext, output_dir = extract_and_generate(sdf_out, confpass_log, cfg)
    elapsed = time.perf_counter() - step6_start
    timings["Extraction & QM Inputs"] = elapsed

    total_elapsed = time.perf_counter() - total_start
    
    summary_text = f"""======================================================================
CONFADS2C v2.2 — WORKFLOW SUMMARY REPORT
======================================================================

Molecule Identifier : {cfg['name']}
SMILES String       : {cfg['smiles']}
Execution Date      : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

CHEMICAL & PHYSICAL PROPERTIES:
  • Formula           : {rdMolDescriptors.CalcMolFormula(Chem.AddHs(cfg['mol']))}
  • Formal Charge     : {cfg['charge']}
  • Spin Multiplicity : {cfg['spin']}

CONFORMATIONAL SEARCH RESULTS:
  • CREST Conformers  : {n_found} total ensemble structures
  • Extracted Top     : {n_ext} conformers for QM re-optimization
  • Energy Window     : {cfg['ewin']} kcal/mol

ORCA & SLURM CONFIGURATION:
  • Level of Theory   : ! {cfg['method']}
  • Parallel Cores    : {cfg['ntasks']} cores (%pal nprocs {cfg['ntasks']})
  • Memory per Core   : {cfg['maxcore'] if cfg['maxcore'] else 'Dynamic default'}
  • SLURM Walltime    : {cfg['walltime']}

TIMING BREAKDOWN:
"""
    for step_name, duration in timings.items():
        percentage = (duration / total_elapsed) * 100
        summary_text += f"  • {step_name:<26} : {duration:7.2f} s ({percentage:5.1f}%)\n"
        
    summary_text += f"  --------------------------------------------------\n"
    summary_text += f"  • TOTAL WORKFLOW TIME       : {total_elapsed:7.2f} s\n\n"
    summary_text += f"CLUSTER SUBMISSION COMMANDS:\n"
    summary_text += f"  cd {output_dir.resolve()}\n"
    summary_text += f"  for f in *.sh; do sbatch \"$f\"; done\n"
    summary_text += f"======================================================================\n"

    summary_file = Path(f"{cfg['name']}_summary.txt")
    with open(summary_file, 'w') as f:
        f.write(summary_text)

    print(f"\n{BOLD}{GREEN}======================================================================{NC}")
    print(f"{BOLD}{GREEN} ✓ WORKFLOW COMPLETED SUCCESSFULLY!{NC}")
    print(f"{BOLD}{GREEN}======================================================================{NC}\n")
    print(summary_text)

if __name__ == "__main__":
    main()
