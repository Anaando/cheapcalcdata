#!/bin/bash
#SBATCH --job-name=INTP2-INTP3-Matching
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --output=INTP2-INTP3-Matching_%j.out
#SBATCH --error=INTP2-INTP3-Matching_%j.err

module purge
module load bluebear
module load bear-apps/2023a/live
module --ignore_cache load OpenMPI/4.1.5-GCC-12.3.0

cd "$SLURM_SUBMIT_DIR"

export PATH="/rds/projects/s/scanlodo-msci-2026/anaandoA/orcaprog/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH"

# Initialise Conda
source /rds/projects/s/scanlodo-msci-2026/anaandoA/miniconda3/etc/profile.d/conda.sh

# Activate Python env
conda activate ase-env

cd /rds/projects/s/scanlodo-msci-2026/anaandoA/cheapcalcdata/AtomOrdering/PTpath/M_INTP2-INTP3 || exit 1

# Run script
python AtomMatchingTemplate.py INT-P-2.xyz INT-P-3.xyz
