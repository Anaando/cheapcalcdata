#!/bin/bash
#SBATCH --job-name=cat1chp
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=12:00:00
#SBATCH --output=cat1chp_%j.out
#SBATCH --error=cat1chp_%j.err

module purge
module load bluebear
module load bear-apps/2023a/live
module --ignore_cache load OpenMPI/4.1.5-GCC-12.3.0

cd "$SLURM_SUBMIT_DIR"

export PATH="/rds/projects/s/scanlodo-msci-2026/anaandoA/orcaprog/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH"

/rds/projects/s/scanlodo-msci-2026/anaandoA/orcaprog/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca TS-P-1_ts_f4s11.inp >  TS-P-1_ts_f4s11.out
