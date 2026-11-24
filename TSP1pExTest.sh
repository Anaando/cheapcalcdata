#!/bin/bash
#SBATCH --job-name=TSP1pExTest
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=9-00:00:00
#SBATCH --output=figS11_TSP1pExTest_%j.out
#SBATCH --error=figS11_TSP1pExTest_%j.err

module purge
module load bluebear
module load bear-apps/2023a/live
module --ignore_cache load OpenMPI/4.1.5-GCC-12.3.0

cd "$SLURM_SUBMIT_DIR"

export PATH="/rds/projects/s/scanlodo-msci-2026/anaandoA/orcaprog/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH"

/rds/projects/s/scanlodo-msci-2026/anaandoA/orcaprog/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca TSP1pExTest.inp > TSP1pExTest.out
