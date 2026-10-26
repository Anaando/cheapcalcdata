#!/bin/bash
#SBATCH --job-name=INTP3-INTP4_chpNEB
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --output=INTP3-INTP4_chpNEB_%j.out
#SBATCH --error=INTP3-INTP4_chpNEB_%j.err

module purge
module load bluebear
module load bear-apps/2023a/live
module --ignore_cache load OpenMPI/4.1.5-GCC-12.3.0

cd "$SLURM_SUBMIT_DIR"

export PATH="/rds/projects/s/scanlodo-msci-2026/anaandoA/orcaprog/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH"

/rds/projects/s/scanlodo-msci-2026/anaandoA/orcaprog/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca INTP3-INTP4_chpNEB.inp > INTP3-INTP4_chpNEB.out
