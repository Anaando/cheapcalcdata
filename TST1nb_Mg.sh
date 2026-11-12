#!/bin/bash
#SBATCH --job-name=TST1nb_Mg
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --output=TST1nb_Mg_%j.out
#SBATCH --error=TST1nb_Mg_%j.err

module purge
module load bluebear
module load bear-apps/2023a/live
module --ignore_cache load OpenMPI/4.1.5-GCC-12.3.0

cd "$SLURM_SUBMIT_DIR"

export PATH="/rds/projects/s/scanlodo-msci-2026/anaandoA/orcaprog/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH"

/rds/projects/s/scanlodo-msci-2026/anaandoA/orcaprog/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca TST1nb_Mg.inp > TST1nb_Mg.out
