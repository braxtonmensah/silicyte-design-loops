#!/bin/bash
#SBATCH -J basesweep
#SBATCH -A students
#SBATCH -p h100-interactive,h100-single
#SBATCH -N 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=2:30:00
#SBATCH -o /N/scratch/bsmensah/silo/logs/basesweep_%j.log
#SBATCH -e /N/scratch/bsmensah/silo/logs/basesweep_%j.log
# BASE SWEEP: keep all hidden layers of many base models, then score them on the honest ruler.
# Usage: sbatch sub_sweep.sh "<shell command to run inside $S/sweep>"
set -uo pipefail
S=/N/scratch/bsmensah/silo
mkdir -p $S/logs $S/sweep $S/emb_sweep
if ! command -v module >/dev/null 2>&1; then
  for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /etc/profile.d/lmod.sh; do
    test ! -f "$f" || { source "$f"; break; }
  done
fi
set +u
module load python/3.13.5 2>/dev/null || true
export LD_LIBRARY_PATH="/N/soft/rhel8/python/gnu/3.13.5/lib:/N/soft/rhel8/gcc/14.2.0/lib64:${LD_LIBRARY_PATH:-}"
export LD_PRELOAD="/N/soft/rhel8/gcc/14.2.0/lib64/libstdc++.so.6${LD_PRELOAD:+:$LD_PRELOAD}"
source /N/slate/bsmensah/arc26/lightning/venv/bin/activate
set -u
V=/N/slate/bsmensah/arc26/lightning/venv/lib/python3*/site-packages
find $V/nvidia $V/torch/lib -name '*.so*' -size +2M -print0 2>/dev/null | xargs -0 -P 6 -n 1 cat > /dev/null || true
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=8 TOKENIZERS_PARALLELISM=false
export HF_HOME=$S/hf HF_HUB_OFFLINE=0
# extra site-packages (scipy etc.) installed with --target so the shared venv is never modified
export PYTHONPATH=$S/pylib${PYTHONPATH:+:$PYTHONPATH}
echo "basesweep ${SLURM_JOB_ID:-manual} $(date) host $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
cd $S/sweep
echo "CMD: $1"
bash -c "$1"
echo "basesweep done $(date)"
