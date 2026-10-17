#!/bin/bash
#SBATCH -J maxoutBASE
#SBATCH -A students
#SBATCH -p h100-single,h100-interactive
#SBATCH -N 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=3:59:00
#SBATCH -o /N/scratch/bsmensah/silo/logs/maxoutbase_%j.log
#SBATCH -e /N/scratch/bsmensah/silo/logs/maxoutbase_%j.log
# MAXOUT rung 6: which BASE MODEL wins once the backbone is allowed to train?
# The frozen ladder said ESM2 saturates (35M -> 650M buys +0.004, 3B is worse). That is a statement
# about a linear probe, not about the family. Fine-tuning is where a better base model can show.
set -uo pipefail
S=/N/scratch/bsmensah/silo
mkdir -p $S/logs $S/maxout
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
echo "maxoutBASE ${SLURM_JOB_ID:-manual} $(date) host $(hostname)"
nvidia-smi --query-gpu=name --format=csv,noheader || true
cd $S
# one seed each: this rung ranks FAMILIES, and the winner gets the 3-seed average afterwards
for M in facebook/esm2_t12_35M_UR50D facebook/esm1b_t33_650M_UR50S ElnaggarLab/ankh-base; do
  echo "######## $M"
  python maxout_finetune.py --folds $S/maxout_folds.json --target mic \
    --model "$M" --seeds 1 --epochs 4 --bs 16 --out $S/maxout || echo "ARM FAILED $M"
done
echo "maxoutBASE done $(date)"
