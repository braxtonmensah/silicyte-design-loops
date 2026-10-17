#!/bin/bash
#SBATCH -J maxoutFT
#SBATCH -A students
#SBATCH -p h100-single,h100-interactive
#SBATCH -N 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=3:59:00
#SBATCH -o /N/scratch/bsmensah/silo/logs/maxoutft_%j.log
#SBATCH -e /N/scratch/bsmensah/silo/logs/maxoutft_%j.log
# MAXOUT rung 4: fine-tune on HONEST folds, seeds averaged, plus the null that was owed.
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
echo "maxoutFT ${SLURM_JOB_ID:-manual} $(date) host $(hostname)"
nvidia-smi --query-gpu=name --format=csv,noheader || true
cd $S
# embed once with the incumbent so the frozen arm is paired on identical folds
python maxout_embed.py --folds $S/maxout_folds.json --out $S/emb --models esm2_35m,esm2_650m
for T in mic hc; do
  python maxout_finetune.py --folds $S/maxout_folds.json --target $T \
    --model facebook/esm2_t33_650M_UR50D --seeds 3 --epochs 4 --bs 16 \
    --emb $S/emb/esm2_650m_mean.npz --out $S/maxout
done
# the null arm, one seed is enough to see whether a leak remains
python maxout_finetune.py --folds $S/maxout_folds.json --target mic \
  --model facebook/esm2_t33_650M_UR50D --seeds 1 --epochs 4 --bs 16 \
  --shuffle-labels --out $S/maxout
echo "maxoutFT done $(date)"
