#!/bin/bash
#SBATCH -J dapt35L
#SBATCH -A students
#SBATCH -p h100-interactive,h100-single
#SBATCH -N 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --nice=100
#SBATCH --time=3:59:00
#SBATCH -o /N/scratch/bsmensah/silo/logs/dapt35L_%j.log
#SBATCH -e /N/scratch/bsmensah/silo/logs/dapt35L_%j.log
# PREREG_DAPT.md addendum 1: 35M, arm B, lr 1e-4, 32 epochs, checkpoints 16 and 32.
set -uo pipefail
S=/N/scratch/bsmensah/silo
mkdir -p $S/logs $S/dapt/emb
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
export HF_HOME=$S/hf HF_HUB_OFFLINE=1
export PYTHONPATH=$S/pylibs:$S/amp/src:$S
echo "dapt35L ${SLURM_JOB_ID:-manual} $(date) host $(hostname)"
cd $S
mkdir -p $S/dapt/scores
E=$S/dapt/emb
python model_dapt.py --stage pretrain --model 35m --exclude near80 --lr 1e-4 --epochs 32 --save-at 16,32 --tag m35_near80_lr1e-4_long
python dapt_score.py --dir $S/dapt/scores --arms m35_near80_lr1e-4_long_ep16=$E/m35_near80_lr1e-4_long_ep16.npz,m35_near80_lr1e-4_long_ep32=$E/m35_near80_lr1e-4_long_ep32.npz
echo "dapt35L done $(date)"
