#!/bin/bash
#SBATCH -J daptft
#SBATCH -A students
#SBATCH -p h100-interactive,h100-single
#SBATCH -N 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --nice=100
#SBATCH --time=3:59:00
#SBATCH -o /N/scratch/bsmensah/silo/logs/daptft_%j.log
#SBATCH -e /N/scratch/bsmensah/silo/logs/daptft_%j.log
# usage: sbatch sub_dapt_ft.sh <checkpoint dir> <tag>
# Fine-tune an ADAPTED trunk on the honest folds with the registered maxout_finetune.py protocol
# (4 epochs, 3 seeds averaged, bs 16, lr 1e-5/1e-3) plus the shuffled-label null for both targets.
set -uo pipefail
S=/N/scratch/bsmensah/silo
CK=$1; TAG=$2
mkdir -p $S/logs $S/dapt/ft/$TAG
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
echo "daptft $TAG ${SLURM_JOB_ID:-manual} $(date) host $(hostname) ckpt $CK"
cd $S
for T in hc mic; do
  python maxout_finetune.py --folds $S/maxout_folds.json --target $T --model $CK --seeds 3 --epochs 4 --bs 16 --out $S/dapt/ft/$TAG
done
for T in hc mic; do
  python maxout_finetune.py --folds $S/maxout_folds.json --target $T --model $CK --seeds 1 --epochs 4 --bs 16 --shuffle-labels --out $S/dapt/ft/$TAG
done
echo "daptft $TAG done $(date)"
