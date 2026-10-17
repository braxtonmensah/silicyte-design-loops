#!/bin/bash
#SBATCH -J daptmt
#SBATCH -A students
#SBATCH -p h100-interactive,h100-single
#SBATCH -N 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --nice=100
#SBATCH --time=3:59:00
#SBATCH -o /N/scratch/bsmensah/silo/logs/daptmt_%j.log
#SBATCH -e /N/scratch/bsmensah/silo/logs/daptmt_%j.log
# usage: sbatch sub_dapt_mt.sh <backbone path or HF id> <tag>
set -uo pipefail
S=/N/scratch/bsmensah/silo
BB=$1; TAG=$2
mkdir -p $S/logs $S/dapt/mt
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
echo "daptmt $TAG ${SLURM_JOB_ID:-manual} $(date) host $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
cd $S
# smoke: one fold, one seed, no nulls, tiny
python dapt_multitask.py --backbone $BB --tag smoke_$TAG --out $S/dapt/mt --seeds 20261014 --max-folds 1 --nulls 0 --targets hc --epochs 1 || { echo SMOKE_FAILED; exit 1; }
python dapt_multitask.py --backbone $BB --tag $TAG --out $S/dapt/mt
echo "daptmt $TAG done $(date)"
