#!/bin/bash
#SBATCH -J dapt35
#SBATCH -A students
#SBATCH -p h100-interactive,h100-single
#SBATCH -N 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --nice=100
#SBATCH --time=3:59:00
#SBATCH -o /N/scratch/bsmensah/silo/logs/dapt35_%j.log
#SBATCH -e /N/scratch/bsmensah/silo/logs/dapt35_%j.log
# PREREG_DAPT.md. ESM2-35M DAPT dose grid: {near80, exact} x lr {1e-5, 1e-4} x checkpoints {1,2,4,8}.
set -uo pipefail
S=/N/scratch/bsmensah/silo
mkdir -p $S/logs $S/dapt/emb $S/out
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
echo "dapt35 ${SLURM_JOB_ID:-manual} $(date) host $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
cd $S
mkdir -p $S/dapt/scores
E=$S/dapt/emb
python model_dapt.py --stage embed --model-id facebook/esm2_t12_35M_UR50D --out $E/plain35_j35.npz
for EX in near80 exact; do
  for LR in 1e-5 1e-4; do
    python model_dapt.py --stage pretrain --model 35m --exclude $EX --lr $LR --epochs 8 --save-at 1,2,4,8 --tag m35_${EX}_lr${LR}
  done
done
ARMS="plain35_j35=$E/plain35_j35.npz"
for EX in near80 exact; do for LR in 1e-5 1e-4; do for EP in 1 2 4 8; do
  ARMS="$ARMS,m35_${EX}_lr${LR}_ep${EP}=$E/m35_${EX}_lr${LR}_ep${EP}.npz"
done; done; done
python dapt_score.py --dir $S/dapt/scores --arms "$ARMS"
echo "dapt35 done $(date)"
