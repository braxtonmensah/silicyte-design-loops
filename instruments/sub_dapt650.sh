#!/bin/bash
#SBATCH -J dapt650
#SBATCH -A students
#SBATCH -p h100-interactive,h100-single
#SBATCH -N 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --nice=100
#SBATCH --time=3:59:00
#SBATCH -o /N/scratch/bsmensah/silo/logs/dapt650_%j.log
#SBATCH -e /N/scratch/bsmensah/silo/logs/dapt650_%j.log
# PREREG_DAPT.md. ESM2-650M DAPT, arm B (near80, honest) and arm A (exact, model_dapt.py as written).
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
echo "dapt650 ${SLURM_JOB_ID:-manual} $(date) host $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
python -c "import torch, transformers, rapidfuzz; print('torch', torch.__version__, 'tf', transformers.__version__, 'cuda', torch.cuda.is_available())"
cd $S
# smoke: whole path on 20 steps (corpus audit, MLM, checkpoint, embed). Abort early if it breaks.
python model_dapt.py --stage pretrain --model 35m --exclude near80 --epochs 1 --max-steps 20 --save-at 1 --tag smoke35 || { echo SMOKE_FAILED; exit 1; }
# controls: plain trunks through the SAME embed code
python model_dapt.py --stage embed --model-id facebook/esm2_t12_35M_UR50D --out $S/dapt/emb/plain35.npz
python model_dapt.py --stage embed --model-id facebook/esm2_t33_650M_UR50D --out $S/dapt/emb/plain650.npz
# the two registered 650M arms (model_dapt.py recipe: lr 1e-5, 2 epochs, bs 16)
python model_dapt.py --stage pretrain --model 650m --exclude near80 --lr 1e-5 --epochs 2 --save-at 1,2 --tag m650_near80_lr1e-5
python model_dapt.py --stage pretrain --model 650m --exclude exact  --lr 1e-5 --epochs 2 --save-at 1,2 --tag m650_exact_lr1e-5
# score: plain controls first so a pipeline mismatch is visible early
mkdir -p $S/dapt/scores
E=$S/dapt/emb
python dapt_score.py --dir $S/dapt/scores --arms plain35=$E/plain35.npz,plain650=$E/plain650.npz
python dapt_score.py --dir $S/dapt/scores --arms m650_near80_ep1=$E/m650_near80_lr1e-5_ep1.npz,m650_near80_ep2=$E/m650_near80_lr1e-5_ep2.npz,m650_exact_ep1=$E/m650_exact_lr1e-5_ep1.npz,m650_exact_ep2=$E/m650_exact_lr1e-5_ep2.npz
echo "dapt650 done $(date)"
