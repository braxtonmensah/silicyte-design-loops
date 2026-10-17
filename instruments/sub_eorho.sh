#!/bin/bash
#SBATCH -J siloEORHO
#SBATCH -A students
#SBATCH -p h100-single,h100-interactive
#SBATCH -N 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3:59:00
#SBATCH -o /N/scratch/bsmensah/silo/logs/eorho_%j.log
#SBATCH -e /N/scratch/bsmensah/silo/logs/eorho_%j.log
# SILO-EO-RHO: all 12 arms in ONE allocation, SEQUENTIALLY.
#
# Sequential is not a convenience here. Concurrent arms sharing one directory truncated the ESM
# embedding cache mid-read and killed job 12243747 with EOFError; the caches came back at 53% and
# 49% of size. Each arm below also gets its OWN --out, so no two arms touch one cache file even if
# this script is ever run twice at once.
#
# One allocation also means every arm sees the same node, the same CUDA libraries and the same
# ESM2 weights, which removes a per-job confound from a comparison designed around a small effect.
#
# Registered in PREREG_SILO_EORHO.md. 2 rungs x 3 seeds x {real, shuffled} = 12 arms, ~6 min each.
set -uo pipefail
S=/N/scratch/bsmensah/silo
mkdir -p $S/logs $S/out_eorho
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
# Warm the big shared libraries before the first CUDA call: job 12243503 sat 7 minutes at 7
# major-faults/sec inside cuDNN plan search off a slow scratch stripe.
V=/N/slate/bsmensah/arc26/lightning/venv/lib/python3*/site-packages
find $V/nvidia $V/torch/lib -name '*.so*' -size +2M -print0 2>/dev/null \
  | xargs -0 -P 6 -n 1 cat > /dev/null || true
echo "libs warmed $(date)"
export PYTHONPATH=$S/pylibs:$S/amp/src:$S
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=8 HF_HOME=$S/hf HF_HUB_OFFLINE=1
echo "eorho pack ${SLURM_JOB_ID:-manual} $(date) host $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

ok=0; fail=0
for FRAC in 1.00 0.20; do
  for SEED in 101 202 303; do
    for ARM in real shuf; do
      EXTRA=""
      [ "$ARM" = "shuf" ] && EXTRA="--shuffled"
      TAG="f${FRAC}_s${SEED}_${ARM}"
      OUT="$S/out_eorho/$TAG"
      # a requeue must not redo finished work, and must not half-overwrite a finished arm either
      if ls "$OUT"/silo2_*.json >/dev/null 2>&1; then
        echo "=== SKIP $TAG (already has json)"; ok=$((ok+1)); continue
      fi
      mkdir -p "$OUT"
      echo "=== ARM $TAG  $(date +%H:%M:%S)"
      if python $S/silo_eorho.py --amp-repo $S/amp --frac "$FRAC" --seed "$SEED" \
           --gens 8 --lib 20000 --k 200 --epochs 2 --xclass --only-class2 --ad-gate 0.40 \
           --lib-meter --tag "_$TAG" --out "$OUT" $EXTRA; then
        ok=$((ok+1))
      else
        fail=$((fail+1)); echo "=== ARM $TAG FAILED"
      fi
    done
  done
done
echo "eorho pack done $(date)  ok=$ok fail=$fail of 12"
