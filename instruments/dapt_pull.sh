#!/bin/bash
# pull small result files (scores json, OOF npz, pretrain logs, multitask json) from Quartz
cd /c/Users/bmens/NQ_local/silicyte-loop || exit 1
mkdir -p dapt_out dapt_out/logs dapt_out/mt
timeout 300 wsl ssh quartz "cd /N/scratch/bsmensah/silo/dapt/scores && tar cf - score_*.json oof_*.npz 2>/dev/null" | tar xf - -C dapt_out
for d in $(timeout 60 wsl ssh quartz "ls /N/scratch/bsmensah/silo/dapt | grep -E '^m(35|650)_'"); do
  timeout 60 wsl ssh quartz "cat /N/scratch/bsmensah/silo/dapt/$d/log.json" > dapt_out/logs/$d.json 2>/dev/null
done
timeout 120 wsl ssh quartz "cd /N/scratch/bsmensah/silo/dapt/mt && tar cf - mt_*.json mt_*_oof.npz 2>/dev/null" | tar xf - -C dapt_out/mt
ls dapt_out | head -50
