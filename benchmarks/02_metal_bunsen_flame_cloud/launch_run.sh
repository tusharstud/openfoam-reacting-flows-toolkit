#!/bin/bash
cd /home/admin123/OpenFOAM/admin123-13/run/MetalCombustion_Case21_BunsenFlame
tmux kill-session -t openfoam_run 2>/dev/null || true
tmux new-session -d -s openfoam_run "cd /home/admin123/OpenFOAM/admin123-13/run/MetalCombustion_Case21_BunsenFlame && source /opt/openfoam13/etc/bashrc && mpirun -np 60 foamRun -parallel > log.foamRun_ignited 2>&1"
echo "foamRun started in tmux session 'openfoam_run' logging to log.foamRun_ignited"