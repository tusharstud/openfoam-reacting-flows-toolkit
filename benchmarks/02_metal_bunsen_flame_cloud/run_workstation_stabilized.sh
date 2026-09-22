#!/bin/bash
# =============================================================================
# AUTOMATED WORKSTATION LAUNCH SCRIPT FOR STABILIZED BUNSEN FLAME (CASE 21)
# =============================================================================
set -e

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -d "$CASE_DIR/MetalCombustion_Case21_BunsenFlame" ]; then
    CASE_DIR="$CASE_DIR/MetalCombustion_Case21_BunsenFlame"
fi
cd "$CASE_DIR"

echo "============================================================================="
echo "   STARTING STABILIZED AL BUNSEN FLAME RUN ON WORKSTATION (60 CORES)"
echo "   Working Directory: $CASE_DIR"
echo "============================================================================="

# 1. Source OpenFOAM 13 environment
if [ -f /opt/openfoam13/etc/bashrc ]; then
    source /opt/openfoam13/etc/bashrc
elif [ -f "$HOME/OpenFOAM/OpenFOAM-13/etc/bashrc" ]; then
    source "$HOME/OpenFOAM/OpenFOAM-13/etc/bashrc"
else
    echo "[!] Warning: /opt/openfoam13/etc/bashrc not found. Assuming environment is already sourced."
fi

# 2. Stop any existing foamRun processes
echo "[1/6] Stopping any old foamRun processes..."
killall -9 foamRun 2>/dev/null || true
sleep 1

# 3. Clean up post-0.3 time directories from previous blow-off run
echo "[2/6] Cleaning time directories > 0.300 from previous runs..."
rm -rf 0.30* 0.31* 0.32* 0.33* 0.34* 0.35* 0.36* 0.37* 0.38* 0.39* 0.4*
rm -rf processor*/0.30* processor*/0.31* processor*/0.32* processor*/0.33* processor*/0.34* processor*/0.35* processor*/0.36* processor*/0.37* processor*/0.38* processor*/0.39* processor*/0.4*

# 4. Decompose the newly updated 0.3 fields (with 1800K pilot and parabolic U)
echo "[3/6] Decomposing updated 0.3 fields across 60 processors..."
decomposePar -fields -time 0.3 > log.decomposePar 2>&1
echo "      decomposePar completed with code: $?"

# 5. Apply localized ignition spark in parallel
echo "[4/6] Applying localized spark at nozzle lip (r ~ 10mm) in parallel..."
mpirun -np 60 setFields -latestTime -parallel > log.setFields 2>&1
echo "      setFields completed with code: $?"

# 6. Launch parallel simulation in background / tmux
LOG_FILE="log.foamRun_stabilized"
echo "[5/6] Launching 60-core parallel simulation logging to $LOG_FILE..."

if command -v tmux >/dev/null 2>&1; then
    tmux kill-session -t bunsen_run 2>/dev/null || true
    tmux new-session -d -s bunsen_run "cd '$CASE_DIR' && mpirun -np 60 foamRun -parallel > '$LOG_FILE' 2>&1"
    echo "      Started in tmux session 'bunsen_run'."
    echo "      To attach and view: tmux attach -t bunsen_run"
else
    nohup mpirun -np 60 foamRun -parallel > "$LOG_FILE" 2>&1 &
    echo "      Started with nohup (PID: $!)."
fi

echo "[6/6] Simulation is RUNNING!"
echo "============================================================================="
echo "   LIVE MONITORING COMMANDS:"
echo "   1. View progress live:         tail -f $LOG_FILE | grep -E 'Time =|Courant'"
echo "   2. Monitor entire log:         tail -f $LOG_FILE"
echo "============================================================================="
