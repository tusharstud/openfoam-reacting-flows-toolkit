#!/bin/bash
# =============================================================================
# WORKSTATION RESULTS RECONSTRUCTION & PACKAGING SCRIPT (PILLAR 2)
# =============================================================================
set -e

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$CASE_DIR"

echo "============================================================================="
echo "   RECONSTRUCTING & PACKAGING PILLAR 2 RESULTS ON WORKSTATION"
echo "   Working Directory: $CASE_DIR"
echo "============================================================================="

# 1. Source OpenFOAM 13
if [ -f /opt/openfoam13/etc/bashrc ]; then
    source /opt/openfoam13/etc/bashrc
elif [ -f "$HOME/OpenFOAM/OpenFOAM-13/etc/bashrc" ]; then
    source "$HOME/OpenFOAM/OpenFOAM-13/etc/bashrc"
else
    echo "[!] Warning: OpenFOAM environment assuming already sourced."
fi

# 2. Check if simulation is still running or finished
echo "[1/4] Checking simulation process status..."
if pgrep -f "foamRun" > /dev/null; then
    echo "      [!] Notice: foamRun process is currently active."
    echo "      Tail of log file:"
    tail -n 10 log.foamRun_stabilized || true
else
    echo "      [OK] foamRun simulation run has completed!"
fi

# 3. Reconstruct the latest stabilized time solution from 60 processors
echo "[2/4] Reconstructing latest time solution from 60 cores..."
reconstructPar -latestTime > log.reconstructPar 2>&1
echo "      reconstructPar finished with code: $?"

LATEST_TIME=$(foamListTimes -latestTime 2>/dev/null || true)
if [ -z "$LATEST_TIME" ]; then
    LATEST_TIME=$(ls -d 0.[0-9]* [1-9]* 2>/dev/null | grep -v "processor" | sort -V | tail -n 1)
fi
echo "      Latest reconstructed time directory: '$LATEST_TIME'"

# 4. Check if we also want to reconstruct recent time steps for particle tracking
echo "[3/4] Reconstructing recent time steps for particle streamline tracking..."
foamListTimes 2>/dev/null | tail -n 5 | while read t; do
    if [ ! -d "$t" ] && [ -n "$t" ]; then
        echo "      Reconstructing time: $t"
        reconstructPar -time "$t" >> log.reconstructPar 2>&1 || true
    fi
done

# 5. Create lightweight tarball archive (~20-40 MB)
ARCHIVE_NAME="Pillar2_Workstation_Completed.tar.gz"
echo "[4/4] Creating lightweight transfer archive: $ARCHIVE_NAME..."

# Gather list of reconstructed time directories (excluding processors)
RECON_TIMES=$(ls -d 0.[0-9]* [1-9]* 2>/dev/null | grep -v "processor" | sort -V | tail -n 10)

tar -czvf "$ARCHIVE_NAME"     $RECON_TIMES     constant     system     log.foamRun_stabilized     log.reconstructPar     run_workstation_stabilized.sh

ARCHIVE_SIZE=$(du -h "$ARCHIVE_NAME" | cut -f1)
echo "============================================================================="
echo "   [SUCCESS] PACKAGE CREATED: $CASE_DIR/$ARCHIVE_NAME ($ARCHIVE_SIZE)"
echo "   "
echo "   HOW TO TRANSFER TO LAPTOP:"
echo "   Windows Path: C:\Users\Admin\Desktop\MetalCombustion_Case21_BunsenFlame\$ARCHIVE_NAME"
echo "   1. Copy this single file ($ARCHIVE_NAME) to a USB Drive, OneDrive, or Google Drive."
echo "   2. Paste it on your laptop in: MetalCombustion_Case21_BunsenFlame/"
echo "   3. Extract: tar -xzvf $ARCHIVE_NAME"
echo "============================================================================="
