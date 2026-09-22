#!/bin/bash
set -e
cd /home/admin123/OpenFOAM/admin123-13/run/MetalCombustion_Case21_BunsenFlame

echo "1. Killing any foamRun processes..."
killall -9 foamRun 2>/dev/null || true

echo "2. Cleaning any micro-step directories..."
rm -rf processor*/0.300*

echo "3. Restoring clean 0.3/T from 0.3_orig/T across all 60 processors..."
for p in processor*; do
    if [ -d "$p/0.3_orig" ]; then
        cp "$p/0.3_orig/T" "$p/0.3/T"
    fi
done

echo "4. Sourcing OpenFOAM 13 environment..."
source /opt/openfoam13/etc/bashrc

echo "5. Re-applying spark with setFields -latestTime -parallel..."
mpirun -np 60 setFields -latestTime -parallel > log.setFields 2>&1
echo "setFields exited with code: $?"

echo "6. Verifying T distribution in 0.3..."
python3 check_T_dist.py

echo "Preparation and spark complete!"