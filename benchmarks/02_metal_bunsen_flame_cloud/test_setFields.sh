#!/bin/bash
source /opt/openfoam13/etc/bashrc
cd /home/admin123/OpenFOAM/admin123-13/run/MetalCombustion_Case21_BunsenFlame
echo "Testing setFields -parallel -latestTime..."
mpirun -np 60 setFields -parallel -latestTime > log.setFields 2>&1
echo "Exit code: $?"
cat log.setFields
