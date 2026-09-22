#!/bin/bash
cd /home/admin123/OpenFOAM/admin123-13/run/MetalCombustion_Case21_BunsenFlame
source /opt/openfoam13/etc/bashrc
mpirun -np 60 foamRun -parallel > test_run.log 2>&1
echo "foamRun finished with code $?"