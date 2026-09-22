#!/bin/bash
set -e
cd /home/admin123/OpenFOAM/admin123-13/run/MetalCombustion_Case21_BunsenFlame
source /opt/openfoam13/etc/bashrc

echo '=== Step 1: Reconstructing Eulerian Gas Fields (0.0 to 0.45 s) ==='
reconstructPar -allTimes > log.reconstructPar 2>&1
echo 'Eulerian fields reconstructed successfully!'

echo '=== Step 2: Reconstructing Lagrangian Particle Cloud ==='
reconstructPar -lagrangian -allTimes > log.reconstructPar_lagrangian 2>&1
echo 'Lagrangian cloud reconstructed successfully!'

echo '=== Step 3: Creating ParaView Dummy File ==='
touch case.foam

echo '=== Step 4: Compressing Reconstructed Case for Laptop Transfer ==='
# Exclude processor* directories to keep the archive lightweight
tar -czvf MetalCombustion_Case21_Reconstructed.tar.gz     --exclude='processor*'     --exclude='*.o'     0* constant system case.foam *.sh log.*

echo '=== Step 5: Copying to Windows Desktop for Easy Access ==='
cp MetalCombustion_Case21_Reconstructed.tar.gz /mnt/c/Users/Admin/Desktop/ 2>/dev/null || true

echo '=== ALL DONE! Archive is ready for transfer ==='
ls -lh MetalCombustion_Case21_Reconstructed.tar.gz
