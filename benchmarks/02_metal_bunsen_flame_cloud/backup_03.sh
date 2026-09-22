#!/bin/bash
set -e
cd /home/admin123/OpenFOAM/admin123-13/run/MetalCombustion_Case21_BunsenFlame
echo "Creating backup of 0.3..."
for p in processor*; do
    if [ -d "$p/0.3" ] && [ ! -d "$p/0.3_orig" ]; then
        cp -r "$p/0.3" "$p/0.3_orig"
    fi
done
echo "Backup successfully created!"
