#!/bin/bash

meson compile -C build || { echo "Compilation failed"; exit 1; }

#mol=coord.xyz
#mol=acetonitrile.xyz
mol=ibu.xyz
#mol=cl.xyz

./build/app/tblite run $mol --cpcm 80.000 --solv-monopoles 
./build/app/tblite run $mol --cpcm   80.000 --solv-monopoles --solv-dipoles  
./build/app/tblite run $mol --cpcm   80.000 --solv-dipoles --solv-quadrupoles 
./build/app/tblite run $mol --cpcm   80.000 --solv-monopoles --solv-dipoles --solv-quadrupoles 
./build/app/tblite run $mol --cpcm   80.000 --solv-full-density 
 
