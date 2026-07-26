#!/bin/bash

meson compile -C build || { echo "Compilation failed"; exit 1; }

mol=ibu.xyz

./build/app/tblite run $mol --ddcosmo water --solv-dipoles --solv-quadrupoles 
./build/app/tblite run $mol --cosmo water --solv-dipoles --solv-quadrupoles 
./build/app/tblite run $mol --cosmo water --solv-full-density 
 
