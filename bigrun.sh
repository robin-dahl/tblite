#!/bin/bash

meson compile -C build || { echo "Compilation failed"; exit 1; }

./build/app/tblite run --method gfn1 xtbopt.xyz --cosmo water  
