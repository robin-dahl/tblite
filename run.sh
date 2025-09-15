#!/bin/bash

meson compile -C build || { echo "Compilation failed"; exit 1; }

./build/app/tblite --method gfn2 ibu.xyz --cpcm water  
