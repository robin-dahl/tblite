#!/bin/bash

meson compile -C build || { echo "Compilation failed"; exit 1; }

./build/app/tblite run xtbopt.xyz --cpcm water 
