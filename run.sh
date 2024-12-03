#!/bin/bash

meson compile -C build

./build/app/tblite run coord.xyz --cpcm water