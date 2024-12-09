#!/bin/bash

meson compile -C build

./build/app/tblite run oh.xyz --cpcm water --charge -1
