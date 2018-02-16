#!/usr/bin/bash

hash em++ &> /dev/null || { echo 'Emscripten toolchain not found.'; }

# Run emscripten commands to generate a DAWG js file
# Run CMake first to get 'config.h' in build/release/src/include/
em++ -std=c++11 -O0 -Wall -fPIC \
    -I ../../boost_1_66_0 \
    -I ../src/include \
    -I build/release/src/include/ \
    ../src/lib/output.cpp ../src/lib/parse.cpp \
    ../src/lib/mutt.cpp ../src/lib/models.cpp \
    ../src/lib/ma.cpp ../src/lib/matic.cpp \
    -o dawg.js
