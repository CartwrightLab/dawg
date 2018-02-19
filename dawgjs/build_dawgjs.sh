#!/usr/bin/bash

hash em++ &> /dev/null || { echo 'Emscripten toolchain not found.'; }

# -I ../../boost_1_66_0 \
# -I ../src/include \
# -I ../build/release/src/include/ \
# ../src/lib/output.cpp \
# ../src/lib/mutt.cpp ../src/lib/models.cpp \
# ../src/lib/ma.cpp ../src/lib/matic.cpp \
# ../src/dawg_node_program.cpp \


# Run emscripten commands to generate a DAWG js file
# Run CMake first to get 'config.h' in build/release/src/include/
# Trick parse functions have 'unresolve symbol' errors from EMCC compiler (spirit?)
#   Fix is to not compile parse.cpp
em++ -std=c++11 -Wall -fPIC --minify 0 -O2 --bind \
    -s DEMANGLE_SUPPORT=1 \
    -s NO_EXIT_RUNTIME=1 \
    -I ~/code/boost_1_66_0 \
    -I ../src/include \
    ../src/lib/mutt.cpp \
    dawg_node_api.cpp -o dawg.js
