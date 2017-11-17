#!/usr/bin/bash

hash swig 2>/dev/null || { echo >&2 "Swig is not installed."; }

# First command generates a *dawg.py* and *dawg_interface_wrap.c* file
swig -python dawg_interface.i
# Compile the dawg lib and SWIG-generated file
clang++ -std=c++11 -Wall -O0 -I /usr/include/python3.6m/ \
    -I ../src/include -c -fPIC dawg.cpp dawg_interface_wrap.c ../src/lib/ma.cpp \
    ../src/lib/matic.cpp ../src/lib/models.cpp ../src/lib/mutt.cpp \
    ../src/lib/output.cpp ../src/lib/parse.cpp
clang++ -shared -Wl,-soname,"dawg.so" -o dawg.so \
    dawg.o dawg_interface_wrap.o ma.o matic.o models.o mutt.o output.o parse.o -lpython3
