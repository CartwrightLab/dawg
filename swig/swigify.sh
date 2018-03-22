#!/usr/bin/bash

hash swig 2>/dev/null || { echo >&2 "Swig is not installed."; }

ls | grep "_dawg.so" && { echo "Removing _dawg.so"; rm _dawg.so; }

# First command generates a *dawg.py* and *dawg_wrap.cxx* file
swig -python -c++ -o dawg_wrap.cpp dawg.i
# Compile the dawg lib and SWIG-generated file
clang++ -std=c++11 -Wall -O0 -I /usr/include/python3.6m/ \
    -I ../src/include -c -fPIC dawg.cpp dawg_wrap.cpp ../src/lib/ma.cpp \
    ../src/lib/matic.cpp ../src/lib/models.cpp ../src/lib/mutt.cpp \
    ../src/lib/output.cpp ../src/lib/parse.cpp
clang++ -shared -Wl,-soname,"_dawg.so" -o _dawg.so \
    dawg.o dawg_wrap.o ma.o matic.o models.o mutt.o output.o parse.o -lpython3
