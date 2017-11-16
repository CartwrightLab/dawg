#!/usr/bin/bash

clang++ -std=c++11 -Wall -O0 -I /usr/include/python3.6m/ \
    -I ../src/include -c -fPIC dawg.cpp ../src/lib/ma.cpp \
    ../src/lib/matic.cpp ../src/lib/models.cpp ../src/lib/mutt.cpp \
    ../src/lib/output.cpp ../src/lib/parse.cpp
clang++ -shared -Wl,-soname,"dawg.so" -o dawg.so dawg.o ma.o matic.o models.o mutt.o output.o parse.o -lpython3 -lboost_python3
