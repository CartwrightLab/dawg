#!/usr/bin/bash

# A very simple and quick way to build the project
clang++ -std=c++11 -Wall -O0 -g \
    -I src/include \ # the location of the dawg headers
    -c src/lib/matic.cpp \
    -o dawg -llibboost_program_options -llibboost_system
