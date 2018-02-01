#!/usr/bin/bash

python setup.py build_ext -i -DDAWG_DEBUG
# cython --embed PyDawg.cpp -o py_dawg_main.c
# clang -std=c11 -Wall -I/usr/include/python3.6m py_dawg_main.c -o py_dawg_main -L/usr/lib -lpython3
