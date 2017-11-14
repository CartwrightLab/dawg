#!/usr/bin/bash

swig -python example.i
clang -c -fpic example.c example_wrap.c -I /usr/include/python3.6m/
ld -shared example.o example_wrap.o -o _example.so
