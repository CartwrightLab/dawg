#!/usr/bin/bash

clang++ -std=c++11 -Wall -L ./ -lPyDawg.cpython-36m-x86_64-linux-gnu.so \
		dawg_park.cpp -o dawg_park