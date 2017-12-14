#!/usr/bin/bash

clang++ -std=c++11 -Wall wood_builder.cpp -o wood_builder -lboost_system

# example newick parse for tree: ((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);
# std::vector of length 5, capacity 8 = {{
#     label = "Man",
#     length = 0.100000001,
#     anc = 2,
#     right = 0
#   }, {
#     label = "Monkey",
#     length = 0.100000001,
#     anc = 1,
#     right = 0
#   }, {
#
#     label = "",
#     length = 0.200000003,
#     anc = 2,
#     right = 2
#   }, {
#     label = "Dawg",
#     length = 0.25,
#     anc = 1,
#     right = 0
#   }, {
#     label = "",
#     length = 0,
#     anc = 0,
#     right = 2
#   }}
