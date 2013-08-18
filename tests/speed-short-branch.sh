#!/bin/sh

dawg=$1

${dawg} - <<EOF > /dev/null
Root.Length = 100000000
Tree.Tree = (A:0.1)B;
Sim.Reps = 100
EOF
