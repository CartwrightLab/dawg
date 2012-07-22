#!/bin/sh

dawg=$1

${dawg} - <<EOF > /dev/null
Root.Length = 1000000
Tree.Tree = (A:1)B;
Sim.Reps = 300
Indel.Model.Del = power
Indel.Rate.Del = 0.1
Indel.Params.Del = 1.5, 10000
EOF
