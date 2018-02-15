#!/usr/bin/bash

# $1 times of iteration
function timerun {
	for (( i = 0; i < $2; i++ ))
	do
		./"$1" &> /dev/null
	done
}

echo "Timing Pydawg - Akita version"
time timerun Akita.py 1000
echo "Timing system dawg (develop)"
time timerun "dawg akita.dawg -o /tmp/akita_pure_cpp.fasta" 1000
