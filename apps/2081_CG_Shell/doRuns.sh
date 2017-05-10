#!/bin/bash

EXEC="$(pwd)/../../cmake-build-release/apps/shellProject"
MESH="thinPlate_"
OUTBASE="output_"


for p in $(echo "1 2 3 4"); do

    for N in $(echo "4 10 20 100"); do

	$EXEC $p "$MESH"$N".msh" "$OUTBASE"$N
	
    done
    
done
