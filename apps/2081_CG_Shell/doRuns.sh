#!/bin/bash

EXEC="$(pwd)/../../cmake-build-release/apps/shellProject"
MESH="thinPlate"
OUTBASE="output"

p=1
$EXEC $p "$MESH""_4.msh" "$OUTBASE"_4
$EXEC $p "$MESH""_10.msh" "$OUTBASE"_10
$EXEC $p "$MESH""_20.msh" "$OUTBASE"_20
$EXEC $p "$MESH""_100.msh" "$OUTBASE"_100

p=2
$EXEC $p "$MESH""_4.msh" "$OUTBASE"_4
$EXEC $p "$MESH""_10.msh" "$OUTBASE"_10
$EXEC $p "$MESH""_20.msh" "$OUTBASE"_20

p=3
$EXEC $p "$MESH""_4.msh" "$OUTBASE"_4
$EXEC $p "$MESH""_10.msh" "$OUTBASE"_10
$EXEC $p "$MESH""_20.msh" "$OUTBASE"_20

p=4
$EXEC $p "$MESH""_4.msh" "$OUTBASE"_4
$EXEC $p "$MESH""_10.msh" "$OUTBASE"_10
$EXEC $p "$MESH""_20.msh" "$OUTBASE"_20

