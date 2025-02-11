#!/bin/sh

echo "plat,nodef,coref,memf,nodebwf,sysbwf,run,Number of Applications,Window Duration (s)"
awk -v OFS=',' -v FS='[ -]+' '$1=="INFO:root:Window" && $3=="duration" {wd=$4} $1=="INFO:root:Window" && $3=="contains" {run=$13; sub(/.txt$/, "", run); print($7,$8,$9,$10,$11,$12,run,$4, wd)}' $*
