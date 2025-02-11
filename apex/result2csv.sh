#!/bin/sh

echo "plat,nodef,coref,memf,nodebwf,sysbwf,run,method,version,iopressure,minyield,efficiency,utilization"

rm -f errors.log

for r in $(find nersc/ trilab/ -name *.res); do
    BASERUN=$(echo $r | awk -v FS='[-/]' '{n=$10; sub(/.res/,"",n); printf("%s%s", $3, n)}')
    NODEF=$(echo $r | awk -v FS='[-/]' '{print $5}')
    COREF=$(echo $r | awk -v FS='[-/]' '{print $6}')
    RAMF=$(echo $r | awk -v FS='[-/]' '{print $7}')
    NODEBWF=$(echo $r | awk -v FS='[-/]' '{print $8}')
    SYSBWF=$(echo $r | awk -v FS='[-/]' '{print $9}')
    PLAT=$(echo $r | awk -v FS='[-/]' '{print $1}')
    awk -v OFS=',' '{if($4>0) print "'$PLAT'", "'$NODEF'", "'$COREF'", "'$RAMF'", "'$NODEBWF'", "'$SYSBWF'", "'$BASERUN'", $1, $2, $3, $4, $5, $6; else print("'$r'") >> "errors.log" }' $r
done
