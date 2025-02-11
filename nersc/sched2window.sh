#!/bin/sh

ROOT=/home/herault/Documents/Recherche/stochastic-io/experiments
APEX=$ROOT/apex-new.py
SIMU=$ROOT/simu
PLAT=nersc
OUTPUT=$ROOT/apex/${PLAT}/all-$$.csv
MYSELF=$ROOT/apex/${PLAT}/sched2window.sh

onerun() {
    baserun=$1
    if [ -f done ]; then
        return
    fi
    for s in *.sched; do
        obase=$(basename $s .sched)
        if [ ! -f ${obase}-020.run ]; then
          python3 $APEX -m window -i $s -o "${obase}.txt"
        fi
        for r in ${obase}*.run; do
            $SIMU $r
        done
    done
    touch done
}

allrun() {
    for d in ${PLAT}-*; do
        baserun=$(echo $d | awk -v FS='[/-]' '{print $2}')
        cd $d
        flock -w 0 . sh $MYSELF onerun $baserun
        cd ../
    done
}

case "x$1" in
  "x" )
    exec nohup sh $MYSELF allrun > $ROOT/apex/${PLAT}/runner.$$ 2>&1 &
  ;;
  "xallrun" )
    allrun 
  ;;
  "xonerun" )
    onerun $2
  ;;
esac
