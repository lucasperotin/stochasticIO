#!/bin/bash

SHELL=$(ps -p $$ | awk '$4~/bash/ {ok=1} $4~/run-apex.sh/ {ok=1} {last=$4} END {if(ok) print "OK"; else print last}')
if [ "x$SHELL" != "xOK" ]; then
    echo "script won't work with $SHELL shell -- please run with bash" >> /dev/stderr
    exit 1
fi

PYTHON=$(which python3)
if [ -z "$PYTHON" ]; then
    PYTHON=$(which python)
    if [ -z "$PYTHON" ]; then
        echo "Couldn't find a python3 interpreter" >> /dev/stderr
        exit 1
    fi
fi

if [ ! -f ./apex.py ]; then
    echo "Couldn't find ./apex.py" >> /dev/stderr
    exit 1
fi

for PLATFORM in nersc trilab; do

    [ -d ./${PLATFORM} ] || mkdir -p ./${PLATFORM}

    BASE="./${PLATFORM}/${PLATFORM}"
    SCHED="${BASE}.sched"
    GANTT="${BASE}.sched.png"
    WINDOW="${BASE}.txt"
    RUN="./${PLATFORM}.run"

    PROCEED="y"
    if [ -f "${SCHED}" ]; then
        read -iN -r -n1 -p "${SCHED} exists do you want to create a new schedule for ${PLATFORM} (y/*)? " PROCEED
        echo
    fi
    if [ "x${PROCEED}" == "xy" ]; then
        echo "Generating Schedule for ${PLATFORM}"
        echo " -- Running $PYTHON ./apex.py -m schedule -p ${PLATFORM} -o ${SCHED}"
        $PYTHON ./apex.py -m schedule -p ${PLATFORM} -o ${SCHED}
        echo "Generating Gantt chart of the schedule"
        echo " -- Running $PYTHON ./apex.py -m draw -i ${SCHED} -o ${GANTT}"
        $PYTHON ./apex.py -m draw -i ${SCHED} -o ${GANTT}
        echo "Removing old windows for ${PLATFORM}"
        rm -f ${BASE}-*.txt
    fi

    PROCEED="y"
    NWINDOWS=$(ls ${PLATFORM}/${PLATFORM}-*.txt | wc -l)
    if [ ${NWINDOWS} -gt 0 ]; then
        read -p "There are ${NWINDOWS} window files, do you want to re-create a new set of windows for ${PLATFORM}? (y/*)" -a PROCEED -i "N" -n 1 -r
        echo
    fi
    if [ "x${PROCEED}" == "xy" ]; then
        echo "Removing previous windows"
        rm -f ${PLATFORM}/${PLATFORM}-*.txt
        echo "Generating windows:"
        echo " -- Running $PYTHON ./apex.py -m window --window 0 -i ${SCHED} -o ${WINDOW}"
        WINSIZE=$($PYTHON ./apex.py -m window --window 0 -i ${SCHED} -o ${WINDOW})
        echo " -- Running $PYTHON ./apex.py -m window --window ${WINSIZE} -i ${SCHED} -o ${WINDOW}"
        $PYTHON ./apex.py -m window --window ${WINSIZE} -i ${SCHED} -o ${WINDOW}
    fi

    echo "Generating ${RUN}. Confirm file suppression if needed."
    [ -f ${RUN} ] && rm -i ${RUN}

    WINDOWS=$(ls ${PLATFORM}/${PLATFORM}-*.txt)

    NSAMPLES=0
    for m in ${WINDOWS}; do
        NSAMPLES=$((NSAMPLES+1))
    	for version in v1 v2; do
            for method in fairShare FCFS greedyYield greedyCom Set10Learn fixedWindow lookAheadGreedyYield nextEvLexMin; do
		for w in 0.2 0.5 0.8 0.9 1.0 1.1; do
	            echo "$m 214748364800 1 $method $version results/results${PLATFORM}/outw$w-$version.txt" >> ${RUN}
		done
            done
	done
    done
    echo "APEX: ${PLATFORM}/${NSAMPLES} traces generated"
    echo "Running ${RUN} with ./simu to generate output"

    if [ -x ./simu ]; then
        ./simu ${RUN}
    else
        echo "Couldn't find ./simu -- did you compile it?" >> /dev/stderr
    fi
done
