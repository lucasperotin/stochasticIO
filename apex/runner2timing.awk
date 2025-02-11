#!/usr/bin/awk

BEGIN {
    OFS=","
    print("plat,machine,method,version,iopressure,time")
}

$1=="TIMING:" { 
    file=$2; 
    plat=file; sub(/-.*$/, "", plat); 
    machine=file; 
    sub(/^[^-]*-/, "", machine); 
    sub(/-.*$/, "", machine); 
    method=$3; 
    version=$4; 
    time=$5 
} 

$1=="Communication" {
    if( plat != "" )
        print(plat, machine, method, version, $6, time);
    plat=""
}
