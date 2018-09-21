#!/bin/bash

if [ $# -lt 5 ]
then
    echo "USAGE: thread_switch.sh <seq> <start_res> <end_res> <input_structure> <username>"
    exit 1
fi

MAXJOBS=20  #number of jobs for your server
seq=$1
start=$2
end=$3
structure=$4
username=$5

for j in $(seq $start $end)
do
    activeJobs=$(ps -U $username | grep -c "rosetta_scripts")
    while [ $activeJobs -ge $MAXJOBS ]
    do
        # wait for 20 seconds
        sleep 20s
        activeJobs=$(ps -U $username | grep -c "rosetta_scripts")
        echo $activeJobs
    done #end_of_while_loop

    echo "Submitting $j"
    rosetta_scripts.linuxgccrelease -s $structure -parser:protocol ~/scripts/rosetta_scripts/thread_relax.xml -parser:script_vars start=$j seq=$seq -out:prefix "$j"_threaded_ &> $j"_log" &
done
