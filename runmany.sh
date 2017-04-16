#!/bin/bash

function myWait {
	local PIDS=("$@")
	local SLEEP=1
	while [ "x$SLEEP" == "x1" ]
	do
		SLEEP=0
		for p in "${PIDS[@]}"
		do
			local run=$(ps -ef | grep "$p" | wc -l)
			if [ "$run" -ge "1" ]
			then
				SLEEP=1
			fi
		done
		
		if [ "x$SLEEP" == "x1" ]
		then
			sleep 1
		fi
	done
}

function testAvgN {
	local N=$1
	local E=$2
	local A=("$@")
	local TEST_NAME="avgNTest"
	local PID=""
	local FILENAMES=""
	for i in "${A[@]:2}"
	do
		local FN="${TEST_NAME}_${N}_${i}.csv"
		while [ -e "./$FN" ]
		do
			FN="_$FN"
		done
		FILENAMES="$FILENAMES $FN"
		local K=$[i+10]
		start bash -c "echo \$\$>tt;python ./simulation.py -n$N -e$E -a$i -k$K -o$FN"
		sleep 0.05
		local P=`cat tt`
		PID="$PID $P"
		rm tt
	done	
	myWait $PID
	echo "Done $TEST_NAME $@"
	python plot.py $FILENAMES &
}

function testN {
	local A=$1
	local K=$[A+10]
	local E=$2
	local N=("$@")
	local TEST_NAME="NTest"	
	local PID=""
	local FILENAMES=""
	for i in "${N[@]:2}"
	do		
		local FN="${TEST_NAME}_${i}_${A}.csv"
		while [ -e "./$FN" ]
		do
			FN="_$FN"
		done
		FILENAMES="$FILENAMES $FN"
		start bash -c "echo \$\$>tt;python ./simulation.py -n$i -e$E -a$i -k$K -o$FN"
		sleep 0.05
		local P=`cat tt`
		PID="$PID $P"
		rm tt
	done
	myWait $PID
	echo "Done $TEST_NAME $@"
	python plot.py $FILENAMES &
}

# configure simulations here

testN 10 200 500 1000 2000 5000 10000
# simultaneously runs 5, 200 frame simulations with avgN=10 and n=500,1000,2000,5000,10000


# it will wait until those 5 simulations have all finished, then it will run
testAvgN 2000 200 10 15 20 25 30 
# this simultaneously runs 5, 200 frame simulations with n=2000 and avgN=10:5:30

wait
echo "Done All"