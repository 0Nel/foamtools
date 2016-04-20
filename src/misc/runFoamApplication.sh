#!/bin/bash

if [ $1 ]
then
	APPLICATION=$1
else
	echo "please provide the application you would like to run as a command line parameter"
	exit
fi

N_RESTART=100
COUNTER=0

echo "number of restarts is set to $N_RESTARTS"

echo "running $APPLICATION"


while [ $COUNTER -lt $N_RESTART  ]
do
	$APPLICATION | tee log.$APPLICATION
	echo "restarting $APPLICATION for the $N_RESTART time"
done

echo "$APPLICATION was restartet $COUNTER times. pleae check the case!"
