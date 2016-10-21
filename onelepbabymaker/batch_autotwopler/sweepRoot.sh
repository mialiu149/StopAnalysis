#!/bin/bash

FILENAME=$1

echo "    sweepRooting: \"$FILENAME\""
RESULT=$(./sweepRoot -b -o "t" $FILENAME | grep "BAD FILES")

if [ [$RESULT == *"BAD FILES"*] ]; then 
    echo $RESULT   
    exit 1
else
    echo 'good file'   
    exit 0
fi
