#!/bin/bash
#Script: testMcCode.sh
#Purpose: sets up the required directory structure on scratch for a test run of the
#         monte carlo code. Used for testing and debugging of new features.
#Author: Reese Horton

if [ ! -d "$MYSCRATCH"/TEST ]; then
   mkdir "$MYSCRATCH"/TEST
fi

if [ ! -d "$MYSCRATCH"/data-1 ]; then
   echo "Copying input data from Acacia"
   mc cp --recursive rhorton/data-1 "$MYSCRATCH"
fi

cd "$MYSOFTWARE"/MCSimulation
cp ./main "$MYSCRATCH"/TEST
cp ./data.in "$MYSCRATCH"/TEST
