#!/bin/bash

# todo: write a script taht sets up everything automatically and detects runnging OpenFOAM versions.

if [ -d $HOME/bin ]
then
	echo "symlinkinking scripts to $HOME/bin"
	ln -s $(pwd)/src/mesh/meshgen.m $HOME/bin/meshgen
	ln -s $(pwd)/src/postprocessing/forceCoeffsAnalysis.py $HOME/bin/forceCoeffsAnalysis
	ln -s $(pwd)/src/postprocessing/forceAnalysis.py $HOME/bin/forceAnalysis
else 
	echo "$HOME/bin does not exists, dying"
fi

source $HOME/.bashrc
