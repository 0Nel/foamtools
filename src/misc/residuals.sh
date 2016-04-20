#!/bin/bash

working_dir="/tmp/residuals"
gnuplot_format_file="$HOME/lib/gnuplot/Residuals"

function usage {
	cat << EOF
$0 caseName casePath

$0 is a tool to monitor residuals on a remote
host by using gnuplot and scp
EOF
	exit
}

if [ $# -ne 1 ]
then
	echo "please provide a case name and path"
	usage
fi

case_name=$1
case_path="/work/asander/"

working_dir="$working_dir/$case_name"

### check for mandatory gnuplot format file
if [ -e $gnuplot_format_file ]
then
	mkdir -p $working_dir;
	cp $gnuplot_format_file $working_dir
	cd $working_dir
	pwd
else 
	echo "could not find plotting template"
	exit
fi

scp -2 $USERNAME\@$HOSTNAME:$case_path/$case_name/$logfile .

gnuplot Residuals

while :
do
	scp -2 $USERNAME\@$HOSTNAME:$case_path/$case_name/$logfile .
	sleep 30s
done
