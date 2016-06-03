#!/usr/bin/python

import sys, getopt, os, numpy, glob, shutil
import matplotlib.pyplot as plt

# parameters for plotting
plotParameters = { 'legend.fontsize' : 18,
					'axes.titlesize' : 20, 	#Title
					'axes.labelsize' : 18	#Axenbeschriftung
	  			 }
plt.rcParams.update(plotParameters)

fluidParameters = 	{
						'Ac' : 1,
						'Uc' : 1,
						'Rho' : 1
					}

def parseLineToFloats(line):
	try:
		tmp = [float(string.replace(',', '.')) for string in line.split()]
		if tmp:
			return tmp
	except:
		return None

def parseLineToDicts(line):
	print "parsing line to dicts"
	# try:


# parses a file and return a 2D array
# if the file is empty or no floats are found the function
# will return None. The same happens if the file can't be opened
def parseFile (filepath):
	# try opening the input file
	try:
		f_in = open(filepath, 'r')
	except(IOError):
		print "could not open file:", filepath
		return None

	tmp = []		# temporary array for parsing lines
	tmp2D = []		# temporary array for building a 2D array
	for line in f_in:
		tmp = parseLineToFloats(line)
		if tmp != None:
			tmp2D.append(tmp)
			tmp = []

	if tmp2D:
		return numpy.array(tmp2D)
	else:
		return None

def parseConfigFile(filepath = '.windkanal.config'):
	try:
		f_in = open(filepath, 'r')
	except(IOError):
		print "could not open file:", filepath
		return None

	#for line in f_in:


#### START OF THE MAIN LOOP
def main (argv):
	file = raw_input("please provide a file: ")
	print parseFile(file)


if __name__ == "__main__":
    main(sys.argv[1:])
