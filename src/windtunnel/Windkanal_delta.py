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

def parseStringToFloat(string):
	try:
		return float(string.replace(',', '.'))
	except:
		return string

def parseLineToFloats(line):
	try:
		tmp = [parseStringToFloat(string) for string in line.split()]
		if tmp:
			return tmp
	except:
		return None

def parseLineToDicts(line):

	if len(line) == 0:
		return

	if line.find('#') == 0:
		return

	tmp = line.split()
	if len(tmp) > 0:
		# the general structure is supposed to be like this:
		# NAME VALUE
		# by name beeing a normal string (which we also don't touch)
		# and value either being a float or a filepath.
		if (os.path.isfile(tmp[1])):
			return {tmp[0]:tmp[1]}
		else:
			return {tmp[0]:parseStringToFloat(tmp[1])}
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

def parseConfigFile(fluidParameters, filepath = '.windkanal.config'):
	try:
		f_in = open(filepath, 'r')
	except(IOError):
		print "could not open file:", filepath
		return

	#for line in f_in:
	for line in f_in:
		try:
			fluidParameters.update(parseLineToDicts(line))
		except:
			return

#### START OF THE MAIN LOOP
def main (argv):
	file = raw_input("please provide a file: ")
	print parseFile(file)
	print fluidParameters
	print parseConfigFile(fluidParameters, '.windkanal.config')
	print fluidParameters


if __name__ == "__main__":
    main(sys.argv[1:])
