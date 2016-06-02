#!/usr/bin/python

import sys, getopt, os, numpy, glob, shutil
import matplotlib.pyplot as plt

# parameters for plotting
params1 = {	'legend.fontsize' : 18,
		'axes.titlesize' : 20, 	#Title 
		'axes.labelsize' : 18	#Axenbeschriftung
	  }
plt.rcParams.update(params1)

# checks if the results directory already exists
def checkdir():
    print "checking directory"

# creates a new directory for the reynolds number results
def mkdir(Re):
    print "creating directory"
# print the combined results into one file
# REFACTOR
def printRes(alpha, Cd, Cl, Cd_std, Cl_std, Re):		 #writes output data files
    print "printing Results"
    

# creates matplotlib plots with the combined data
def plot_fig (alpha, Cd, Cl, Cd_std, Cl_std, Re):
    print "plotting figures"
    
### reference function to convert mV/V to N
# REFACTOR to generic regression analysis function
# see commentary above
def mVtoC(param, Cl, Ac, Uc):
    print "mv to Cx"
    
# REFACTOR with get reference lift
def getCalibration(files):
    print "getting calibration references"

def parseLine (line):
    print "parsing lines"

# refactor the parameters -> instead of allfiles the correct file path
def average(allfiles, key, Re):                  #param = Dictionary and Key
    print "averaging"

def getKey(item):
    print "getting key"
    
def merge(files, Re):
    print "merging files"
    
def getData(filepaths):
    print "collecting all input data"

def getCalibration():

    

    print "calculating calibration"
    
def getMount():
    print "getting mount offset"

def getZero():
    print "getting zero offset"

def getDisk():
    print "getting disk offset"

#### START OF THE MAIN LOOP
def main (argv):

    cal_files = True;
    mount_files = False;
    zero_files = False;
    disk_files = False;

    print "getting user input"

    getData("/tmp/test")

    if (cal_files):
        getCalibration()

    if (mount_files):
        getMount()
        
    if (zero_files):
        getZeor()
    
    if (disk_files):
        getDisk()
    
    
    
    

if __name__ == "__main__":
    main(sys.argv[1:])
