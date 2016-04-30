#!/usr/bin/python

import sys, getopt, os, numpy, math, glob
# import re
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})

def nextpow2 (i):
    n = 1
    while n < i: n *= 2
    return n

def plot_fig (alpha, Cl, Cl_std, Cd, Cd_std):

    plt.figure(3)
    plt.title('Auftriebsbeiwert')
    plt.plot(alpha, Cl)
    plt.grid()
    plt.xlabel("alpha")
    plt.ylabel("Cl")
    plt.savefig("Cl.png", dpi=300)

    plt.figure(4)
    plt.title('Widerstandsbeiwert')
    plt.plot(alpha, Cd)
    plt.grid()
    plt.xlabel("alpha")
    plt.ylabel("Cd")
    plt.savefig("Cd.png", dpi=300)

    plt.figure(5)
    plt.title('Lilienthalpolare')
    plt.plot(Cd, Cl)
    plt.grid()
    plt.xlabel("Cd")
    plt.ylabel("Cl")
    plt.savefig("lilienthal.png", dpi=300)

### reference function to convert mV/V to N
def mVtoDrag(param, Cd):
    tmp = param[0] * Cd + param[1]
    
def mVtoLift(param, Cl):
    return param[0] * Cl + param[1]

def getReferenceDrag(filename):

    weight = []
    force = []
    force_std = []

    f_in = open(filename, 'r')
    for line in f_in:
        tmp = line.split()
        if not tmp[0] == '#':
            weight.append(float(tmp[0]))
            force.append(float(tmp[1]))
            force_std.append(float(tmp[2]))
    
    weight2D = numpy.array([weight, numpy.ones(len(weight))])
    
    param = numpy.linalg.lstsq(weight2D.T, force)[0]
   
    ## example code for plotting
    plt.figure(1)
    plt.title('linear regression drag')
    plt.plot(weight, force)
    plt.grid()
    plt.xlabel("weight")
    plt.ylabel("force")
    plt.savefig("drag_reference.png", dpi=300)
    
    return param

def getReferenceLift(filename):

    weight = []
    force = []
    force_std = []

    f_in = open(filename, 'r')
    for line in f_in:
        tmp = line.split()
        if not tmp[0] == '#':
            weight.append(float(tmp[0]))
            force.append((tmp[1]))
            force_std.append(tmp[2])
    
    weight2D = numpy.array([weight, numpy.ones(len(weight))])
    
    param = numpy.linalg.lstsq(weight2D.T, force)[0]
    
    ## example code for plotting
    plt.figure(2)
    plt.title('linear regression lift')
    plt.plot(weight, force)
    plt.grid()
    plt.xlabel("weight")
    plt.ylabel("force")
    plt.savefig("lift_reference.png", dpi=300)

    return param

#### START OF THE MAIN LOOP
def main (argv):

## positions of the data in the files
    ALPHA = 0

    CD = 1
    CD_STDV = 2

    CL = 3
    CL_STDV = 4

    inputfiles = []         # array to save the file path
    outputfile = ''         # output file path
    parameterfile = 'windtunnel.inputrc'
    drag_reference = ''
    lift_reference = ''

    coeffs = False          # calculate coeffs?

    cd_axis = 'x'
    cl_axis = 'y'
    cm_axis = 'z'

    Ac = 1.0                # Reference area, default value
    Uc = 1.0                # Reference speed
    Lc = 1.0                # Reference length

    conversion_drag = [1.0, 1.0]
    conversion_lift = [1.0, 1.0]

    # 2 * Fw / rho * u_c^2  * A

    ### now get the command line parameters provided by the user

    try:
        opts, args = getopt.getopt(argv, "i:o:v:a:d:l:", ["input-files=", "output-file=", "inlet-velocity=", "reference-area", "coeffs", "drag-reference", "lift-reference"])
    except getopt.GetoptError:
        print 'windtunnel.py -i <inputfiles> [-o <outputfile> -v <velocity> -r <reference area> -c -d <drag reference file> -l <lift reference file>]'
        sys.exit(-1)

    for opt, arg in opts:
        if opt in ("-i", "--input-file"):
            if (os.path.isfile(arg)):
                print "using ", arg, " as input"
                inputfiles.append(arg)
            elif (len(arg) > 0):
                inputfiles = glob.glob(arg)
                for file in inputfiles:
                    if (not os.path.isfile(file)):
                        print "not a valid file: ", file
                        sys.exit(-1)
            else:
                print "no valid input.. dying."
                sys.exit(-1)
        elif opt in ("-o", "--output-file"):
            outputfile = arg
        elif opt in ("-v", "--inlet-velocity"):
            print "reference velocity set to: ", arg
            Uc = arg
        elif opt in ("-a", "--reference-area"):
            print "reference area set to ", arg
            Ac = arg
        elif opt in ("-d", "--drag-reference"):
            if (os.path.isfile(arg)):
                print "using ", arg, " as drag reference file"
                conversion_drag = getReferenceDrag(arg)
                print "conversion drag: ", conversion_drag
            else:
                print "not a valid file: ", arg, "exiting"
                sys.exit(-1)
        elif opt in ("-l", "--lift-reference"):
            if (os.path.isfile(arg)):
                print "using ", arg, " as lift reference file"
                conversion_lift = getReferenceLift(arg)
                print "conversion lift: ", conversion_lift
            else:
                print "not a valid file: ", arg, "exiting"
                sys.exit(-1)

    raw2D = []    # leeres array initialisieren
    raw3D = []
    file_counter = 0

    for file in inputfiles:
        f_in = open(file, 'r')
        for line in f_in:
            # tmp = [x.strip('(').strip(')') for x in line.split()]
            tmp = [x.strip('(').strip(')') for x in line.split()]
            if tmp[0] == '#':           # skip comments
                # print "skipping line: ", line
                next
            else:
                raw2D.append(map(float, tmp))
        raw3D.append(raw2D)
        raw2D = []

    # now that all the data is imported from the inutfiles, convert the array to a numpy array
    raw3D = numpy.array(raw3D)

    # get the references for lift and drag and calculate the linear regression
    

    '''
    the data now needs to be cooked with some voodo:
    at the moment the 3rd dimension is depended on the different input filesi
    the array needs to be transposed in order to easily iterate over the same alphas for
    different files
    '''

    raw3DT = numpy.transpose(raw3D)     # its not that much voodo. actually its no vodoo at all,..

    averaged = []
    averaged2D = []
    for file in raw3DT:
        for alpha in file:
            averaged.append(numpy.average(alpha))
        averaged2D.append(averaged)
        averaged = []
    averaged2D = numpy.array(averaged2D)

    '''
    now we need to transpose it back in order to print it out
    '''

    averaged2DT = numpy.transpose(averaged2D)

    f_out = open('averaged_data_raw.dat', 'w')

    for alpha in averaged2DT:
        for value in alpha:
            f_out.write(str(value))
            f_out.write(" ")
        f_out.write ("\n")
    f_out.close()

    alpha = averaged2DT[:,0]

    Cd = mVtoDrag(conversion_drag, averaged2DT[:,1])
    Cd_std = mVtoDrag(conversion_drag, averaged2DT[:,2])
    
    Cl = mVtoLift(conversion_drag, averaged2DT[:,3])
    Cl_std = mVtoLift(conversion_drag, averaged2DT[:,4])

    plot_fig(alpha, Cd, Cd_std, Cl, Cl_std)

if __name__ == "__main__":
    main(sys.argv[1:])
