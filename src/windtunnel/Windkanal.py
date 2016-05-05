#!/usr/bin/python

import sys, getopt, os, numpy, glob
# import re
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})

def nextpow2 (i):
    n = 1
    while n < i: n *= 2
    return n
 
def plot_fig (alpha, Cd, Cd_std, Cl, Cl_std):

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

    plt.figure(4)
    plt.title('Lilienthalpolare')
    plt.plot(Cd, Cl)
    plt.grid()
    plt.xlabel("Cd")
    plt.ylabel("Cl")
    plt.savefig("lilienthal.png", dpi=300)

### reference function to convert mV/V to N
def mVtoDrag(param, Cd, Ac, Uc):
    '''TODO: Catch Cd vs. Cd_std'''
    try:
        Uc = float(Uc)
        Ac = float(Ac)
        return (((param[0] * Cd + param[1])*2*0.009810)/(Ac*Uc**2))
    except(TypeError):
        print('Invalid reference Area or Velocity!')
        sys.exit(-1)

def mVtoLift(param, Cl, Ac, Uc):
    try:
        Uc = float(Uc)
        Ac = float(Ac)
        return (((param[0] * Cl + param[1])*2*0.009810)/(Ac*Uc**2))
    except(TypeError):
        print('Invalid reference Area or Velocity!')
        sys.exit(-1)

def getReferenceDrag(filename):

    weight = []
    force = []
    force_std = []

    f_in = open(filename, 'r')
    for line in f_in:
        tmp = line.split()
        if len(tmp) > 0 and not tmp[0] == '#':
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

def checkLine(line):
    if 'Stepper' in line:                   #besser bei allen moeglichen buchstaben?
        return True
    else:
        return False

def average(allfiles, key):                  #param = Dictionary and Key
    try:
        files = allfiles[key]    
        raw2D = []    # leeres array initialisieren
        raw3D = []
        for file in files:
            f_in = open(file, 'r')            
            for line in f_in:
                skip = checkLine(line)
                if skip:
                    next
                else:
                    line = line.replace(',','.')
                    tmp = line.split()
                    raw2D.append(map(float, tmp))
            raw3D.append(raw2D)
            raw2D = []

        # now that all the data is imported from the inutfiles, convert the array to a numpy array
        raw3D = numpy.array(raw3D)
        raw3DT = numpy.transpose(raw3D)     # its not that much voodo. actually its no vodoo at all,..
    
        averaged = []
        averaged2D = []
        for file in raw3DT:
            for alpha in file:
                averaged.append(numpy.average(alpha))
            averaged2D.append(averaged)
            averaged = []
        averaged2D = numpy.array(averaged2D)
        
        averaged2DT = numpy.transpose(averaged2D)
    
        f_out = open(key + '_averaged_data_raw.dat', 'w')
    
        for alpha in averaged2DT:
            for value in alpha:
                f_out.write(str(value))
                f_out.write(" ")
            f_out.write ("\n")
        f_out.close()
        
        return averaged2DT
        
    except KeyError :
        print key + ' is not a valid key! Please provide ' + key + ' offset file!'    
        sys.exit(-1)

#### START OF THE MAIN LOOP
def main (argv):

## positions of the data in the files
    ALPHA = 0

    CD = 3
    CD_STDV = 4

    CL = 1
    CL_STDV = 2

    inputfiles = []         # array to save the file path
    zerofiles = []         # array to save the file path
    mountfiles = []         # array to save the file path
    allfiles = {}         # dictionary to save the all file arrays and set a keyword
    
    outputfile = ''         # output file path

    Ac = 1.0                # Reference area, default value
    Uc = 1.0                # Reference speed
    Lc = 1.0                # Reference length

    conversion_drag = [1.0, 1.0]
    conversion_lift = [1.0, 1.0]

    # 2 * Fw / rho * u_c^2  * A

    ### now get the command line parameters provided by the user

    try:
        opts, args = getopt.getopt(argv, "i:o:v:a:d:l:z:m:", ["input-files=", "output-file=", "inlet-velocity=", "reference-area", "coeffs", "drag-reference", "lift-reference", "zero-offset", "mount-offset"])
    except getopt.GetoptError:
        print 'windtunnel.py -i <inputfiles> [-o <outputfile> -v <velocity> -a <reference area> -c -d <drag reference file> -l <lift reference file> -z <zero offset> -m <mount offset>]'
        sys.exit(-1)

    for opt, arg in opts:
        if opt in ("-i", "--input-file"):
            if (os.path.isfile(arg)):                       #checks existence of inputfile
                print "using ", arg, " as input file"
                inputfiles.append(arg)
            elif (len(arg) > 0):                            #catches and appends all inputfiles
                inputfiles = glob.glob(arg)
                if (len(inputfiles) == 0):                  #catch Schreibfehler im Pfad
                    print "not a valid file: ", arg
                    sys.exit(-1)                
                for i in range(0, len(inputfiles)):
                    print 'using', inputfiles[i], " as input file"
                for file in inputfiles:
                    if (not os.path.isfile(file)):
                        print "not a valid file: ", file
                        sys.exit(-1)
            else:
                print "no valid input.. dying."
                sys.exit(-1)
            allfiles['input'] = inputfiles
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
        elif opt in ("-z", "--zero-offset"):
            if (os.path.isfile(arg)):                       #checks existence of inputfile
                print "using ", arg, " as zero offset file"
                zerofiles.append(arg)
            elif (len(arg) > 0):                            #catches and appends all inputfiles
                zerofiles = glob.glob(arg)
                if (len(zerofiles) == 0):
                    print "not a valid file: ", arg
                    sys.exit(-1)
                for file in zerofiles:
                    if (not os.path.isfile(file)):
                        print "not a valid file: ", arg
                        sys.exit(-1)
                for i in range(0, len(zerofiles)):
                    print 'using', zerofiles[i], " as zero offset file"
            else:
                print "not a valid file: ", arg, "exiting"
                sys.exit(-1)
            allfiles['zero'] = zerofiles
        elif opt in ("-m", "--mount-offset"):
            if (os.path.isfile(arg)):                       #checks existence of inputfile
                print "using ", arg, " as mount offset file"
                mountfiles.append(arg)
            elif (len(arg) > 0):                            #catches and appends all inputfiles
                mountfiles = glob.glob(arg)
                if (len(mountfiles) == 0):
                    print "not a valid file: ", arg
                    sys.exit(-1)                
                for i in range(0, len(mountfiles)):
                    print 'using', mountfiles[i], " as mount offset file"
                for file in mountfiles:
                    if (not os.path.isfile(file)):
                        print "not a valid file: ", file
                        sys.exit(-1)
            else:
                print "not a valid file: ", arg, "exiting"
                sys.exit(-1)
            allfiles['mount'] = mountfiles
            
    data = average(allfiles, 'input')
    zero = average(allfiles, 'zero')
    mount = average(allfiles, 'mount')

#### Mehrere Kalibrierdateien einlesen, sortieren und dann erst verarbeiten

    if numpy.average(zero.T[0]) == 0:
        zeroLift = numpy.average(zero.T[1])
        zeroDrag = numpy.average(zero.T[3])
        data.T[1] = numpy.subtract(data.T[1], zeroLift)
        data.T[3] = numpy.subtract(data.T[3], zeroDrag)
    else:
        print 'Mount Offset has an angle of attack greater than zero.'

    if numpy.average(mount.T[0]) == 0:
        mountLift = numpy.average(mount.T[1])
        mountDrag = numpy.average(mount.T[3])
        data.T[1] = numpy.subtract(data.T[1], mountLift)
        data.T[3] = numpy.subtract(data.T[3], mountDrag) 
    else:
        print 'Mount Offset has an angle of attack greater than zero.'

    alpha = data[:,ALPHA]

    Cd = mVtoDrag(conversion_drag, data[:,CD],Ac,Uc)
    Cd_std = mVtoDrag(conversion_drag, data[:,CD_STDV],Ac,Uc)
    Cl = mVtoLift(conversion_drag, data[:,CL], Ac, Uc)
    Cl_std = mVtoLift(conversion_drag, data[:,CL_STDV], Ac, Uc)

    plot_fig(alpha, Cd, Cd_std, Cl, Cl_std)

if __name__ == "__main__":
    main(sys.argv[1:])
