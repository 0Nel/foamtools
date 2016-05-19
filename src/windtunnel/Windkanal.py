#!/usr/bin/python

import sys, getopt, os, numpy, glob
# import re
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})

def nextpow2 (i):
    n = 1
    while n < i: n *= 2
    return n
 
def plot_fig (alpha, Cd, Cl):	# Cd_std, Cl, Cl_std):

    plt.figure(3)
    plt.title('Auftriebsbeiwert')
    plt.plot(alpha, Cl, color = 'red')
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
def mVtoDrag(param, Cd, Ac, Uc):
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

    f_in = open(filename, 'r')
    for line in f_in:
        tmp = line.split()
        if len(tmp) > 0 and not tmp[0] == '#':
            weight.append(float(tmp[0]))
            force.append(float(tmp[1]))

    force2D = numpy.array([force, numpy.ones(len(force))])
    
    param = numpy.linalg.lstsq(force2D.T, weight)[0]

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

    f_in = open(filename, 'r')
    for line in f_in:
        tmp = line.split()
        if not tmp[0] == '#':
            weight.append(float(tmp[0]))
            force.append((tmp[1]))

    force2D = numpy.array([force, numpy.ones(len(force))])

    param = numpy.linalg.lstsq(force2D.T, weight)[0]

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
    line = line.replace(',','.')
    tmp = line.split()
    try:
        float(tmp[0])
        return False
    except (ValueError):
        return True

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

def getKey(item):
	return item[0]
    
def merge(files):
    if len(files) > 2:
        print "Too many reference files given: ", files
    elif len(files[0]) == len(files[1]):
        raw0 = []
        raw1 = []
        count = False
        for file in files:
            f_in = open(file, 'r')            
            if 'LIFT' in f_in.readline():
                name = 'LIFT'
            else:
                name = 'DRAG'
            for line in f_in:
                skip = checkLine(line)
                if skip:
                    next
                else:
                    line = line.replace(',','.')
                    tmp = line.split()
                    if (count == False):
                        raw0.append(map(float, tmp))
                    else:
                        raw1.append(map(float, tmp))
            count = True       

        raw0 = numpy.array(sorted(raw0, key = getKey))
        raw1 = numpy.array(sorted(raw1, key = getKey))
       
        a = raw0.T[1]
        b = raw1.T[1]  

        t = raw0.T[0]
        f_out = open(name.lower() + '_calibration_averaged_data_raw.dat', 'w')
        for i in range(0, len(raw0)):
            f_out.write(str(t[i]) + ' ' + str(((a[i]+b[i])/2)) +'\n')
        f_out.close()
        
        return (name.lower() + '_calibration_averaged_data_raw.dat')
    else:
          print "Reference files have inconsistent lengths"
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
    dragfiles = []      #array to save drag calibration measurements
    liftfiles = []      #array to save drag calibration measurements
    allfiles = {}         # dictionary to save the all file arrays and set a keyword

    Ad = 1.0                # Reference area drag, default value
    Al = 1.0                # Reference area lift, default value
    Uc = 1.0                # Reference speed
#    Lc = 1.0                # Reference length

    conversion_drag = [1.0, 1.0]
    conversion_lift = [1.0, 1.0]

    # 2 * Fw / rho * u_c^2  * A

    ### now get the command line parameters provided by the user

    try:
        opts, args = getopt.getopt(argv, "i:v:s:a:d:l:z:m:", ["input-files=", "inlet-velocity=", "reference-area-drag", "reference-area-lift", "drag-reference", "lift-reference", "zero-offset", "mount-offset"])
    except getopt.GetoptError:
        print 'windtunnel.py -i <inputfiles> [-v <velocity> -s <reference area drag> -a <reference area lift> -c -d <drag reference file> -l <lift reference file> -z <zero offset> -m <mount offset>]'
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
        elif opt in ("-d", "--drag-reference"):
            if (os.path.isfile(arg)):                       #checks existence of inputfile
                print "using ", arg, " as drag reference file"
                dragfiles.append(arg)
            elif (len(arg) > 0):                            #catches and appends all inputfiles
                dragfiles = glob.glob(arg)
                if (len(dragfiles) == 0):                  #catch Schreibfehler im Pfad
                    print "not a valid file: ", arg
                    sys.exit(-1)                
                else:
                    for i in range(0, len(dragfiles)):
                        print 'using', dragfiles[i], " as drag reference file"
                    dragfiles = merge(dragfiles)
            else:
                print "not a valid file: ", arg, "exiting"
                sys.exit(-1)
            conversion_drag = getReferenceDrag(dragfiles)
            print "conversion drag: ", conversion_drag
            allfiles['drag'] = dragfiles             
        elif opt in ("-l", "--lift-reference"):
            if (os.path.isfile(arg)):                       #checks existence of inputfile
                print "using ", arg, " as lift reference file"
                liftfiles.append(arg)
            elif (len(arg) > 0):                            #catches and appends all inputfiles
                liftfiles = glob.glob(arg)
                if (len(liftfiles) == 0):                  #catch Schreibfehler im Pfad
                    print "not a valid file: ", arg
                    sys.exit(-1)                
                else:
                    for i in range(0, len(liftfiles)):
                        print 'using', liftfiles[i], " as lift reference file"
                    liftfiles = merge(liftfiles)
            else:
                print "not a valid file: ", arg, "exiting"
                sys.exit(-1)
            conversion_lift = getReferenceLift(liftfiles)
            print "conversion lift: ", conversion_lift
            allfiles['lift'] = liftfiles  
        elif opt in ("-v", "--inlet-velocity"):
            print "reference velocity set to: ", arg
            Uc = arg
        elif opt in ("-s", "--reference-area-lift"):
            print "reference area drag set to ", arg
            Ad = arg
        elif opt in ("-a", "--reference-area-lift"):
            print "reference area lift set to ", arg
            Al = arg
            
    data = average(allfiles, 'input')
    zero = average(allfiles, 'zero')
    mount = average(allfiles, 'mount')

#### Mehrere Kalibrierdateien einlesen, sortieren und dann erst verarbeiten
    if  numpy.std(zero.T[0]) == 0:
        zeroLift = numpy.average(zero.T[1])
        zeroDrag = numpy.average(zero.T[3])
        data.T[1] = numpy.subtract(data.T[1], zeroLift)
        data.T[3] = numpy.subtract(data.T[3], zeroDrag)
    else:
        print 'Zero Offset has an inconsistent angle of attack.'

    if numpy.std(mount.T[0]) == 0:
        mountLift = numpy.average(mount.T[1])
        mountDrag = numpy.average(mount.T[3])
        data.T[1] = numpy.subtract(data.T[1], mountLift)
        data.T[3] = numpy.subtract(data.T[3], mountDrag) 
    else:
        print 'Mount Offset has an inconsistent angle of attack.'

    alpha = data[:,ALPHA]

    Cd = mVtoDrag(conversion_drag, data[:,CD],Ad,Uc)
#    Cd_std = mVtoDrag(conversion_drag, data[:,CD_STDV],Ad,Uc)
    Cl = mVtoLift(conversion_lift, data[:,CL], Al, Uc)
#    Cl_std = mVtoLift(conversion_lift, data[:,CL_STDV], Al, Uc)

    plot_fig(alpha, Cd, Cl) # Cd_std, Cl, Cl_std)

if __name__ == "__main__":
    main(sys.argv[1:])
