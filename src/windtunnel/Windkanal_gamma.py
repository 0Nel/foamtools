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
    # get all the directories that might contain single reynolds number results
	dirs = glob.glob('Results_RE*')
	if len(dirs) > 0:
	        # get permission to overwrite the results files if they already exist
    		check = raw_input("Results directory already exists! Do you want to overwrite? (y/n)\t")
		if check != 'y':
			print 'Quitting job!'
			sys.exit(-1)
		else:
			for dir in dirs:
				print 'Deleting ' + dir + ' Directory'
				shutil.rmtree(dir)  # because of reasons we need to make a syscall

# creates a new directory for the reynolds number results
def mkdir(Re):
	if os.path.exists('Results_RE-'):
		os.rename('Results_RE-', 'Results_RE'+Re)
	else:
		os.makedirs('Results_RE'+Re)

# print the combined results into one file
# REFACTOR
def printRes(alpha, Cd, Cl, Cd_std, Cl_std, Re):		 #writes output data files
    deg = u"\u00b0"
    minCd = str(numpy.min(Cd))
    minCl = str(numpy.min(Cl))
    pos = int(len(alpha)/2)
    if alpha[pos] != 0:
	print 'Number of negative and positive angles is not equal!'
	sys.exit()
    Cl0 = str(Cl[pos])

    #specific values
    f_out = open('Results_RE'+Re+'/min_max_data_raw.txt', 'w')
    f_out.write('Min Cd\t' + minCd + '\n')
    f_out.write('Min Cl\t' + minCl + '\n')
    f_out.write('Cl0 ' + Cl0 + '\n')
    f_out.close()

    #final values to file
    f_out = open('Results_RE'+Re+'/final_averaged_data_Re' + Re + '.txt', 'w')
    f_out.write(Re + '\n')
    f_out.write('alpha\tCl\tCl_Std\tCd\tCd_Std\n')
    for i in range(0, len(alpha)):
	f_out.write(str(alpha[i])+'\t'+str(Cl[i])+'\t'+str(Cl_std[i])+'\t'+str(Cd[i])+'\t'+str(Cd_std[i])+'\n')
    f_out.close()

# creates matplotlib plots with the combined data
def plot_fig (alpha, Cd, Cl, Cd_std, Cl_std, Re):

    deg = u"\u00b0"
    minAngle = numpy.min(alpha)
    maxAngle = numpy.max(alpha)
    minCd = str(numpy.min(Cd))
    minCl = str(numpy.min(Cl))
    maxCl = numpy.max(Cl)
    legend = 'Re' + Re

    plt.figure(3)
    plt.title('Auftriebsbeiwert')
    plt.plot([minAngle, maxAngle], [0, 0], color = 'black', linewidth = 0.5)
    plt.plot(alpha, Cl, 'ro', label = legend, linewidth = 1.0)
    plt.errorbar(alpha, Cl, Cl_std, color='red')
    plt.plot(alpha, Cl, color = 'red', linewidth = 1.0)
    plt.grid()
    plt.xlabel("Anstellwinkel [" + deg + "]")
    plt.ylabel("$C_l$")
    plt.legend(loc = 4)
    plt.savefig("Results_RE"+Re+"/Cl.png", dpi=300)

    plt.figure(4)
    plt.title('Widerstandsbeiwert')
    plt.plot([minAngle, maxAngle], [minCd, minCd], color = 'black', linewidth = 0.5)
    plt.plot(alpha, Cd, 'bo', label = legend, linewidth = 1.0)
    plt.errorbar(alpha, Cd, Cd_std, color='blue')
    plt.plot(alpha, Cd, color = 'blue', linewidth = 1.0)
    plt.grid()
    plt.ylim(0, max(Cd))
    plt.xlabel("Anstellwinkel [" + deg + "]")
    plt.ylabel("$C_d$")
    plt.legend()
    plt.savefig("Results_RE"+Re+"/Cd.png", dpi=300)

    plt.figure(5)
    plt.title('Lilienthalpolare')
    plt.plot(Cd, Cl, 'ko', linewidth = 1.0)
    plt.plot(Cd, Cl, color = 'black', linewidth = 1.0)
    plt.grid()
    plt.xlabel("$C_d$")
    plt.ylabel("$C_l$")
    plt.savefig("Results_RE"+Re+"/lilienthal.png", dpi=300)

### reference function to convert mV/V to N
# REFACTOR to generic regression analysis function
def mVtoDrag(param, Cd, Ac, Uc):

    try:
        Uc = float(Uc)
        Ac = float(Ac)
        return (((param[0] * Cd + param[1])*2*0.009810)/(Ac*Uc**2))
    except(TypeError):
        print('Invalid reference Area or Velocity!')
        sys.exit(-1)

# see commentary above
def mVtoLift(param, Cl, Ac, Uc):
    try:
        Uc = float(Uc)
        Ac = float(Ac)
        return (((param[0] * Cl + param[1])*2*0.009810)/(Ac*Uc**2))
    except(TypeError):
        print('Invalid reference Area or Velocity!')
        sys.exit(-1)

# REFACTOR with get reference lift
def getReferenceDrag(filename, Re):

    weight = []
    force = []

    f_in = open(filename, 'r')
    for line in f_in:
        tmp = line.split()
        if len(tmp) > 0 and not tmp[0] == '#':
            weight.append(float(tmp[0]))
            force.append(float(tmp[1]))
    numpy.subtract(weight, weight[0])

    force2D = numpy.array([force, numpy.ones(len(force))])
    
    param = numpy.linalg.lstsq(force2D.T, weight)[0]

    ## example code for plotting
    plt.figure(1)
    plt.title('linear regression drag')
    plt.plot(force, weight)
    plt.grid()
    plt.xlabel("weight")
    plt.ylabel("force")
    plt.savefig("Results_RE"+Re+"/drag_reference.png", dpi=300)

    return param

def getReferenceLift(filename, Re):

    weight = []
    force = []

    f_in = open(filename, 'r')
    for line in f_in:
        tmp = line.split()
        if not tmp[0] == '#':
            weight.append(float(tmp[0]))
            force.append((tmp[1]))
    numpy.subtract(weight, weight[3])

    force2D = numpy.array([force, numpy.ones(len(force))])


    param = numpy.linalg.lstsq(force2D.T, weight)[0]

    ## example code for plotting
    plt.figure(2)
    plt.title('linear regression lift')
    plt.plot(weight, force)
    plt.grid()
    plt.xlabel("weight")
    plt.ylabel("force")
    plt.savefig("Results_RE"+Re+"/lift_reference.png", dpi=300)

    return param

def toDot (string):
    tmp = string.replace(',', '.')
    return float(tmp)

def parseLine (line):
    try:
        return [toDot(string) for string in line.split()]
    except:
        return None

# deprecated: replace with parseLine
def checkLine(line):
    line = line.replace(',','.')
    tmp = line.split()
    try:
        float(tmp[0])
        return False
    except (ValueError):
        return True

# refactor the parameters -> instead of allfiles the correct file path
def average(allfiles, key, Re):                  #param = Dictionary and Key
    try:
        files = allfiles[key]
        raw2D = []    # leeres array initialisieren
        raw3D = []
        for file in files:
            f_in = open(file, 'r')            
            for line in f_in:
                # new: 
                # raw2d.append(parseLine(line))
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
        f_out = open('Results_RE'+Re+'/'+ key + '_averaged_data_raw.txt', 'w')
    
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
    
def merge(files, Re):
#    if len(files) == 1:
#	print "Only using one reference file"
#	files.append(files[0])
    if len(files) > 2:
        print "Too many reference files given: ", files
    elif len(files[0]) == len(files[1]):		#len ist hier voellig Banane
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
            
                # replace the whole block with parseLine
                
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
        f_out = open('Results_RE'+Re+'/'+name.lower() + '_calibration_averaged_data_raw.txt', 'w')
        for i in range(0, len(raw0)):
            f_out.write(str(t[i]) + ' ' + str(((a[i]+b[i])/2)) +'\n')
        f_out.close()
        
        return (name.lower() + '_calibration_averaged_data_raw.txt')
    else:
          print "Reference files have inconsistent lengths"
          sys.exit(-1)        

#### START OF THE MAIN LOOP
def main (argv):
    print argv

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
    Re = '-'		    #default Reynoldsnumber
    mkdir(Re)		    #default Resultsdirectory
#    Lc = 1.0                # Reference length

    conversion_drag = [1.0, 1.0]
    conversion_lift = [1.0, 1.0]

    # 2 * Fw / rho * u_c^2  * A

    ### now get the command line parameters provided by the user

    try:
        opts, args = getopt.getopt(argv, "i:v:s:a:d:l:z:m:y:", ["input-files=", "inlet-velocity=", "reference-area-drag", "reference-area-lift", "drag-reference", "lift-reference", "zero-offset", "mount-offset", "reynoldsnumber"])
    except getopt.GetoptError:
        print 'windtunnel.py -i <inputfiles> [-v <velocity> -s <reference area drag> -a <reference area lift> -c -d <drag reference file> -l <lift reference file> -z <zero offset> -m <mount offset> -y <reynoldsnumber>]'
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
        elif opt in ("-y", "--reynoldszahl"):
            print "reynoldsnumber set to ", arg
            Re = str(arg)
	    mkdir(Re)
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
                    dragfiles = merge(dragfiles, Re)
            else:
                print "not a valid file: ", arg, "exiting"
                sys.exit(-1)
            conversion_drag = getReferenceDrag("Results_RE"+Re+"/" + dragfiles, Re)
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
                    liftfiles = merge(liftfiles, Re)
            else:
                print "not a valid file: ", arg, "exiting"
                sys.exit(-1)
            conversion_lift = getReferenceLift("Results_RE"+Re+"/" + liftfiles, Re)
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

    data = average(allfiles, 'input', Re)
    zero = average(allfiles, 'zero', Re)
    mount = average(allfiles, 'mount', Re)


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
    Cd_std = mVtoDrag(conversion_drag, data[:,CD_STDV],Ad,Uc)
    Cl = mVtoLift(conversion_lift, data[:,CL], Al, Uc)
    Cl_std = mVtoLift(conversion_lift, data[:,CL_STDV], Al, Uc)

    #Schreibt absolute Minimalwerte und Cl-Wert bei Anstellwinkel null


    printRes(alpha, Cd, Cl, Cd_std, Cl_std, Re)

    plot_fig(alpha, Cd, Cl, Cd_std, Cl_std, Re)

if __name__ == "__main__":
    checkdir()
    main(sys.argv[1:])
