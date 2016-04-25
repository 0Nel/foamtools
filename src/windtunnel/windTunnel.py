#!/usr/bin/python

import sys, getopt, os, numpy, math, glob
import re
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})

def nextpow2 (i):
    n = 1
    while n < i: n *= 2
    return n

# def plot_fig (alpha, Cl, Cl_stdv, Cd, Cd_stdv):

    ### plot an Lift polar, drag polar and a lilienthalpolar

    ## example code for plotting
    # plt.figure(1)
    # plt.plot(time, Fpx)
    # plt.grid()
    # plt.xlabel("t")
    # plt.ylabel("Fpx")
    # plt.savefig("plots/fpx.png", dpi=300)


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

    coeffs = False          # calculate coeffs?

    cd_axis = 'x'
    cl_axis = 'y'
    cm_axis = 'z'

    Ac = 1.0                # Reference area, default value
    Uc = 1.0                # Reference speed
    Lc = 1.0                # Reference length

    # 2 * Fw / rho * u_c^2  * A

    ### now get the command line parameters provided by the user

    try:
        opts, args = getopt.getopt(argv, "i:o:v:p:f:s:e:", ["input-files=", "output-file=", "inlet-velocity=", "reference-area", "coeffs"])
    except getopt.GetoptError:
        print 'foamFancy_averageForces.py -i <inputfiles> [-o <outputfile> -v <velocity> -r <reference area> -c]'
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
            Uc = arg
        elif opt in ("-c", "--coeffs"):
            print "adding coefficients to the workflow"
            calc_coeffs = True
        elif opt in ("-r", "--reference-area"):
            print "setting reference area to ", arg
            Ac = arg

    sys.exit(0)

    raw = []    # leeres array initialisieren

    f_in = open(inputfile, 'r')

    for line in f_in:
        tmp = [x.strip('(').strip(')') for x in line.split()]
        if tmp[0] == '#':           # skip comments
            next
        else:
            if float(tmp[0]) < starttime:
                next
            elif float(tmp[0]) > endtime:
                break
            else:
                raw.append(map(float, tmp))

    raw = numpy.array(raw)

    time = numpy.array(raw[:,TIME])

    Fpx = numpy.array(raw[:,FORCE_PRESSURE_X])
    Fpy = numpy.array(raw[:,FORCE_PRESSURE_Y])
    Fpz = numpy.array(raw[:,FORCE_PRESSURE_Z])

    Fvx = numpy.array(raw[:,FORCE_VISCOUS_X])
    Fvy = numpy.array(raw[:,FORCE_VISCOUS_Y])
    Fvz = numpy.array(raw[:,FORCE_VISCOUS_Z])

    Fp_avg = [ numpy.average(Fpx), numpy.average(Fpy), numpy.average(Fpz) ]
    Fp_stdv = [ numpy.std(Fpx), numpy.std(Fpy), numpy.std(Fpz) ]

    Fv_avg = [ numpy.average(Fvx), numpy.average(Fvy), numpy.average(Fvz) ]
    Fp_stdv = [ numpy.std(Fvx), numpy.std(Fvy), numpy.std(Fvz) ]

    print "Average Forces due to Pressure: ", Fp_avg
    print "Average Froces due to Viscous Effect: ", Fv_avg

    if calc_cd == True:
        Cd = (2 * Fp_avg[0]) / ( u_c * u_c * coeff_area)
        print "Average Drag Coefficient: ", Cd

    if not os.path.isdir('./plots'):
        os.mkdir('plots')

    f_out = open('plots/summary.txt', 'w')
    f_out.writelines('# Average Forces due to pressure:\n'
        +'Fpx: '
        + str(Fp_avg[0])
        +  ' Fpy: '
        + str(Fp_avg[1])
        + ' Fpz: '
        + str(Fp_avg[2])
        + '\n' )

    f_out.writelines('# Average Forces due to viscous effects:\n'
        + 'Fpx: '
        + str(Fv_avg[0])
        +  ' Fpy: '
        + str(Fv_avg[1])
        + ' Fpz: '
        + str(Fv_avg[2])
        + '\n' )

    if (calc_cd == True):
        f_out.writelines('# Average Drag Coefficient:\n'
            + 'Cd: '
            + str(Cd)
            + '\n')

    f_out.close()

    plot_fig(time, Fpx, Fpy, Fpz, Fvx, Fvy, Fvz)

    if calc_fft == True:
        fft_analysis(time, Fpx, Fpy, Fpz, Fvx, Fvy, Fvz, Fp_avg, Fv_avg)

    # nach Dr. Baars
    # Nfft = int(math.pow(2, math.ceil(math.log(len(Fx),2))))

    # Fx_fft = numpy.fft.fft(Fx, Nfft)
    # Fx_fft_freq = numpy.fft.fftfreq(Nfft, d=sf)

    #Fx_fft_freq = fs/2 * numpy.linspace(0, 1, Nfft)

    #Fx_fft.abs

    # Fx_fft_abs = abs(Fx_fft)
    #plt.plot(Fx_fft_freq, Fx_fft_abs, 'rs')
    # plt.plot(Fx_fft_freq, Fx_fft_abs)
    #plt.hist(Fx_fft_abs)
    # plt.xlim(0)
    # plt.yscale('log')

if __name__ == "__main__":
    main(sys.argv[1:])
