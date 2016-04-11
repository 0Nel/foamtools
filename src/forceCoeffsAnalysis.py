#!/usr/bin/python

import sys, getopt, os, numpy, math
import re
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})

def movingAverage (values, window):
    weights = numpy.repeat(1.0, window)/window
    return numpy.convolve(values, weights, 'valid')

def hann_filter (g, n):
    for i in range(0,len(g)):
        g[i] = g[i] * 0.5 * (1.0 - math.cos((2.0 * math.pi * i) / (n - 1)))

def nextpow2 (i):
    n = 1
    while n < i: n *= 2
    return n

def plot_fig (time, Cd, Cl, Cm, Cl_f, Cl_r):
    ## all the data is now loaded into the RAM
    plt.figure(1)
    plt.plot(time, Cd)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Cd")
    plt.savefig("plots/Cd.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(2)
    plt.plot(time, Cl)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Cl")
    plt.savefig("plots/Cl.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(3)
    plt.plot(time, Cm)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Cm")
    plt.savefig("plots/Cm.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(4)
    plt.plot(time, Cl_f)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Cl(f)")
    plt.savefig("plots/Cl_f.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(5)
    plt.plot(time, Cl_r)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Cl(r)")
    plt.savefig("plots/Cl_r.png", dpi=300)

def plot_fig_fft (f, Cd_fft_abs, Cl_fft_abs, Cm_fft_abs):
    ## all the data is now loaded into the RAM
    plt.figure(7)
    plt.plot(f, Cd_fft_abs)
    plt.axis([0, 20, 0, 500])
    plt.grid()
    plt.xlabel("f")
    plt.ylabel("Amplitude FFT Cd")
    plt.savefig("plots/Cd_fft.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(8)
    plt.plot(f, Cl_fft_abs)
    plt.axis([0, 20, 0, 500])
    plt.grid()
    plt.xlabel("f")
    plt.ylabel("Amplitude FFT Cl")
    plt.savefig("plots/Cl_fft.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(9)
    plt.plot(f, Cm_fft_abs)
    plt.axis([0, 20, 0, 500])
    plt.grid()
    plt.xlabel("f")
    plt.ylabel("Amplitude FFT Cm")
    plt.savefig("plots/Cm_fft.png", dpi=300)

def fft_analysis (time, Cd, Cl, Cm):

    N = len(time)
    T = (max(time) - min(time))/(N - 1)
    f = numpy.linspace(0.0, 1.0/(2*T), N/2)

    hann_filter((Cd - numpy.average(Cd)), len(Cd))
    hann_filter((Cl - numpy.average(Cl)), len(Cl))
    hann_filter((Cm - numpy.average(Cm)), len(Cm))

    Cd_fft_abs = numpy.abs(numpy.fft.fft(Cd))[0:N/2]
    Cl_fft_abs = numpy.abs(numpy.fft.fft(Cl))[0:N/2]
    Cm_fft_abs = numpy.abs(numpy.fft.fft(Cm))[0:N/2]

    print "deltaT: ", T

    f = numpy.linspace(0.0, 1.0/(2*T), N/2)

    plot_fig_fft(f, Cd_fft_abs, Cl_fft_abs, Cm_fft_abs)

    maxCd = (-Cd_fft_abs[0:N/2]).argsort()[:20]
    print "maxima are: "
    for i in maxCd:
        print f[i]

#### START OF THE MAIN LOOP
def main (argv):

    TIME = 0

    CM = 1
    CD = 2
    CL = 3

    CL_f = 4
    CL_r = 5

    inputfile = ''
    outputfile = ''

    starttime = 0.0
    endtime = 10.0

    fft = False

    try:
        opts, args = getopt.getopt(argv, "i:o:v:p:f:s:e:", ["input-file=", "output-file=", "inlet-velocity=", "propeller-diameter=", "rotational-frequency=", "start-time=", "end-time=", "fft"])
    except getopt.GetoptError:
        print 'foamFancy_averageForces.py -i <inputfile> [-o <outputfile> -v velocity -p diameter -f frequency -s starttime -e endtime -f]'
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-i", "--input-file"):
            inputfile = arg
            if (not os.path.isfile(inputfile)):
                print "not a valid input file, exiting..."
                sys.exit(-1)
        elif opt in ("-o", "--output-file"):
            outputfile = arg
        elif opt in ("-v", "--inlet-velocity"):
            Vc = arg
        elif opt in ("-p", "--propeller-diameter"):
            D = arg
        elif opt in ("-f", "--rotational-frequency"):
            F = arg
        elif opt in ("-s", "--start-time"):
            starttime = float(arg)
            print "start time set to: ", starttime
        elif opt in ("-e", "--end-time"):
            print "end time set to: ", endtime
            endtime = float(arg)
        elif opt in ("--fft"):
            print "adding fast fourier transformation analysis"
            fft = True

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

    Cm = numpy.array(raw[:,CM])
    Cd = numpy.array(raw[:,CD])
    Cl = numpy.array(raw[:,CL])

    Cl_f = numpy.array(raw[:,CL_f])
    Cl_r = numpy.array(raw[:,CL_r])

    Cm_avg = numpy.average(Cm)
    Cm_std = numpy.std(Cm)
    
    Cd_avg = numpy.average(Cd)
    Cd_std = numpy.std(Cd)
    
    Cl_avg = numpy.average(Cl)
    Cl_std = numpy.std(Cl)

    Cl_f_avg = numpy.average(Cl_f)
    Cl_r_avg = numpy.average(Cl_r)

    print "Average Drag Coefficient: ", Cd_avg, " +- ", Cd_std
    print "Average Lift Coefficient: ", Cl_avg, " +- ", Cl_std
    print "Average Moment Coefficient: ", Cm_avg, " +- ", Cm_std

    if not os.path.isdir('./plots'):
        os.mkdir('plots')

    f_out = open('plots/summary.txt', 'w')
    f_out.writelines('# Beiwerte:\n'
        + 'Cd: '
        + str(Cd_avg)
        + '\n'
        + 'Cl: '
        + str(Cl_avg)
        + '\n'
        + ' Cm: '
        + str(Cm_avg)
        + '\n' )


    f_out.close()

    plot_fig(time, Cd, Cl, Cm, Cl_f, Cl_r)

    if fft == True:
        fft_analysis(time, Cd, Cl, Cm)


if __name__ == "__main__":
    main(sys.argv[1:])
