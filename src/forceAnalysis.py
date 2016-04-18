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

def plot_fig (time, Fpx, Fpy, Fpz, Fvx, Fvy, Fvz):
    ## all the data is now loaded into the RAM
    plt.figure(1)
    plt.plot(time, Fpx)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fpx")
    plt.savefig("plots/fpx.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(2)
    plt.plot(time, Fpy)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fpy")
    plt.savefig("plots/fpy.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(3)
    plt.plot(time, Fpz)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fpz")
    plt.savefig("plots/fpz.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(4)
    plt.plot(time, Fvx)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fvx")
    plt.savefig("plots/fvx.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(5)
    plt.plot(time, Fvy)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fvy")
    plt.savefig("plots/fvy.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(6)
    plt.plot(time, Fvz)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fvz")
    plt.savefig("plots/fvz.png", dpi=300)

def plot_fig_fft (f, Fpxfft, Fpyfft, Fpzfft, Fvxfft, Fvyfft, Fvzfft):
    ## all the data is now loaded into the RAM
    plt.figure(7)
    plt.plot(f, Fpxfft)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fpxfft")
    plt.savefig("plots/fpxfft.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(8)
    plt.plot(f, Fpyfft)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fpyfft")
    plt.savefig("plots/fpyfft.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(9)
    plt.plot(f, Fpzfft)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fpzfft")
    plt.savefig("plots/fpzfft.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(10)
    plt.plot(f, Fvxfft)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fvxfft")
    plt.savefig("plots/fvxfft.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(11)
    plt.plot(f, Fvyfft)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fvyfft")
    plt.savefig("plots/fvyfft.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(12)
    plt.plot(f, Fvzfft)
    plt.grid()
    plt.xlabel("t")
    plt.ylabel("Fvzfft")
    plt.savefig("plots/fvzfft.png", dpi=300)

def fft_analysis (time, Fpx, Fpy, Fpz, Fvx, Fvy, Fvz, Fp_avg, Fv_avg):

    N = len(time)
    T = (max(time) - min(time))/(N - 1)
    f = numpy.linspace(0.0, 1.0/(2*T), N/2)

    hann_filter((Fpx - Fp_avg[0]), len(Fpx))
    hann_filter((Fpy - Fp_avg[1]), len(Fpy))
    hann_filter((Fpz - Fp_avg[2]), len(Fpz))

    hann_filter((Fvx - Fv_avg[0]), len(Fvx))
    hann_filter((Fvy - Fv_avg[1]), len(Fvy))
    hann_filter((Fvz - Fv_avg[2]), len(Fvz))

    Fpxfft_abs = numpy.abs(numpy.fft.fft(Fpx))[0:N/2]
    Fpyfft_abs = numpy.abs(numpy.fft.fft(Fpy))[0:N/2]
    Fpzfft_abs = numpy.abs(numpy.fft.fft(Fpz))[0:N/2]

    Fvxfft_abs = numpy.abs(numpy.fft.fft(Fvx))[0:N/2]
    Fvyfft_abs = numpy.abs(numpy.fft.fft(Fvy))[0:N/2]
    Fvzfft_abs = numpy.abs(numpy.fft.fft(Fvz))[0:N/2]

    N = len(time)
    T = (max(time) - min(time))/(N - 1)

    print "deltaT: ", T

    f = numpy.linspace(0.0, 1.0/(2*T), N/2)

    plot_fig_fft(f, Fpxfft_abs, Fpyfft_abs, Fpzfft_abs, Fvxfft_abs, Fvyfft_abs, Fvzfft_abs)

    maxisFpx = (-Fpxfft_abs[0:N/2]).argsort()[:20]
    print "maxima are: "
    for i in maxisFpx:
        print f[i]

#### START OF THE MAIN LOOP
def main (argv):

    TIME = 0

    FORCE_PRESSURE_X = 1
    FORCE_PRESSURE_Y = 2
    FORCE_PRESSURE_Z = 3

    FORCE_VISCOUS_X  = 4
    FORCE_VISCOUS_Y  = 5
    FORCE_VISCOUS_Z  = 6

    FORCE_POROUS_X   = 7
    FORCE_POROUS_Y   = 8
    FORCE_POROUS_Z   = 9

    TORQUE_PRESSURE_X = 10
    TORQUE_PRESSURE_Y = 11
    TORQUE_PRESSURE_Z = 12

    TORQUE_VISCOUS_X  = 13
    TORQUE_VISCOUS_Y  = 14
    TORQUE_VISCOUS_Z  = 15

    TORQUE_POROUS_X   = 16
    TORQUE_POROUS_Y   = 17
    TORQUE_POROUS_Z   = 18

    inputfile = ''
    outputfile = ''

    starttime = 0.0
    endtime = 10.0
    coeffs = False
    calc_cd = False
    calc_cl = False
    calc_cm = False
    calc_fft = False

    cd_axis = 'x'
    cl_axis = 'y'
    cm_axis = 'z'

    Cd = 0.0
    Cl = 0.0
    Cm = 0.0

    # coeff_area = 0.0051
    coeff_area = 0.0008
    u_c = 1

    # 2 * Fw / rho * u_c^2  * A

    try:
        opts, args = getopt.getopt(argv, "i:o:v:p:f:s:e:", ["input-file=", "output-file=", "inlet-velocity=", "propeller-diameter=", "rotational-frequency=", "start-time=", "end-time=", "fft", "Cd"])
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
            calc_fft = True
        elif opt in ("--coeffs"):
            print "adding coefficients to the workflow"
            calc_coeffs = True
        elif opt in ("--Cd"):
            print "adding drag coefficient calculation to the workflow"
            calc_cd = True
        elif opt in ("--Cl"):
            print "adding lift coefficient calculation to the workflow"
            calc_cl = True
        elif opt in ("-Cm"):
            print "adding moment coefficient calculation to the workflow"
            calc_cm = True



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
