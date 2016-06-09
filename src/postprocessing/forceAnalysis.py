#!/usr/bin/python

import sys, getopt, os, numpy, math
import re
import matplotlib.pyplot as plt

plt.rcParams.update({	'legend.fontsize' : 18,
		'axes.titlesize' : 20, 	#Title
		'axes.labelsize' : 18	#Axenbeschriftung
	  })

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

def plot_Cx (time, Cd, Cl):

    plt.figure(13)
    plt.grid()
    plt.xlabel("$t^{*}$")
    plt.ylabel("$C_{D}$")
    # plt.title("Druckkraft in $X$-Richtung")
    plt.title("Zeitlicher Verlauf des Widerstandsbeiwerts")
    plt.plot(time, Cd, linestyle = 'dashed')
    plt.savefig("plots/cd.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(14)
    plt.grid()
    plt.xlabel("$t^{*}$")
    plt.ylabel("$C_{L}$")
    plt.title("Zeitlicher Verlauf des Auftriebsbeiwerts")
    # plt.title("Druckkraft in $Y$-Richtung")
    plt.plot(time, Cl, linestyle = 'dashed')
    plt.savefig("plots/cl.png", dpi=300)


def plot_fig (time, Fpx, Fpy, Fpz, Fvx, Fvy, Fvz):

    minTime = numpy.min(time)
    maxTime = numpy.max(time)

    ## all the data is now loaded into the RAM
    plt.figure(1)
    plt.grid()
    plt.xlabel("$t [s]$")
    plt.title("Druckkraft in $X$-Richtung")
    plt.title("Druckkraft in X-Richtung")
    plt.plot(time, Fpx, linestyle = 'dashed')
    plt.savefig("plots/fpx.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(2)
    plt.grid()
    plt.xlabel("$t [s]$")
    plt.ylabel("$F_{p,y} [N]$")
    plt.title("Druckkraft in $Y$-Richtung")
    plt.plot(time, Fpy, linestyle = 'dashed')
    plt.savefig("plots/fpy.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(3)
    plt.grid()
    plt.xlabel("$t [s]$")
    plt.ylabel("$F_{p,z} [N]$")
    plt.title("Druckkraft in $Z$-Richtung")
    plt.plot(time, Fpz, linestyle = 'dashed')
    plt.savefig("plots/fpz.png", dpi=300)

'''
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
'''

def plot_fig_fft (f, Fpxfft, Fpyfft, Fpzfft, Fvxfft, Fvyfft, Fvzfft, fmax):
    ## all the data is now loaded into the RAM

    # voodo: test the average of amplitudes against the max. if >> than factor 100, make it a logarithmic plot

    Fpxfft[0] = 0.0
    Fpyfft[0] = 0.0
    Fpzfft[0] = 0.0

    plt.figure(7)
    if numpy.max(Fpxfft) / 100.0 - numpy.average(Fpxfft) > 100:
        plt.loglog(f, Fpxfft)
    else:
        plt.semilogx(f, Fpxfft)
    plt.grid()
    plt.title("Spektrum $F_{px}$")
    plt.xlabel("$f [ s^{-1} ]$")
    plt.ylabel("$A_{F_{p,x},fft}$")
    plt.savefig("plots/fpxfft.png", dpi=300)

    print f[numpy.argmax(Fpxfft)] 

    ## all the data is now loaded into the RAM
    plt.figure(8)
    if numpy.max(Fpyfft) / 100.0 - numpy.average(Fpyfft) > 100:
        plt.loglog(f, Fpxfft)
    else:
        plt.semilogx(f, Fpxfft)
    plt.grid()
    plt.title("Spektrum $F_{py}$")
    plt.xlabel("$f  [ s^{-1} ]$")
    plt.ylabel("$A_{F_{p,y},fft}$")
    plt.savefig("plots/fpyfft.png", dpi=300)

    ## all the data is now loaded into the RAM
    plt.figure(9)
    # plt.xlim(0, fmax)
    plt.semilogx(f, Fpzfft)
    plt.grid()
    plt.title("Spektrum $F_{pz}$")
    plt.xlabel("$f  [ s^{-1} ]$")
    plt.ylabel("$A_{F_{p,z},fft}$")
    plt.savefig("plots/fpzfft.png", dpi=300)

'''
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
'''

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

    plot_fig_fft(f, Fpxfft_abs, Fpyfft_abs, Fpzfft_abs, Fvxfft_abs, Fvyfft_abs, Fvzfft_abs, 1000.0)

#    maxisFpx = (-Fpxfft_abs[0:N/2]).argsort()[:20]
#    print "maxima are: "
#    for i in maxisFpx:
#        print f[i]

def toDot (string):
    tmp = string.replace(',', '.')
    return float(tmp)

def parseLine (line):
    try:
        return [toDot(string.strip('(').strip(')')) for string in line.split()]
    except:
        return None

# just coz itz coolz
def argmax(L):
    return sorted((item, index) for index, item in enumerate(L))[-1][1]

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
    density = 1.0
    coeffs = False
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
    l_c = 1

    # 2 * Fw / rho * u_c^2  * A

    try:
        opts, args = getopt.getopt(argv, "i:o:v:a:s:e:d:", ["input-file=", "output-file=", "inlet-velocity=", "reference-area", "start-time=", "end-time=", "density", "fft", "coeffs"])
    except getopt.GetoptError:
        print 'foamFancy_averageForces.py -i <inputfile> [-o <outputfile> -v velocity -s starttime -e endtime -a Area -d density --Coeffs --fft]'
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
            u_c = float(arg)
        elif opt in ("-a", "--reference-area"):
            coeff_area = float(arg)
        elif opt in ("-s", "--start-time"):
            starttime = float(arg)
            print "start time set to: ", starttime
        elif opt in ("-e", "--end-time"):
            print "end time set to: ", endtime
            endtime = float(arg)
        elif opt in ("-d", "--density"):
            print "correcting for density:", arg
            density = float(arg)
        elif opt in ("--fft"):
            print "adding fast fourier transformation analysis"
            calc_fft = True
        elif opt in ("--coeffs"):
            print "adding coefficients to the workflow"
            coeffs = True


    raw = []    # leeres array initialisieren

    f_in = open(inputfile, 'r')

    for line in f_in:
        tmp = parseLine(line)
        if tmp == None:
            continue
        if tmp[TIME] < starttime:
            continue
        elif tmp[TIME] > endtime:
            break
        else:
            raw.append(tmp)

    raw = numpy.array(raw)

    time = numpy.array(raw[:,TIME])

    Fpx = numpy.array(raw[:,FORCE_PRESSURE_X] * density)
    Fpy = numpy.array(raw[:,FORCE_PRESSURE_Y] * density)
    Fpz = numpy.array(raw[:,FORCE_PRESSURE_Z] * density)

    # do we need a density correction for the forces due to wall shear stress?
    Fvx = numpy.array(raw[:,FORCE_VISCOUS_X] * density)
    Fvy = numpy.array(raw[:,FORCE_VISCOUS_Y] * density)
    Fvz = numpy.array(raw[:,FORCE_VISCOUS_Z] * density)

    Fp_avg = [ numpy.average(Fpx), numpy.average(Fpy), numpy.average(Fpz) ]
    Fp_stdv = [ numpy.std(Fpx), numpy.std(Fpy), numpy.std(Fpz) ]

    Fv_avg = [ numpy.average(Fvx), numpy.average(Fvy), numpy.average(Fvz) ]
    Fp_stdv = [ numpy.std(Fvx), numpy.std(Fvy), numpy.std(Fvz) ]

    print "Average Forces due to Pressure: ", Fp_avg
    print "Average Froces due to Viscous Effect: ", Fv_avg

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

    if coeffs == True:
        Cd = (2 * (Fp_avg[0] + Fv_avg[0])) / ( density * u_c * u_c * coeff_area)
        print "Average Drag Coefficient: ", Cd
        f_out.writelines('Average Drag Coefficient: ' + str(Cd) + '\n')
        Cl = (2 * (Fp_avg[1] + Fv_avg[1])) / ( density * u_c * u_c * coeff_area)
        print "Average Lift Coefficient: ", Cl
        f_out.writelines('Average Lift Coefficient: ' + str(Cl) + '\n')

    f_out.close()

    plot_fig(time, Fpx, Fpy, Fpz, Fvx, Fvy, Fvz)

    if calc_fft == True:
        fft_analysis(time, Fpx, Fpy, Fpz, Fvx, Fvy, Fvz, Fp_avg, Fv_avg)

    if coeffs == True:
        plot_Cx( (time * u_c * l_c) , ((2 * (Fpx + Fvx)) / ( density * u_c * u_c * coeff_area)), ((2 * (Fpy + Fvy)) / ( density * u_c * u_c * coeff_area)) )

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
