# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 12:35:36 2015

@author: asander
"""

import sys, numpy, math
import matplotlib.pyplot as plt

def hann_filter (g, n):
    for i in range(0,len(g)):  
        g[i] = g[i] * 0.5 * (1.0 - math.cos((2.0 * math.pi * i) / (n - 1)))
        
def main (argv):

    N = 600.0
    
    T = 1.0/800.0
    
    x = numpy.linspace(0.0, N*T, N)
    
    y = numpy.sin(50.0 * 2.0 * numpy.pi * x) + 0.5 * numpy.sin(80.0 * 2.0 * numpy.pi * x)
    
    yf = numpy.fft.fft(y)
    
    xf = numpy.linspace(0.0, 1.0/(2.0*T), N/2)
    
    print "length yf: ", len(yf), " length xf: ", len(xf)
    
    plt.plot(xf, 2.0/N * numpy.abs(yf[0:N/2]))
    
    plt.grid()
    
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])    
    