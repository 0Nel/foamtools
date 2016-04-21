#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 09:22:09 2013
@author: pete
"""

import matplotlib.pyplot as plt
import re
import numpy as np
import sys

forceRegex = r"([0-9.eE\-+]+)\s+\(+([0-9.eE\-+]+)\s([0-9.eE\-+]+)\s([0-9.eE\-+]+)\)"
forceRegex += r"\,\(([0-9.eE\-+]+)\s([0-9.eE\-+]+)\s([0-9.eE\-+]+)\)"
forceRegex += r"\,\(([0-9.eE\-+]+)\s([0-9.eE\-+]+)\s([0-9.eE\-+]+)\)+"
forceRegex += r"\s+\(+([0-9.eE\-+]+)\s([0-9.eE\-+]+)\s([0-9.eE\-+]+)\)"
forceRegex += r"\,\(([0-9.eE\-+]+)\s([0-9.eE\-+]+)\s([0-9.eE\-+]+)\)"
forceRegex += r"\,\(([0-9.eE\-+]+)\s([0-9.eE\-+]+)\s([0-9.eE\-+]+)\)+"

t = []
fpx = []; fpy = []; fpz = []
fpox = []; fpoy = []; fpoz = []
fvx = []; fvy = []; fvz = []
mpx = []; mpy = []; mpz = []
mpox = []; mpoy = []; mpoz = []
mvx = []; mvy = []; mvz = []

pipefile = open('./forces.dat','r')
lines = pipefile.readlines()

for line in lines:
        match = re.search(forceRegex,line)
        if match:
                print "match"
                t.append(float(match.group(1)))
                fpx.append(float(match.group(2)))
                fpy.append(float(match.group(3)))
                fpz.append(float(match.group(4)))
                fvx.append(float(match.group(5)))
                fvy.append(float(match.group(6)))
                fvz.append(float(match.group(7)))
                fpox.append(float(match.group(8)))
                fpoy.append(float(match.group(9)))
                fpoz.append(float(match.group(10)))
                mpx.append(float(match.group(11)))
                mpy.append(float(match.group(12)))
                mpz.append(float(match.group(13)))
                mvx.append(float(match.group(14)))
                mvy.append(float(match.group(15)))
                mvz.append(float(match.group(16)))
                mpox.append(float(match.group(17)))
                mpoy.append(float(match.group(18)))
                mpoz.append(float(match.group(19)))

print len(t)
# Convert to numpy arrays
t = np.asarray(t)
torque = np.asarray(np.asarray(mpz) + np.asarray(mvz))

sys.exit(0)

# Import turbine angular velocity in degrees per second
f = open("constant/dynamicMeshDict", "r")
dpsRegex = r"\(([0-9.eE\-+]+)\s([0-9.eE\-+]+)\s([0-9.eE\-+]+)\)"
for line in f.readlines():
    if "radialVelocity" in line:
        dps = float(re.search(dpsRegex, line).group(3))
f.close()

# Create a theta vector
theta = t*dps

# Compute tip speed ratio
r = 0.5
U = 1.0
omega = dps/360*2*np.pi
tsr = omega*r/U
print "Mean tsr:", tsr

# Pick an index to start from for mean calculations and plotting
# (allow turbine to reach steady state)
i = 4000

# Compute power coefficient
area = 1.0*0.05
power = torque*omega
rho = 1000.0
cp = power/(0.5*rho*area*U**3)
print "Mean cp:", np.mean(cp[i:])

# Compute drag coefficient
drag = np.asarray(np.asarray(fpx) + np.asarray(fvx))
cd = drag/(0.5*rho*area*U**2)
print "Mean cd:", np.mean(cd[i:])

plt.close('all')
plt.plot(theta[i:], torque[i:])
plt.title(r"Torque at $\lambda = %1.1f$" %tsr)
plt.xlabel(r"$\theta$ (degrees)")
plt.ylabel("Torque (Nm)")
plt.show()
