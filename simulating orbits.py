# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 13:52:35 2021

@author: Erin
"""
import binary_orbit as bo
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

# testing Mike's code for planet simulation
'''
bo.binary_orbit([time, period, semi-major axis, eccentricity, orbital angle,
little omega, inclination], dates compute pos, number times eccentricity
anomoly, do deriv?)
'''
rho, theta, vr = bo.binary_orbit([0,365,1,0,0,0,0],np.array(np.linspace(0,365,1000)))

plt.axis('equal')
plt.plot(rho*np.cos(np.radians(theta)), rho*np.sin(np.radians(theta)))

dd = Table.read('HGCA_vEDR3.fits')

