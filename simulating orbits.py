# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 13:52:35 2021

@author: Erin
"""
import binary_orbit as bo
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import model_parameters as models
from mpl_toolkits import mplot3d

# testing Mike's code for planet simulation
'''
bo.binary_orbit([time, period, semi-major axis, eccentricity, orbital angle,
little omega, inclination], dates compute pos, number times eccentricity
anomoly, do deriv?)
'''
rho, theta, vr = bo.binary_orbit([0,365,1,0,0,0,0],np.array(np.linspace(0,365,1000)))

plt.figure()
plt.clf()
plt.axis('equal')
plt.plot(rho*np.cos(np.radians(theta)), rho*np.sin(np.radians(theta)))

# importing the data
data = Table.read('HGCA_vEDR3.fits')
data_in_30pc = [data[i] for i in np.where((data['parallax_gaia']>33))[0]]

a = models.interpolate_semimajoraxis()
T = models.interpolate_period(a)
omega = models.interpolate_omega()
Omega = models.interpolate_Omega()
i = models.interpolate_inclination()
e = models.interpolate_eccentricity()

def simulate_orbit(time= 0, T=0, a=0, e=0, Omega=0, omega=0, i=0, dates=0, 
                   random = False):
    '''

    Parameters
    ----------
    time : TYPE, optional
        DESCRIPTION. The default is 0.
    T : TYPE, optional
        DESCRIPTION. The default is 0.
    a : TYPE, optional
        DESCRIPTION. The default is 0.
    e : TYPE, optional
        DESCRIPTION. The default is 0.
    Omega : TYPE, optional
        DESCRIPTION. The default is 0.
    omega : TYPE, optional
        DESCRIPTION. The default is 0.
    i : TYPE, optional
        DESCRIPTION. The default is 0.
    dates : TYPE, optional
        DESCRIPTION. The default is 0.
    random : TYPE, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    if random == False and a == 0:
        raise ValueError('No values given and not a random orbit.')
     
    if random == False:
        rho, theta, vr = bo.binary_orbit([time,T,a,e,Omega,omega,i],dates)
    else:
        j = np.random.randint(0,10000)
        
        a = models.interpolate_semimajoraxis()[j]
        T = models.interpolate_period([a])[0]
        omega = models.interpolate_omega()[np.random.randint(0,10000)]
        Omega = models.interpolate_Omega()[np.random.randint(0,10000)]
        i = models.interpolate_inclination()[np.random.randint(0,10000)]
        e = models.interpolate_eccentricity()[np.random.randint(0,10000)]
        
        rho, theta, vr = bo.binary_orbit([0,T,a,e,Omega,omega,i],
                                        np.array(np.linspace(0,T,1000)))
        
    print('The time of periastron passage:',time)
    print('Period:', T, 'days.')
    print('Semi-major axis:', a,'au.')
    print('Eccentricity:',e)
    print('Omega:', Omega, 'degrees.')
    print('omega:', omega, 'degrees.')
    print('Inclination:', i, 'degrees.')
        
    plt.figure()
    plt.clf()
    plt.axis('equal')
    plt.plot(rho*np.cos(np.radians(theta)), rho*np.sin(np.radians(theta)),'k')
    
    plt.figure()
    plt.clf()
    ax = plt.axes(projection='3d')
    ax.plot(rho*np.cos(np.radians(theta)), rho*np.sin(np.radians(theta)),'k')
    ax.plot([0],[0],[0],'o', color = 'orange')
    
    