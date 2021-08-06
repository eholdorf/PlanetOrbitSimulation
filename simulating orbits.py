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
import re
import astropy.units as u


# testing Mike's code for planet simulation
'''
bo.binary_orbit([time, period, semi-major axis, eccentricity, orbital angle,
little omega, inclination], dates compute pos, number times eccentricity
anomoly, do deriv?)
'''
rho, theta, vr = bo.binary_orbit([2021,365,1,0,0,0,0],np.array([1991.25,2015,2017]))

plt.figure()
plt.clf()
plt.axis('equal')
plt.plot(rho*np.cos(np.radians(theta)), rho*np.sin(np.radians(theta)))

# importing the data
data = Table.read('HGCA_vEDR3.fits')
data_in_30pc = [data[i] for i in np.where((data['parallax_gaia']>33))[0]]
data_30pc_chisq = [data[i] for i in np.where((data['parallax_gaia']>33) & 
                            (data['chisq']>16) & (data['chisq']<100))[0]]

bp_rp = []
f = open('bp_rps.txt','r')
for line in f.readlines():
    bp_rp.append(float(line))   
f.close()

a = models.interpolate_semimajoraxis()
m = models.interpolate_mass(bp_rp, len(bp_rp))
T = models.interpolate_period(a,m)
omega = models.interpolate_omega()
Omega = models.interpolate_Omega()
i = models.interpolate_inclination()
e = models.interpolate_eccentricity()

def simulate_orbit(time= 0, T=0, a=0, e=0, Omega=0, omega=0, i=0, dates=0, m=0, 
                   m_p = 0, random = False, plot = False):
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
    plot : TYPE, optional
        DESCRIPTION. The defult is False.

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
        rho_planet, theta, vr = bo.binary_orbit([time,T,a,e,Omega,omega,i],dates)
    else:
        j = np.random.randint(0,10000)
        
        a = models.interpolate_semimajoraxis()[j]
        T = models.interpolate_period([a],[1])[0]
        omega = models.interpolate_omega()[np.random.randint(0,10000)]
        Omega = models.interpolate_Omega()[np.random.randint(0,10000)]
        i = models.interpolate_inclination()[np.random.randint(0,10000)]
        e = models.interpolate_eccentricity()[np.random.randint(0,10000)]
        
        rho_planet, theta, vr = bo.binary_orbit([0,T,a,e,Omega,omega,i],
                                        np.array(np.linspace(0,T,1000)))
        rho_star = rho_planet * m/m_p
        dec = rho_star*np.cos(np.radians(theta))
        ra = rho_star*np.sin(np.radians(theta))
        
    if plot:    
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
        
    rho_star = rho_planet * m_p/m
    dec = rho_star*np.cos(np.radians(theta))
    ra = rho_star*np.sin(np.radians(theta))
    return [[rho_star,rho_planet], [theta, ra, dec], vr,[time,T,a,e,Omega,omega,i]]

# trial plot simulation with first star
all_data = []
for val in bp_rp[0:2]:
    dats = []
    for i in range(2): 
        m_planet = models.interpolate_planet_mass()[np.random.randint(0,10000)]
        m_star = models.interpolate_mass(val,1)
        t_k = np.array([1991.25,2015.0,2017.0])
        a = models.interpolate_semimajoraxis()[np.random.randint(0,10000)]
        T = models.interpolate_period([a],[m_star])[0]
        omega = models.interpolate_omega()[np.random.randint(0,10000)]
        Omega = models.interpolate_Omega()[np.random.randint(0,10000)]
        i = models.interpolate_inclination()[np.random.randint(0,10000)]
        e = models.interpolate_eccentricity()[np.random.randint(0,10000)]
        
        dats.append(simulate_orbit(2021,T,a,e,Omega,omega,i,t_k,m_star,m_planet))
    all_data.append(dats)
    
def proper_motion(RA1,DEC1,t1,RA2,DEC2,t2):
    RA1 = (RA1 * u.deg).to(u.mas)
    RA2 = (RA2 * u.deg).to(u.mas)
    DEC1 = (DEC1 * u.deg).to(u.mas)
    DEC2 = (DEC2 * u.deg).to(u.mas)
    delta_RA = RA2 - RA1
    delta_DEC = DEC2 - DEC2
    delta_t = round((t2 - t1)*365.25)*u.yr
    # ask mike which DEC to use in cos
    mu_RA = delta_RA/delta_t # * np.cos(np.radians(DEC1.to(u.deg)))
    mu_DEC = delta_DEC/delta_t
    
    mu = np.sqrt(mu_RA**2 - mu_DEC**2)
    
    return [mu, mu_RA, mu_DEC]

def chisq(mu_S1,mu_S2,mu_E1,mu_E2):
    
    mu_S = abs(abs(mu_S1) - abs(mu_S2))
    mu_E = abs(abs(mu_E1) - abs(mu_E2))
    
    return (mu_S - mu_E)**2 /mu_E

def calc_chi_sq(data,mu_RA_hip, mu_RA_gaia, mu_DEC_hip,mu_DEC_gaia,i,j):
    
    # i is the star iterator
    # j is the simulation iterator
    
    # simulation values of RA and DEC between gaia dates
    a = proper_motion(data[i][j][1][1][1],data[0][0][1][2][1],2015.0,
                      data[i][j][1][1][2],data[0][0][1][2][2],2017.0)

    # simulation values of RA and DEC between hipparcos and gaia
    b = proper_motion(data[i][j][1][1][0],data[0][0][1][2][0],1991.25,
                      data[i][j][1][1][2],all_data[0][0][1][2][2],2017.0)

    chi_RA = chisq(b[1].value,a[1].value,mu_RA_hip,mu_RA_gaia)
    chi_DEC = chisq(b[2].value,a[2].value,mu_DEC_hip,mu_DEC_gaia)
    
    return [chi_RA,chi_DEC]