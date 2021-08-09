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

def simulate_orbit(time= 0, T=0, a=0, e=0, Omega=0, omega=0, i=0, dates=0, m=1, 
                   m_p = 0.0029, parallax=100, random = False, plot = False):
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
    The RA and Dec arrays for the star, in milli-arcsec from the center of mass. 
    NB if you want to save the input data, do this separately.
    '''
    if random == False and a == 0:
        raise ValueError('No values given and not a random orbit.')
     
    if random == False:
        rho_planet, theta, vr = bo.binary_orbit([time,T,a,e,Omega,omega,i],dates)
    else:
        
        m_planet = models.interpolate_planet_mass(num=1)[0]
        m_star = models.interpolate_mass(this_bprp,ddir=data_dir)
        t_k = np.array([1991.25,2015.0,2017.0])
        a = models.interpolate_semimajoraxis(num=1)[0]
        T = models.interpolate_period([a],[m_star])[0]
        omega = models.interpolate_omega(num=1)[0]
        Omega = models.interpolate_Omega(num=1)[0]
        i = models.interpolate_inclination(num=1)[0]
        e = models.interpolate_eccentricity(num=1)[0]
        
        rho_planet, theta, vr = bo.binary_orbit([0,T,a,e,Omega,omega,i],t_k)
        rho_star = rho_planet * (m_planet/m_star)/parallax
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
        
    #Convert planet separation in au to star separation from center of mass
    #in milli-arcsec.
    rho_star = rho_planet * (m_p/m) * parallax #ask mike
    # rho_star = rho_planet * (m_p/m) / parallax
    dec = rho_star*np.cos(np.radians(theta))
    ra = rho_star*np.sin(np.radians(theta))
    #return [[rho_star,rho_planet], [theta, ra, dec], vr,[time,T,a,e,Omega,omega,i]]
    return np.array([ra,dec])
    
def proper_motion(RA1,DEC1,t1,RA2,DEC2,t2):
    #Returns proper motion in mas/year
    #Inputs must be in mas and years!
    delta_RA = RA2 - RA1
    delta_DEC = DEC2 - DEC1
    delta_t = round((t2 - t1))
    # ask mike which DEC to use in cos
    mu_RA = delta_RA/delta_t # * np.cos(np.radians(DEC1.to(u.deg)))
    mu_DEC = delta_DEC/delta_t
    
    #mu = np.sqrt(mu_RA**2 - mu_DEC**2)
    
    return [mu_RA, mu_DEC]


def calc_chi_sq(sim,data_row,t_k):
    """
    Asks the question: what would the chi-squared have been for this planet if it had
    the same uncertainties as the Gaia and HIPPARCOS data?
    """
    
    # simulation values of RA and DEC between gaia dates
    sim_pm_gaia = proper_motion(sim[0,1], sim[1,1], t_k[1], sim[0,2], sim[1,2], 
                                t_k[2])
    
    # simulation values of RA and DEC between hipparcos and gaia
    sim_pm_hg = proper_motion(sim[0,0], sim[1,0], t_k[0], sim[0,2], sim[1,2], 
                              t_k[2])

    chi2 = 0
    chi2 += (sim_pm_gaia[0]-sim_pm_hg[0])**2/(data_row['pmra_hg_error']**2 +
                                              data_row['pmra_gaia_error']**2)
    chi2 += (sim_pm_gaia[1]-sim_pm_hg[1])**2/(data_row['pmdec_hg_error']**2 + 
                                              data_row['pmdec_gaia_error']**2)
    
    return chi2

if __name__=="__main__":
    #For testing, you may not want all stars!
    test_binary_orbit=False
    data_dir = './' #For this directory, use './'
    # data_dir = '/Users/mireland/Google Drive/EDR3_Planets/'
    n_sims = 1000
    n_stars = 1 #Use -1 for all stars.
    
    if test_binary_orbit:
        rho, theta, vr = bo.binary_orbit([2021,1,1,0,0,0,0],np.array([1991.25,2015,2017]))

        plt.figure()
        plt.clf()
        plt.axis('equal')
        plt.plot(rho*np.cos(np.radians(theta)), rho*np.sin(np.radians(theta)))

    # Import the data
    data = Table.read(data_dir + 'HGCA_vEDR3.fits')
    data_in_30pc = data[np.where((data['parallax_gaia']>33))[0]]
    data_30pc_chisq = data[np.where((data['parallax_gaia']>33) & 
                                (data['chisq']>16) & (data['chisq']<100))[0]]
    
    # here need to check that want to be using the right data
    bp_rp = []
    f = open(data_dir + 'bp_rps.txt','r')
    for line in f.readlines():
        bp_rp.append(float(line))   
    f.close()
    if n_stars==-1:
        n_stars = len(bp_rp)

    #a = models.interpolate_semimajoraxis()
    #m = models.interpolate_mass(bp_rp, len(bp_rp))
    #T = models.interpolate_period(a,m)
    #omega = models.interpolate_omega()
    #Omega = models.interpolate_Omega()
    #i = models.interpolate_inclination()
    #e = models.interpolate_eccentricity()

    # trial plot simulation with first star
    all_star_radecs = []
    for this_bprp, parallax in zip(bp_rp[:n_stars], 
                                   data_30pc_chisq[:n_stars]['parallax_gaia']):
        radecs = []
        for i in range(n_sims): 
            m_planet = models.interpolate_planet_mass(num=1)[0]
            m_star = models.interpolate_mass(this_bprp,ddir=data_dir)
            t_k = np.array([1991.25,2015.0,2017.0])
            a = models.interpolate_semimajoraxis(num=1)[0]
            T = models.interpolate_period([a],[m_star+m_planet])[0]
            omega = models.interpolate_omega(num=1)[0]
            Omega = models.interpolate_Omega(num=1)[0]
            i = models.interpolate_inclination(num=1)[0]
            e = models.interpolate_eccentricity(num=1)[0]
        
            radecs.append(simulate_orbit(2021,T,a,e,Omega,omega,i,t_k,m_star,
                                         m_planet,parallax))
        all_star_radecs.append(radecs)

    chis = np.empty((n_stars, n_sims))
    for i in range(n_stars):
        for j in range(n_sims):
            chis[i,j] = calc_chi_sq(radecs[j],data_30pc_chisq[i],t_k)
            
      #Should search for "significant" chi-squared values     
      #  for elem in chi:
      #      if elem > 9.5:
      #          print(elem, j, data_30pc_chisq[i])