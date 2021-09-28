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
import astropy.units as u
plt.ion()

# testing Mike's code for planet simulation
'''
bo.binary_orbit([time, period, semi-major axis, eccentricity, orbital angle,
little omega, inclination], dates compute pos, number times eccentricity
anomoly, do deriv?)
'''

def simulate_orbit(time= 0, T=0, a=0, e=0, Omega=0, omega=0, i=0, dates=0, m=0, 
                   m_p = 0, parallax=100, random = False, plot = False):
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
    The RA and Dec arrays for the star, in milli-arcsec from center of mass
    NB if you want to save the input data, do this separately.
    '''
    if random == False and a == 0:
        raise ValueError('No values given and not a random orbit.')
     
    if random == False:
        rho_planet, theta, vr = bo.binary_orbit(
            [time,T,a,e,Omega,omega,i],dates)
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
        plt.plot(rho*np.cos(np.radians(theta)),rho*np.sin(np.radians(theta)),
                 'k')
        
        plt.figure()
        plt.clf()
        ax = plt.axes(projection='3d')
        ax.plot(rho*np.cos(np.radians(theta)),rho*np.sin(np.radians(theta)),
                'k')
        ax.plot([0],[0],[0],'o', color = 'orange')
        
    #Convert planet separation in au to star separation from center of mass
    #in milli-arcsec.
    rho_star = rho_planet * (m_p/m) * parallax
    dec = rho_star*np.cos(np.radians(theta))
    ra = rho_star*np.sin(np.radians(theta))
    return np.array([ra,dec])
    
def proper_motion(RA1,DEC1,t1,RA2,DEC2,t2):
    #Returns proper motion in mas/year
    #Inputs must be in mas and years!
    delta_RA = RA2 - RA1
    delta_DEC = DEC2 - DEC1
    delta_t = t2 - t1
    # ask mike which DEC to use in cos
    mu_RA = (delta_RA/delta_t) #* np.cos(((delta_DEC*u.mas).to(u.rad)))
    mu_DEC = delta_DEC/delta_t
    
    return [mu_RA, mu_DEC]

def calc_chi_sq(sim,data_row):
    """
    Asks the question: what would the chi-squared have been for this planet if 
    it had the same uncertainties as the Gaia and HIPPARCOS data?
    """
    
    # simulation values of RA and DEC between gaia dates
    sim_pm_gaia = proper_motion(sim[0,1], sim[1,1], 2015.0, sim[0,2], sim[1,2], 
                                2017.0)
    
    # simulation values of RA and DEC between hipparcos and gaia
    sim_pm_hg = proper_motion(sim[0,0], sim[1,0], 1991.25, sim[0,2], sim[1,2], 
                              2017.0)

    chi2 = 0
    chi2 += (sim_pm_gaia[0]-sim_pm_hg[0])**2/(data_row['pmra_hg_error']**2 + 
                                              data_row['pmra_gaia_error']**2)
    chi2 += (sim_pm_gaia[1]-sim_pm_hg[1])**2/(data_row['pmdec_hg_error']**2+ 
                                              data_row['pmdec_gaia_error']**2)
    
    return chi2

if __name__=="__main__":
    #For testing, you may not want all stars!
    test_binary_orbit=False
    data_dir = './'
    #data_dir = '/Users/mireland/Google Drive/EDR3_Planets/'
    n_sims = 1000
    chi2_threshold=10

    #Use -1 for all stars.
    n_stars = -1
    
    if test_binary_orbit:
        rho, theta, vr = bo.binary_orbit([2021,365,1,0,0,0,0],
                                         np.array([1991.25,2015,2017]))

        plt.figure()
        plt.clf()
        plt.axis('equal')
        plt.plot(rho*np.cos(np.radians(theta)), rho*np.sin(np.radians(theta)))

    # hip ids for binaries
    # binary_hip_id = np.array([1292, 1598,3588, 5110,12158,16846,17666,17749,
    #                          18512,21482,26501,29568,31246,32723,41211,54155,
    #                          56242,56838,60370,64583,77358,79607,86282,92919,
    #                          93926,94336,96895,97640,97944,109812,113421])
    # # hip ids for binaries, and likely binaries
    # binary_hip_id = np.array([1292, 1598,3588, 5110,12158,16846,17666,17749,
    #                           18512,21482,26501,29568,31246,32723,41211,54155,
    #                           56242,56838,60370,64583,77358,79607,86282,92919,
    #                           93926,94336,96895,97640,97944,109812,113421,
    #                           22431,43410,45836,51138,54922,54952,56280,
    #                           60074,95730,111686,113701,30314])
    
    # Import the data, all stars in 30pc
    data_30pc = Table.read(data_dir + 'data_30_pc.fits')
    bp_rp = np.array([data_30pc[i]['bp_rp'] for i in range(len(data_30pc))])
    
    #Cut down the list of stars based on a mass cut
    in_mass_range = (bp_rp >  0.33) & (bp_rp < 1.84)
    bp_rp = bp_rp[in_mass_range]
    data_30pc = data_30pc[in_mass_range]
    # not_binary=np.array([True if data_30pc[i]['hip_id'] not in binary_hip_id 
    #                      else False for i in range(len(data_30pc))])
    
    # data_30pc = data_30pc[not_binary]
    
    # choose data with desired chisq range, and also with reasonale
    data_30pc_chisq_10_100 = data_30pc[np.where((data_30pc['parallax_gaia']>33) & 
        (data_30pc['chisq']>chi2_threshold) & (data_30pc['chisq']<100) &
        (data_30pc['pmra_gaia_error']<0.1) & 
        (data_30pc['pmdec_gaia_error']<0.1))[0]]
    
    data_30pc_chisq_100_1000 = data_30pc[np.where((data_30pc['parallax_gaia']>33) & 
        (data_30pc['chisq']>100) & (data_30pc['chisq']<1000) &
        (data_30pc['pmra_gaia_error']<0.1) & 
        (data_30pc['pmdec_gaia_error']<0.1))[0]]
    # have this line if want to restrict the data set
   # bp_rp = np.array([data_30pc_chisq[i]['bp_rp']                             # !!
   #                   for i in range(len(data_30pc_chisq))])

    if n_stars==-1:
        n_stars = len(data_30pc)                                                #!!

    # trial plot simulation with first star
    all_star_radecs = []
    all_star_params = []
    j = 0
    for this_bprp, parallax in zip(bp_rp[:n_stars], data_30pc[:n_stars]  # !!
                                   ['parallax_gaia']):
        
        radecs = []
        params = []
        m_star = models.interpolate_mass(this_bprp,ddir=data_dir)
        t_k = np.array([1991.25,2015.0,2017.0])
        for i in range(n_sims): 
            m_planet = models.interpolate_planet_mass(num=1)[0]
            a = models.interpolate_semimajoraxis(num=1)[0]
            T = models.interpolate_period([a],[m_star])[0]
            omega = models.interpolate_omega(num=1)[0]
            Omega = models.interpolate_Omega(num=1)[0]
            i = models.interpolate_inclination(num=1)[0]
            e = models.interpolate_eccentricity(num=1)[0]
            T0 = 2000 + np.random.uniform()*T
            radec_sim = simulate_orbit(T0,T,a,e,Omega,omega,i,t_k,
                                         m_star,m_planet,parallax)
            radecs.append(radec_sim)

            chi_for_sim= calc_chi_sq(radec_sim,data_30pc[j])
            param = np.array([chi_for_sim,a,m_planet])
            
            params.append(param)
        all_star_radecs.append(radecs)
        all_star_params.append(params) 
        j +=1
    all_star_radecs = np.array(all_star_radecs)
    all_star_params = np.array(all_star_params, dtype = object)
    np.save('sim_params.npy',all_star_params)                       #!!
    f = open('chisq.txt','w')                                       #!!
    chis = np.empty((n_stars, n_sims))
    for i in range(n_stars):
        for j in range(n_sims):
            chis[i,j] = calc_chi_sq(all_star_radecs[i,j],data_30pc[i])#!!
            f.write(str(chis[i,j])+'\n')
    f.close()
    
planet_frequency = 0.1968

all_star_planet_detections = []

for i in range(n_stars):
    number_detections = np.sum(chis[i,:]>chi2_threshold)/n_sims * planet_frequency
    all_star_planet_detections.append(number_detections)
    