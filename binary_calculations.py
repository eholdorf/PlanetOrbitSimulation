# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 19:06:31 2021

@author: Erin
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u
import astropy.constants as c
import model_parameters
plt.ion()

def position_angle(ra_pm_diff, dec_pm_diff):
    '''

    Parameters
    ----------
    ra_pm_diff : units - angle
        The difference in right ascension.
    dec_pm_diff : units - angle
        The difference in declination.

    Returns
    -------
    units - degrees
        The position angle of the system

    '''
    ra_pm_diff = ra_pm_diff.to(u.rad/u.yr)
    dec_pm_diff = dec_pm_diff.to(u.rad/u.yr)
    return np.degrees(np.arctan2(ra_pm_diff, dec_pm_diff))

def acceleration(ra_pm_diff, dec_pm_diff, dt, distance):
    ra_pm_diff = ra_pm_diff.to(u.rad/u.s)/u.rad
    dec_pm_diff = dec_pm_diff.to(u.rad/u.s)/u.rad
    return (distance * np.sqrt(ra_pm_diff**2 + dec_pm_diff**2)/(dt)).to(
        u.nm/u.s**2)

def newton_acceleration(M,separation,d):
    return (c.G*M/(separation*d)**2).to(u.nm/(u.rad**2 * u.s**2))

def compare_acceleration(hip_id,bp_rp,sep):
    star = data_30pc[np.where(data_30pc['hip_id']==hip_id)]
    ra_pm_diff = (star['pmra_gaia_error'] - star['pmra_hg_error'])*u.mas/u.yr
    dec_pm_diff = (star['pmdec_gaia_error'] - star['pmdec_hg_error'])*u.mas/u.yr
    dt = (2016.0 - np.mean([1991.25,2016.0]))*u.year
    distance = (1000/star['parallax_gaia'])*u.pc
    accel = acceleration(ra_pm_diff, dec_pm_diff,dt, distance)
    M = model_parameters.interpolate_mass([bp_rp])[0]*u.M_sun
    accel_N = newton_acceleration(M,sep,distance)
    
    print(accel,accel_N)
    
if __name__=="__main__":
    data_dir = './'
    #data_dir = '/Users/mireland/Google Drive/EDR3_Planets/'
    n_sims = 1000
    chi2_threshold=9.5
    
    # Import the data, all stars in 30pc
    data_30pc = Table.read(data_dir + 'data_30_pc.fits')
    bp_rp = np.array([data_30pc[i]['bp_rp'] for i in range(len(data_30pc))])
     
    #Cut down the list of stars based on a mass cut
    in_mass_range = (bp_rp >  0.33) & (bp_rp < 1.84)
    bp_rp = bp_rp[in_mass_range]
    data_30pc = data_30pc[in_mass_range]
    # choose data with desired chisq range, and also with reasonale
    data_30pc_chisq = data_30pc[np.where((data_30pc['parallax_gaia']>33) & 
        (data_30pc['chisq']>chi2_threshold) & (data_30pc['chisq']<100) &
        (data_30pc['pmra_gaia_error']<0.1) & 
        (data_30pc['pmdec_gaia_error']<0.1))[0]]
    # have this line if want to restrict the data set
    # bp_rp = np.array([data_30pc_chisq[i]['bp_rp']                             # !!
    #                   for i in range(len(data_30pc_chisq))])
    
    dt = (2016.0 - np.mean([1991.25,2016.0]))*u.year
    
    accelerations = []
    position_angles = []
    
    for star in data_30pc_chisq:
        ra_pm_diff = (star['pmra_gaia_error'] - star['pmra_hg_error'])*u.mas/u.yr
        dec_pm_diff = (star['pmdec_gaia_error'] - star['pmdec_hg_error'])*u.mas/u.yr
        
        distance = (1000/star['parallax_gaia'])*u.pc
        
        pos_angle = position_angle(ra_pm_diff,dec_pm_diff)
        accel = acceleration(ra_pm_diff, dec_pm_diff,dt, distance)
        position_angles.append(pos_angle)
        accelerations.append(accel)
