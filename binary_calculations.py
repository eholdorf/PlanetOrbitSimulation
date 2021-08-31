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
    ra_pm_diff = ra_pm_diff.to(u.rad)
    dec_pm_diff = dec_pm_diff.to(u.rad)
    return np.degrees(np.arctan2(ra_pm_diff, dec_pm_diff))

def acceleration(ra_pm_diff, dec_pm_diff, dt, distance):
    ra_pm_diff = ra_pm_diff.to(u.rad)/u.rad
    dec_pm_diff = dec_pm_diff.to(u.rad)/u.rad
    return (distance * np.sqrt(ra_pm_diff**2 + dec_pm_diff**2)/(dt**2)).to(
        u.nm/u.s**2)

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
        ra_pm_diff = (star['pmra_gaia_error'] - star['pmra_hg_error'])*u.mas
        dec_pm_diff = (star['pmdec_gaia_error'] - star['pmdec_hg_error'])*u.mas
        
        distance = (1000/star['parallax_gaia'])*u.pc
        
        pos_angle = position_angle(ra_pm_diff,dec_pm_diff)
        accel = acceleration(ra_pm_diff, dec_pm_diff,dt, distance)
        position_angles.append(pos_angle)
        accelerations.append(accel)
