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
    pa =  np.degrees(np.arctan2(ra_pm_diff, dec_pm_diff))
    return pa

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
    
    return (accel,accel*2*np.pi,accel_N)
    
if __name__=="__main__":
    data_dir = './'
    #data_dir = '/Users/mireland/Google Drive/EDR3_Planets/'
    n_sims = 1000
    chi2_threshold=10
    
    binary_hip_id = np.array([1292,1598,3588,5110,5260,12158,13081,13642,16846,
                         17666,17749,18512,20552,21482,22431,22715,24783,
                         26501,29568,31246,32723,38931,39826,41211,42173,
                         42697,43233,43410,43422,43557,44295,45343,45836,
                         47080,48133,51138,51248,51525,54155,54922,54952,
                         56242,56280,56452,56809,56838,59707,60074,60370,
                         60759,64532,64583,65343,66238,67308,69965,73184,
                         73633,73941,79492,79607,79702,80008,80925,83020,
                         86282,91605,92919,93926,94336,95730,96895,97640,
                         97944,98677,103768,104214,108162,109812,110109,
                         111686,113421,113701,117542,118310,77358])
 
    
    # Import the data, all stars in 30pc
    data_30pc = Table.read(data_dir + 'data_30_pc.fits')
    
    data_30pc_all = Table.read(data_dir + 'data_30_pc.fits')
    
    data_all = Table.read(data_dir + 'HGCA_vEDR3.fits')
    
    bp_rp = np.array([data_30pc[i]['bp_rp'] for i in range(len(data_30pc))])
    
    is_binary=np.array([True if data_30pc[i]['hip_id'] in binary_hip_id 
                    else False for i in range(len(data_30pc))])
     
    #Cut down the list of stars based on a mass cut
    in_mass_range = (bp_rp >  0.33) & (bp_rp < 1.84)
    bp_rp = bp_rp[in_mass_range]
    data_30pc = data_30pc[in_mass_range]
    # choose data with desired chisq range, and also with reasonale
    data_30pc_chisq = data_30pc[np.where((data_30pc['parallax_gaia']>33) & 
        (data_30pc['chisq']>chi2_threshold) & (data_30pc['chisq']<1000) &
        (data_30pc['pmra_gaia_error']<0.1) & 
        (data_30pc['pmdec_gaia_error']<0.1))[0]]
    # have this line if want to restrict the data set
    # bp_rp = np.array([data_30pc_chisq[i]['bp_rp']                             # !!
    #                   for i in range(len(data_30pc_chisq))])
    is_binary=np.array([True if data_30pc[i]['hip_id'] in binary_hip_id 
                    else False for i in range(len(data_30pc))])
    
    data_30pc_binary = data_30pc[is_binary]
    
    dt = (2016.0 - np.mean([1991.25,2016.0]))*u.year
    
    accelerations = []
    position_angles = []
    # (hip_id_host,bp_rp companion,separation angle)
    companion_bp_rp = [(1292,1.2718449,4),(1598,0.823,22),
                       (3588,1.6613917,10),(5110,1.7336645,6.6),
                       (5260,5.004,1.8),(12158,0.823,0),
                       (13081,3.1132164,23.4),(13642,1.737998,9),
                       (16846,1.1802888,6.7),(17666,1.1327901,8.2),
                       (17749,2.1908627,17),(18512,2.4244537,11),
                       (20552,0.8566065,5.3),(21482,0.823,7.7),
                       (22431,1.1047196,8.28),(22715,2.1229582,4.6),
                       (24783,2.7237701,4.9),(26501,1.3389177,4.8),
                       (29568,2.7866783,24.9),(31246,2.2648344,10.1),
                       (32723,2.7953777,6.8),(38931,2.4318743,4.5),
                       (39826,2.4997482,10.8),(41211,3.111949,35.5),
                       (42173,0.6800256,12.5),(42697,0.823,5.7),
                       (43233,3.16,2.6),(43410,0.823,8.3),
                       (43422,1.9531021,1.9),(43557,3.255,3.4),
                       (44295,1.6499186,5),(45343,1.8464499,21.1),
                       (45836,1.609067,5.6),(47080,2.7309647,6.7),
                       (48133,2.756097561,1.2),(51138,0.823,5.5),
                       (51248,1.6217604,4.5),(51525,2.8133726,4.8),
                       (54155,0.823,4.7),(54922,3.0115,2.9),
                       (54952,2.4039812,8),(56242,1.4227114,15.7),
                       (56280,0.69677544,9.8),(56452,0.93699646,17),
                       (56809,1.1792297,12.5),(56838,1.3197966,3.5),
                       (59707,2.0616837,4.3),(60074,0.823,6.3),
                       (60370,3.3300438,17.8),(60759,1.363409,2.8),
                       (64532,0.823,1.6),(64583,0.823,4.6),
                       (65343,1.7998562,1.7),(66238,12.66,3.3),
                       (67308,1.7068429,9.3),(69965,1.1960459,3.6),
                       (73184,0.823,26.2),(73633,2.5188332,4.5),
                       (73941,1.823801,4.4),(79492,0.9383755,4),
                       (79607,0.7817273,7.2),(79702,1.4445858,3.5),
                       (80008,2.45,2.1),(80925,1.0724525,4.1),
                       (83020,2.1911469,5.1),(86282,1.757102,4.2),
                       (91605,0.823,9.3),(92919,0.823,0),
                       (93926,0.823,0),(94336,0.8281636,9),
                       (95730,1.7640142,2.1),(96895,3.95,3.5),
                       (97640,0.823,0),(97944,0.823,0),
                       (98677,3.43,3.2),(103768,3.69,2.7),
                       (104214,1.7153406,31.8),(108162,2.428772,7.6),
                       (109812,1.58,3),(110109,0.823,2.5),
                       (111686,1.8069353,1.4),(113421,0.823,0),
                       (113701,2.57617,8.1),(117542,3.398,2.9),
                       (118310,0.823,4.9),(77358,0.104569435,15.2)]
    non_companion_bp_rp = [(13513,2.5174942,37),(22498,	2.3427963,44),
                           (26369,1.0027103,8.8),(34950,3.55,73),
                           (38594,1.1156178,400),(98767,3.1665573,179)]
    acceleration_comparison = []
    for companion in non_companion_bp_rp:
        accel_comparison_single = (
        compare_acceleration(companion[0],companion[1], companion[2]*u.arcsec))
        acceleration_comparison.append(accel_comparison_single)
        
    to_check = data_30pc_binary
    for star in to_check:
        ra_pm_diff=(star['pmra_gaia_error']-star['pmra_hg_error'])*u.mas/u.yr
        dec_pm_diff=(star['pmdec_gaia_error']-star['pmdec_hg_error'])*u.mas/u.yr
        
        distance = (1000/star['parallax_gaia'])*u.pc
        
        pos_angle = position_angle(ra_pm_diff,dec_pm_diff)
        accel = acceleration(ra_pm_diff, dec_pm_diff,dt, distance)
        position_angles.append((star['hip_id'],pos_angle))
        accelerations.append((star['hip_id'],accel))
