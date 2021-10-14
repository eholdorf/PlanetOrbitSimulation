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
from astropy.coordinates import SkyCoord
plt.ion()

def position_angle(ra_pm_diff, dec_pm_diff,sky_coord = False):
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
    if sky_coord == False:
        ra_pm_diff = ra_pm_diff.to(u.rad/u.yr)
        dec_pm_diff = dec_pm_diff.to(u.rad/u.yr)
        pa =  np.degrees(np.arctan2(ra_pm_diff,dec_pm_diff))
    else:
        pa = position_angle(ra_pm_diff,dec_pm_diff).to(u.deg)
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
    
    return (hip_id,accel,accel_N,accel.value-accel_N.value)

#def companion_pa(ra_diff,dec_diff,pa2,T2):
def companion_pa(ra_host,ra_companion,dec_host,dec_companion,sky_coord = False):
    dec_diff = dec_companion - dec_host
    ra_diff = (ra_companion - ra_host)*np.cos(np.radians(dec_host))
    if sky_coord == False:
        pa = np.degrees(np.arctan2(np.radians(ra_diff),np.radians(dec_diff)))
    else:
        host = SkyCoord(ra_host*u.deg, dec_host*u.deg,frame='icrs')
        companion = SkyCoord(ra_companion*u.deg, dec_companion*u.deg,frame='icrs')
        pa = host.position_angle(companion).degree

    return pa
    
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
                         111686,113421,113701,117542,118310,77358,97944,29208])
 
    
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
                       (42173,0.6800256,12.5),(42697,2.4031868,5.7),
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
                       (65343,1.7998562,1.7),(66238,3.5,3.3),
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
                       (118310,0.823,4.9),(77358,0.104569435,15.2),
                       (29208,4.16,45.20),(64583,1.8238201,4.6),
                       (97944,3.4028444,41.6)]
    
    non_companion_bp_rp = [(13513,2.5174942,37),(22498,	2.3427963,44),
                           (26369,1.0027103,8.8),(34950,3.55,73),
                           (38594,1.1156178,400),(98767,3.1665573,179)]
    
    companions_ra_dec = [(1292,4.062898388828774,-79.8503334744852),
                         (3588,11.454623039010272,-41.91009936192158),
                         (5110,16.373992107158337,15.387623150174809),
                         (5260,16.80653181333637,-32.42990302588197),
                         (13081,42.04194929266767,27.07333447155978),
                         (13642,43.91308317587121,26.87135846439827),
                         (16846,54.195031019527484,0.587145201977596),
                         (17666,56.76446366591412,41.42279214714805),
                         (17749,57.008656122272384,68.67859060565098),
                         (18512,59.36949755306553,-1.1571447425364605),
                         (20552,66.04746621445787,-57.072217259228566),
                         (22431,72.40652087854936,-53.8811039087136),
                         (22715,73.26963747632014,22.23384233042776),
                         (24783,79.69741566123787,-21.39366198032197),
                         (26501,84.54018837077523,-46.10782739017265),
                          (29568,93.43901141965836,-23.868215811159086),
                          (31246,98.30084010323947,5.464334446898048),
                          (32723,102.40124070017366,35.1374611688146),
                          (38931,119.4912788146194,-0.8135016353709931),
                          (39826,122.0553664165581,21.101076538256066),
                          (41211,126.13981809120203,-3.7429508429656004),
                          (42173,128.9618235438334,6.61962210473276),
                          (43233,132.11000732066893,6.465671940367521),
                          (43422,132.67602494116778,7.8641244595586155),
                          (43557,133.06866908011042,8.062715034354722),
                          (44295,135.32226274260103,15.262939477049024),
                          (45343,138.5913117966687,52.68342944184038),
                          (45836,140.1836971614705,51.26549428582397),
                          (47080,143.9127091741506,35.809996593021495),
                          (48133,147.19283909351617,-52.614387010197234),
                          (51248,157.01736219007225,48.782002086829245),
                          (51525,157.8491324699124,45.52339445319689),
                          (54922,168.69872764763568,-23.107179608230158),
                          (54952,168.78985102981798,73.47721553167881),
                          (56242,172.9336262061169,14.367437934390086),
                          (56280,173.06672103365204,-29.262681479002637),
                          (56452,173.62327851191714,-32.83031725241797),
                          (56809,174.68010645411127,45.107484782201254),
                          (56838,174.7863861531965,-27.696627775942417),
                          (59707,183.67530780555754,-24.77566230824648),
                          (60370,185.6848359514601,-39.17392714850606),
                          (60759,186.80768310744068,27.02434012080881),
                          (65343,200.8847186546332,29.238599001019995),
                          (66238,203.6374056692441,-38.90901257846076),
                          (67308,206.92897157434948,-32.431131168201624),
                          (69965,214.7511116706934,-25.814436322669373),
                          (73633,225.77622957721366,-41.99237273574981),
                          (73941,226.6453620961607,36.45733473004558),
                          (79492,243.32764025304397,13.526167955854836),
                          (79607,243.66672417598053,33.857177734562256),
                          (79702,243.98882744491237,7.355867188848203),
                          (80008,244.97996875599816,39.708051701644976),
                          (80925,247.87332100384344,-39.015645055130406),
                          (83020,254.47248277202303,47.36853477385246),
                          (86282,264.45198486411437,22.95376522107837),
                          (94336,288.01821557984215,49.85668878637743),
                          (95730,292.06409274077214,12.536099964498035),
                          (96895,295.45226847798733,50.523559441942325),
                          (98677,300.6424929932435,15.589658963597557),
                          (103768,315.4150952449052,-32.52512368646832),
                          (104214,316.753662752556,38.75607277205679),
                          (108162,328.71099285330774,-77.33745756697532),
                          (109812,333.63077733160594,27.857545343181602),
                          (111686,339.37125759857753,-67.1040210920978),
                          (113701,345.42363652870927,-51.46672871975605),
                          (117542,357.563573253572,-29.40099238299216),
                          (77358,236.87270043635263,-37.919946060994015),
                          (42697,130.52741860853445,-42.94233325218802),
                          (29208,92.39986752504565,5.668524106172984),
                          (64583,198.56189135258217,-59.10272664847473),
                          (97944,298.5853870969822,-23.946408599564357)]
    
    acceleration_comparison = []
    for companion in companion_bp_rp:
        accel_comparison_single = (
        compare_acceleration(companion[0],companion[1], companion[2]*u.arcsec))
        acceleration_comparison.append(accel_comparison_single)
        
    to_check = data_30pc_binary
    for star in to_check:
        ra_pm_diff=(star['pmra_gaia']-star['pmra_hg'])*u.mas/u.yr
        dec_pm_diff=(star['pmdec_gaia']-star['pmdec_hg'])*u.mas/u.yr
        
        distance = (1000/star['parallax_gaia'])*u.pc
        
        pos_angle = position_angle(ra_pm_diff,dec_pm_diff,sky_coord = True)
        accel = acceleration(ra_pm_diff, dec_pm_diff,dt, distance)
        position_angles.append((star['hip_id'],pos_angle))
        accelerations.append((star['hip_id'],accel))
        
    pa_comparison = []    
    for companion in companions_ra_dec:
        star = data_30pc[np.where(data_30pc['hip_id']==companion[0])]
        pa2 =[angle[1] for angle in position_angles if angle[0]==companion[0]]
        #pa1 = companion_pa(ra_diff,dec_diff,pa2[0].value,T2)
        #pa1 = companion_pa(ra_diff,dec_diff,companion[4],companion[3]) 
        pa1 = companion_pa(star['gaia_ra'],companion[1],star['gaia_dec'],
                           companion[2],sky_coord = True)
        print(companion[0],pa1,pa2)
        pa_comparison.append((companion[0],(pa1[0])%360,pa2[0].value%360,
                             ((pa2[0].value%360-(pa1[0])%360)+180)%360-180))
        
        