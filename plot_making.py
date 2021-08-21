# -*- coding: utf-8 -*-
"""
Created on Sat Aug 14 15:56:22 2021

@author: Erin
"""
import model_parameters as models
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from astropy.table import Table

# Import the data, all stars in 30pc
data_dir = './'
chi2_threshold = 9.5

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
#bp_rp = np.array([data_30pc_chisq[i]['bp_rp']                                   # !!
#                  for i in range(len(data_30pc_chisq))])
    
x = np.linspace(0,20,1000)
plt.figure(1)
plt.clf()
plt.plot(x, x+np.sin(x),'r', label = 'Proper & Companion Motion')
plt.ylabel('Right Ascension or Declination')
plt.xlabel('Time')
plt.plot(x,x,'royalblue', label = 'Proper Motion')
plt.plot([np.pi,np.pi],[0,20],'k--')
plt.plot([3*np.pi,3*np.pi],[0,20],'k--')
plt.annotate("", xy=(np.pi, 15), xytext=(3*np.pi, 15), 
             arrowprops=dict(arrowstyle="<->"))
plt.annotate("Companion Period", (np.pi+0.2,16))
plt.legend(loc = 'best')

plt.figure(2)
plt.clf()
plt.plot(x[0:500], np.sin(x[0:500]),'k', label = 'Companion Motion')
plt.plot([3],[np.sin(3)],'r.')
plt.annotate('Hipparcos',xy = (3,np.sin(3)+0.07))
plt.plot([6],[np.sin(6)],'r.')
plt.annotate('Gaia edr3',xy = (6,np.sin(6)-0.11))
plt.annotate('',xy=(3,np.sin(3)), xytext=(6,np.sin(6)),
             arrowprops = dict(arrowstyle='-'))
plt.ylabel('Right Ascension or Declination')
plt.xlabel('Time')

distances = [1000/data_30pc[i]['parallax_gaia']                           #!!               
             for i in range(len(data_30pc))]                                #!!
masses = models.interpolate_mass(bp_rp)

plt.figure(3)
plt.clf()
plt.plot(distances,masses,'k.')
plt.xlabel('Distance (pc)')
plt.ylabel(r'Mass ($M_\odot$)')

f = open('chisq.txt','r')
# change these values depending on what simulation run, assuming all stars
# and 1000 simulations run on each
n_stars = len(bp_rp)
n_sims = 1000
chis = np.zeros((n_stars, n_sims))
i = 0
j = 0
for line in f.readlines():
    chis[i,j] = float(line) 
    j += 1
    if j >= n_sims:
        i += 1
        j = 0
f.close()

plt.figure(4)
plt.clf()
plt.hist(chis[0],bins = 200)
plt.xlim(0,50)

nums = []
n_sims = 1000
for val in chis:
    a = 0
    for elem in val:
        if elem > chi2_threshold:
            a += 1
    nums.append(a/n_sims)  
nums = np.array(nums)

plt.figure(5)
plt.clf()
plt.plot(1000/data_30pc['parallax_gaia'],nums,'k.')                             #!!
plt.xlabel('distance')
plt.ylabel('proportion of chisq > 9.5')

plt.figure(6)
plt.clf()
plt.scatter(1000/data_30pc['parallax_gaia'], masses, c = 'k',                     #!!
            s=nums*30)
plt.xlabel('Distance (pc)')
plt.ylabel('Stellar Mass') 


planet_frequency = 0.1968

all_star_planet_detections = []

for i in range(n_stars):
    number_detections = np.sum(chis[i,:]<9.5)/n_sims * planet_frequency 
    all_star_planet_detections.append(number_detections)

all_star_planet_detections = np.array(all_star_planet_detections)

plt.figure(7)
plt.clf()
plt.plot(distances,all_star_planet_detections,'k.')
plt.xlabel('Distance (pc)')
plt.ylabel('Number of Detections')

exoplanet_distributions = Table.read('PS_2021.08.20_22.02.50.csv')

transit_distribution = exoplanet_distributions[np.where(
    exoplanet_distributions['discoverymethod']=='Transit')]

microlensing_distribution = exoplanet_distributions[np.where(
    exoplanet_distributions['discoverymethod']=='Microlensing')]

direct_imaging_distribution = exoplanet_distributions[np.where(
    exoplanet_distributions['discoverymethod']=='Imaging')]

radial_velocity_distribution = exoplanet_distributions[np.where(
    exoplanet_distributions['discoverymethod']=='Radial Velocity')]

plt.figure(8)
plt.clf()
plt.hist(transit_distribution['pl_orbsmax'],range=(0,1),bins = 75)
plt.xlabel('Semi-Major Axis (AU)')
plt.ylabel('Count')
#plt.savefig('transit_distribution.pdf')

plt.figure(9)
plt.clf()
plt.hist(microlensing_distribution['pl_orbsmax'],bins = 17)
plt.xlabel('Semi-Major Axis (AU)')
plt.ylabel('Count')
#plt.savefig('microlensing_distribution.pdf')

plt.figure(10)
plt.clf()
plt.hist(direct_imaging_distribution['pl_orbsmax'],range = (0,30),bins = 10)
plt.xlabel('Semi-Major Axis (AU)')
plt.ylabel('Count')
#plt.savefig('direct_imaging_distribution.pdf')

plt.figure(11)
plt.clf()
plt.hist(radial_velocity_distribution['pl_orbsmax'],range = (0,5),bins = 30)
plt.xlabel('Semi-Major Axis (AU)')
plt.ylabel('Count')
#plt.savefig('radial_velocity_distribution.pdf')

fig, axs = plt.subplots(2,2)
axs[0,0].hist(direct_imaging_distribution['pl_orbsmax'],range = (0,30),
              bins = 10)
axs[0,0].set_title('Direct Imaging')
axs[0,1].hist(radial_velocity_distribution['pl_orbsmax'],range = (0,5),
              bins = 10)
axs[0,1].set_title('Radial Velocity')
axs[1,0].hist(transit_distribution['pl_orbsmax'],range=(0,1),bins = 75)
axs[1,0].set_title('Transit')
axs[1,1].hist(microlensing_distribution['pl_orbsmax'],bins = 17)
axs[1,1].set_title('Microlensing')

for ax in axs.flat:
    ax.set(xlabel='Semi-Major Axis (AU)', ylabel='Count')
fig.tight_layout()

#fig.savefig('semi_major_axis_distn.pdf')

