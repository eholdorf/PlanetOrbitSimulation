# -*- coding: utf-8 -*-
"""
Created on Sat Aug 14 15:56:22 2021

@author: Erin
"""
import model_parameters as models
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u
import astropy.constants as c

# Import the data, all stars in 30pc
data_dir = './'
chi2_threshold = 10
all_params = np.load('sim_params.npy',allow_pickle = True)

# hip ids for binaries (wide, close, wide (look at acceleration))
binary_hip_id = np.array([1292,3588,5110,5260,13642,16846,17666,17749,18512,
                          20552,22431,22715,24783,26501,29568,31246,32723,
                          38931,39826,43233,43422,43557,44295,45343,45836,
                          47080,48133,51248,51525,54922,54952,56452,56809,
                          56838,59707,60759,64583,65343,66238,67308,69965,
                          73633,73941,77358,79492,79607,79702,80008,80925,
                          83020,86282,94336,95730,96895,98677,103768,104214,
                          108162,109812,111686,117542,
                          12158,13081,13513,21482,24205,26369,42173,55052,
                          60816,92919,97640,102119,118162,
                          29208,41211])
 
#Import the data, all stars in 30pc
data_30pc = Table.read(data_dir + 'data_30_pc.fits')
# not_binary=np.array([True if data_30pc[i]['hip_id'] not in binary_hip_id 
#                       else False for i in range(len(data_30pc))])
 
# data_30pc = data_30pc[not_binary]


# have this line if want to restrict the data set
# bp_rp = np.array([data_30pc_chisq[i]['bp_rp']                             # !!
#                   for i in range(len(data_30pc_chisq))])
bp_rp = np.array([data_30pc[i]['bp_rp'] for i in range(len(data_30pc))])
 
#Cut down the list of stars based on a mass cut
in_mass_range = (bp_rp >  0.33) & (bp_rp < 1.84)
bp_rp = bp_rp[in_mass_range]
data_30pc = data_30pc[in_mass_range]

# choose data with desired chisq range, and also with reasonale
data_30pc_chisq_10_100 = data_30pc[np.where((data_30pc['parallax_gaia']>33) & 
    (data_30pc['chisq']>chi2_threshold) & (data_30pc['chisq']<100) &
    (data_30pc['pmra_gaia_error']<0.1) & 
    (data_30pc['pmdec_gaia_error']<0.1))[0]]

data_30pc_chisq_100_1000 = data_30pc[np.where((data_30pc['parallax_gaia']>33) & 
    (data_30pc['chisq']>100) & (data_30pc['chisq']<1000) &
    (data_30pc['pmra_gaia_error']<0.1) & 
    (data_30pc['pmdec_gaia_error']<0.1))[0]]

data_30pc_chisq = data_30pc[np.where((data_30pc['parallax_gaia']>33) & 
    (data_30pc['chisq']>chi2_threshold) & (data_30pc['chisq']<1000) &
    (data_30pc['pmra_gaia_error']<0.1) & 
    (data_30pc['pmdec_gaia_error']<0.1))[0]]

x = np.linspace(0,20,1000)
plt.figure(1)
plt.clf()
plt.plot(x, x+np.sin(x),'r', label = 'Proper & Companion Motion')
plt.ylabel('Right Ascension or Declination (mas)')
plt.xlabel('Time (yrs)')
plt.plot(x,x,'royalblue', label = 'Proper Motion')
plt.plot([np.pi,np.pi],[0,20],'k--')
plt.plot([3*np.pi,3*np.pi],[0,20],'k--')
plt.annotate("", xy=(np.pi, 15), xytext=(3*np.pi, 15), 
             arrowprops=dict(arrowstyle="<->"))
plt.annotate("Companion Period", (np.pi+0.2,16))
plt.legend(loc = 'best')
#plt.savefig('proper_motion_mine.pdf')

plt.figure(2)
plt.clf()
plt.plot(x[0:500], np.sin(x[0:500]),'k', label = 'Companion Motion')
plt.plot([3],[np.sin(3)],'r.')
plt.annotate('Hipparcos',xy = (3,np.sin(3)+0.07))
plt.plot([6],[np.sin(6)],'r.')
plt.annotate('Gaia edr3',xy = (6,np.sin(6)-0.11))
plt.annotate('',xy=(3,np.sin(3)), xytext=(6,np.sin(6)),
             arrowprops = dict(arrowstyle='-'))
plt.ylabel('Right Ascension or Declination (mas)')
plt.xlabel('Time (yrs)')
#plt.savefig('planet_motion_mine.pdf')

distances = [1000/data_30pc[i]['parallax_gaia']                           #!!               
             for i in range(len(data_30pc))]                                #!!
masses = models.interpolate_mass(bp_rp)

plt.figure(3)
plt.clf()
plt.plot(distances,masses,'k.')
plt.xlabel('Distance (pc)')
plt.ylabel(r'Mass ($M_\odot$)')
#plt.savefig('mass_distance_stars_used.pdf')
                                                      #!!
# change these values depending on what simulation run, assuming all stars
# and 1000 simulations run on each
n_stars = len(data_30pc)
n_sims = 1000
chis = np.zeros((n_stars, n_sims))
i = 0
for star_param in all_params:
    j = 0
    for sim_param in star_param:
        chis[i,j] = sim_param[0]
        j+= 1
    i += 1

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
plt.ylabel(r'Stellar Mass ($M_\odot$)') 
#plt.savefig('mass_distance_stars_used.pdf')


planet_frequency = 0.2233 

all_star_planet_detections = []

for i in range(n_stars):
    number_detections = np.sum(chis[i,:]>chi2_threshold)/n_sims * planet_frequency 
    all_star_planet_detections.append(number_detections)

all_star_planet_detections = np.array(all_star_planet_detections)

plt.figure(7)
plt.clf()
plt.plot(distances,all_star_planet_detections/planet_frequency,'k.')
plt.xlabel('Distance (pc)')
plt.ylabel('Likihood of Planet Detection')
#plt.savefig('number_detections.pdf')

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

# #fig.savefig('semi_major_axis_distn.pdf')

not_binary=np.array([True if data_30pc_chisq[i]['hip_id'] not in binary_hip_id 
                        else False for i in range(len(data_30pc_chisq))])
 
#data_30pc_chisq = data_30pc_chisq[not_binary]
chi_sq_gaia = data_30pc_chisq['chisq']
chi_sq_sim = []
for star in chis:
    for sim in star:
        if sim > chi2_threshold:
            chi_sq_sim.append(sim)

chi_sq_gaia = np.log10(chi_sq_gaia)
chi_sq_sim = np.log10(chi_sq_sim)
          
plt.figure(12)
plt.clf()
plt.hist(chi_sq_gaia, alpha = 0.5,density = True)
plt.hist(chi_sq_sim, alpha = 0.5,density = True)
plt.legend(['Gaia','Simulation'])

plt.figure(13)
plt.clf()
plt.hist([chi_sq_gaia,chi_sq_sim],density = True)
plt.plot([2,2],[0,1],'k')
plt.legend(['$\chi^2$ = 100','Gaia','Simulation'])
plt.xlabel('log$_{10}$($\chi^2$)')
plt.ylabel('Density')
plt.title('$\chi^2$ Distribution With Binaries')
plt.xlim(1,4)
#plt.savefig('chi_sq_distn_with_binary.pdf')

params = all_params
a = []
m_p = []

# only keep the stars that have a min chi sq that is greater than 9.5
above_threshold = [True if max(chis[i])>chi2_threshold else False for i in 
                   range(len(chis))]
# limit data
detections = data_30pc[above_threshold]
# limit chi sq
detection_chis = chis[above_threshold]
# find the index's of the maximum chi sq
max_chi_index = [list(detection_chis[i]).index(max(detection_chis[i])) 
                 for i in range(len(detection_chis))]
#find the max chi for each star
max_chi = [max(detection_chis[i]) for i in range(len(detection_chis))]

d = [1000/detections[i]['parallax_gaia'] for i in range(len(detections))]
params = params[above_threshold]

for i in range(len(params)):
    a.append(params[i][max_chi_index[i]][1])
    m_p.append(params[i][max_chi_index[i]][2])
        
m_p = np.array(m_p)
plt.figure(14)
plt.clf()
plt.scatter(d,a,s=np.exp(m_p)/5000,c=detections['bp_rp'],cmap = 'RdYlBu')
plt.xlabel('Distance (pc)')
plt.ylabel('Semi-Major Axis (AU)')
plt.colorbar(label =r'$B_p - R_p$')

#plt.savefig('distance_vs_a.pdf')

# plt.figure(15)
# plt.clf()
# i = 0
# used = []
# points = []
# for point in nums:
#     plt.scatter(1000/data_30pc['parallax_gaia'][i], masses[i], c = 'k',                     #!!
#             s=point*30+1,)
#     if point*100 in range(11) and 10 not in used:
#         points.append((i,10))
#         used.append(10)
#     elif point*100 in range(11,21) and 20 not in used:
#         points.append((i,20))
#         used.append(20)
#     elif point*100 in range(21,31) and 30 not in used:
#         points.append((i,30))
#         used.append(30)
#     elif point*100 in range(31,41) and 40 not in used:
#         points.append((i,40))
#         used.append(40)
#     elif point*100 in range(41,51) and 50 not in used:
#         points.append((i,50))
#         used.append(50) 
#     elif point*100 in range(51,61) and 60 not in used:
#         points.append((i,60))
#         used.append(60) 
#     elif point*100 in range(61,71) and 70 not in used:
#         points.append((i,70))
#         used.append(70) 
#     elif point*100 in range(71,81) and 80 not in used:
#         points.append((i,80))
#         used.append(80) 
#     elif point*100 in range(81,91) and 90 not in used:
#         points.append((i,90))
#         used.append(90) 
#     elif point*100 in range(91,101) and 100 not in used:
#         points.append((i,100))
#         used.append(100)
#     i += 1
# for elem in [10,20,40,80,90]:
#     i = [point[0] for point in points if point[1] == elem]
#     plt.scatter(1000/data_30pc['parallax_gaia'][i], masses[i], c = 'k',                     #!!
#                 s=nums[i]*30+1,label = str(nums[i]*100)+'%')
# plt.legend(loc='upper left',ncol = 3)
# plt.xlabel('Distance (pc)')
# plt.ylabel(r'Stellar Mass ($M_\odot$)') 
# plt.savefig('mass_distance_stars_used.pdf')
params_for_plot = np.load('sim_params_for_plot.npy',allow_pickle = True)
n_stars_plot = len(data_30pc)
n_sims_plot = 1
chis_plot = np.zeros((n_stars_plot, n_sims_plot))
i = 0
for star_param in params_for_plot:
    j = 0
    for sim_param in star_param:
        chis_plot[i,j] = sim_param[0]
        j+= 1
    i += 1

a = []
m_p = []

# only keep the stars that have a min chi sq that is greater than 10
#above_threshold = [True if max(chis[i])>chi2_threshold else False for i in 
#                   range(len(chis))]
best_planet = np.argmax(chis_plot, axis=1)

#Lets only consider the first planet here.
above_threshold = np.max(chis_plot, axis=1)>chi2_threshold
above_threshold = chis_plot[:,0]>chi2_threshold

# limit data
best_planet = best_planet[above_threshold]
detections = data_30pc[above_threshold]
# limit chi sq
detection_chis = chis_plot[above_threshold]

d = [1000/detections[i]['parallax_gaia'] for i in range(len(detections))]
params_for_plot = params_for_plot[above_threshold]

#Here we plot all stars *if* the planet is above the threshold
for i in range(len(params_for_plot)):
    a.append(params_for_plot[i][0][1])
    m_p.append(params_for_plot[i][0][2])

#!!! This is the mysterious figure !!!        
m_p = np.array(m_p)
plt.figure(16)
plt.clf()
plt.scatter(d,a,s=(m_p*10000),c=detections['bp_rp'],cmap = 'RdYlBu_r')
plt.xlabel('Distance (pc)')
plt.ylabel('Semi-Major Axis (AU)')
plt.colorbar(label =r'$B_p - R_p$')
#plt.savefig('distance_vs_a.pdf')
#!!! MJI Lets have a look at some of those mysterious systems.
#data_30pc is all planets
#above_threshold is True if at least 1 planet is above the threshold for that
#star. 
#Do these match? params_for_plot is only a and m_p from the simulation.
#plt.semilogy(params_for_plot[:,:,0].flatten(), chis_plot.flatten(),'.')
#No change between 3 and 12 AU...

#Here is a messay way to get a plot.
fig = plt.figure(19)
plt.clf()
ax = plt.gca()
im = ax.scatter(d,m_p*c.M_sun.to(u.M_jup).value, c=detections['bp_rp'],
                alpha = 0.8,edgecolors='none', cmap='RdYlBu_r')
dplot = np.arange(10,31)
ax.plot(dplot, dplot/20)
ax.axis([3,31,0.5,20])
ax.set_yscale('log')
plt.xlabel('Distance (pc)')
plt.ylabel(r'Mass (M$_{Jup}$)')
fig.colorbar(im,label =r'$B_p - R_p$',ax=ax)
plt.show()

plt.figure(17)
plt.clf()
plt.scatter(1000/data_30pc['parallax_gaia'],data_30pc['pmdec_gaia_error'],
            c = 'red',s=20,label=r'$\mu_{dec}$ Error')
plt.scatter(1000/data_30pc['parallax_gaia'],data_30pc['pmra_gaia_error'],
            c= 'black',s=20,alpha = 0.5,label=r'$\mu_{RA}$ Error')
plt.legend(loc='best')
plt.xlabel('Distance (pc)')
#plt.savefig('proper_motion_error.pdf')
plt.ylabel(r'Error (mas yr$^{-1}$)')

plt.figure(18)
plt.clf()
binary=np.array([True if data_30pc_chisq[i]['hip_id'] in binary_hip_id 
                        else False for i in range(len(data_30pc_chisq))])
 
data_30pc_chisq_binary = data_30pc_chisq[binary]
plt.hist(np.log10(data_30pc_chisq_binary['chisq']))
plt.plot([2,2],[0,13],'k',label = r'$\chi^2$ = 100')
plt.xlabel('log$_{10}$($\chi^2$)')
plt.ylabel('Number of Binaries')
plt.legend()
#plt.savefig('removed_binary_chisq.pdf')

C = [130,200,250,300,350,375,400]
C_detections = [70.26,72.39,74.02,75.88,77.82,77.91,78.47]
C_erry = [1/np.sqrt(C_detections[i]) for i in range(len(C))]
a_0 = [0.5,1.0,1.5,2.0,2.5,3.0,3.6]
a_0_detections = [75.80,76.06,76.94,77.74,78.26,78.15,77.82]
a_0_erry = [1/np.sqrt(a_0_detections[i]) for i in range(len(a_0))]

plt.figure(20)
plt.errorbar(C,C_detections,yerr = C_erry,fmt='k.',capsize = 3.0)
plt.xlabel('C')
plt.ylabel('Number Detections')
#plt.savefig('c_versus_detections.pdf')

plt.figure(21)
plt.errorbar(a_0,a_0_detections,yerr = a_0_erry,fmt='k.',capsize = 3.0)
plt.xlabel('a0')
plt.ylabel('Number Detections')
#plt.savefig('a0_versus_detections.pdf')