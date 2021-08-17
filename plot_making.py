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

bp_rp = []
f = open('bp_rps.txt','r')
for line in f.readlines():
    bp_rp.append(float(line))   
f.close()

dd = Table.read('HGCA_vEDR3.fits')
dd.add_column(np.linspace(0,0,len(dd)),name = 'bp_rp')

ww = np.where((dd['parallax_gaia']>33) & (dd['chisq']>16) & (dd['chisq']<100))[0]
ddw = [dd[i] for i in ww]

for i in range(len(bp_rp)):
    ddw[i]['bp_rp'] = bp_rp[i]
    
x = np.linspace(0,20,1000)
plt.figure(1)
plt.clf()
plt.plot(x, x+np.sin(x),'r', label = 'Proper & Companion Motion')
plt.ylabel('Right Ascension or Declination')
plt.xlabel('Time')
plt.plot(x,x,'royalblue', label = 'Proper Motion')
plt.plot([np.pi,np.pi],[0,20],'k--')
plt.plot([3*np.pi,3*np.pi],[0,20],'k--')
plt.annotate("", xy=(np.pi, 15), xytext=(3*np.pi, 15), arrowprops=dict(arrowstyle="<->"))
plt.annotate("Companion Period", (np.pi+0.2,16))
plt.legend(loc = 'best')

plt.figure(2)
plt.clf()
plt.plot(x[0:500], np.sin(x[0:500]),'k', label = 'Companion Motion')
plt.plot([3],[np.sin(3)],'r.')
plt.annotate('Hipparcos',xy = (3,np.sin(3)+0.07))
plt.plot([6],[np.sin(6)],'r.')
plt.annotate('Gaia edr3',xy = (6,np.sin(6)-0.11))
plt.annotate('',xy=(3,np.sin(3)), xytext=(6,np.sin(6)),arrowprops = dict(arrowstyle='-'))
plt.ylabel('Right Ascension or Declination')
plt.xlabel('Time')

distances = [1000/ddw[i]['parallax_gaia'] for i in range(len(ddw))]
masses = models.interpolate_mass(bp_rp)

plt.figure(3)
plt.clf()
plt.plot(distances,masses,'k.')
plt.xlabel('Distance (pc)')
plt.ylabel(r'Mass ($M_\odot$)')

f = open('chis.txt','r')
chis = []
for line in f.readlines():
    chis.append(float(line))  
f.close()

plt.figure()
plt.clf()
plt.hist(chis,bins = 200)
#plt.xlim(0,75)