# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 13:52:35 2021

@author: Erin
"""
import binary_orbit as bo
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as scipygamma

rho, theta, vr = bo.binary_orbit([0,365,1,0,0,0.3,45],np.array(np.linspace(0,1000)),50)
plt.plot(rho*np.cos(np.radians(theta)), rho*np.sin(np.radians(theta)))


def Fulton_model(a, C = 350, beta = -0.86,a0 = 3.6, gamma = 1.59):
    '''

    Parameters
    ----------
    a : number in units of au
        the size of the semi-major axis.
    C : number, optional
        normalisation constant. The default is 350.
    beta : number, optional
        occurence power law index after breaking point. The default is -0.86.
    a0 : number, optional
        semi-major axis location of the breaking point. The default is 3.6.
    gamma : number, optional
        power law index within breaking point - beta. The default is 1.59.

    Returns
    -------
    Number - occurence of planet

    '''
    occurence= C * a**beta * (1-np.exp(-(a/a0)**gamma))
    return occurence

SemiMajorAxis = np.linspace(0.01,30,200)
occurence = Fulton_model(SemiMajorAxis)

plt.figure()
plt.clf()
plt.plot(SemiMajorAxis,occurence)
plt.xlabel('Semi Major Axis (au)')
plt.ylabel('Occurence')


# =============================================================================
# def cumulative_Fulton_distribution(a, C = 350, beta = -0.86,a0 = 3.6, gamma = 1.59):
#     '''
# 
#     Parameters
#     ----------
#     a : list with entries as numbers in units of au
#         the size of the semi-major axis.
#     C : number, optional
#         normalisation constant. The default is 350.
#     beta : number, optional
#         occurence power law index after breaking point. The default is -0.86.
#     a0 : number, optional
#         semi-major axis location of the breaking point. The default is 3.6.
#     gamma : number, optional
#         power law index within breaking point - beta. The default is 1.59.
# 
#     Returns
#     -------
#     Number - occurence of planet
# 
#     '''
#     distn = []
#     for num in a:
#         val = Fulton_model(num,C,beta,a0,gamma)
#         distn.append(distn[+val)
#     return distn
# =============================================================================

CumulativeDistn = np.cumsum(occurence)
plt.figure()
plt.clf()
plt.plot(SemiMajorAxis,CumulativeDistn)
plt.xlabel('Semi Major Axis (au)')
plt.ylabel('Cummulative Occurence')
plt.title('Cumulative Distribution')

plt.figure()
plt.clf()
plt.plot(SemiMajorAxis,CumulativeDistn/max(CumulativeDistn))
plt.xlabel('Semi Major Axis (au)')
plt.ylabel('Normalised Cummulative Occurence')
plt.title('Normalised Cumulative Distribution')


def eccentricity_model(e,a = 0.867,b=3.03):
    top = scipygamma(a+b)
    bottom= scipygamma(a)*scipygamma(b)
    return top/bottom* e**(a-1) * (1-e)**(b-1)

eccentricities = np.linspace(0.01,1,100)
EccentricityDistn = eccentricity_model(eccentricities)
plt.figure()
plt.clf()
plt.plot(eccentricities,EccentricityDistn )
plt.xlabel('Eccentricities')
plt.ylabel('Occurence')

EccentricitesCumulative = np.cumsum(EccentricityDistn)

plt.figure()
plt.clf()
plt.plot(eccentricities, EccentricitesCumulative/max(EccentricitesCumulative))
plt.xlabel('Eccentricities')
plt.ylabel('Cumulative Probability of Eccentricity')

# inclination distn propto sin(i) if edge on 90^o

inclinations = np.linspace(0,np.pi,1000)
InclinationDistn = np.sin(inclinations)
plt.figure()
plt.clf()
plt.plot(inclinations,InclinationDistn )
plt.xlabel('Inclination')
plt.ylabel('Occurence')

InclinationCumulative = np.cumsum(InclinationDistn)

plt.figure()
plt.clf()
plt.plot(inclinations, InclinationCumulative/max(InclinationCumulative))
plt.xlabel('Inclinations')
plt.ylabel('Cumulative Probability of Eccentricity')

# python lines for making a cumulative distribution 

test = np.arange(10) # some pdf
cumdist = np.cumsum(test).astype(float) # cumulative distn
cumdist /= cumdist[-1] # normalise the values
vals = np.interp(np.random.uniform(size=10000),cumdist,test) # interpolate 10 000 random values onto this distribution
plt.figure()
plt.clf()
plt.hist(vals) # verify 
