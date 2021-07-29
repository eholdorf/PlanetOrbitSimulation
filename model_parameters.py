# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 13:55:17 2021

@author: Erin
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as scipygamma
import astropy.units as u


def Fulton_model(a, C = 350, beta = -0.86,a0 = 3.6, gamma = 1.59):
    '''

    Parameters
    ----------
    a : type - number in units of au
        the size of the semi-major axis.
    C : type - number, optional
        normalisation constant. The default is 350.
    beta : type - number, optional
        occurence power law index after breaking point. The default is -0.86.
    a0 : type - number, optional
        semi-major axis location of the breaking point. The default is 3.6.
    gamma : type - number, optional
        power law index within breaking point - beta. The default is 1.59.

    Returns
    -------
    type - number 
        occurence of  with the given semi-major axis size.

    '''
    occurence= C * a**beta * (1-np.exp(-(a/a0)**gamma))
    return occurence

# set of semi major axis lengths to look at 
SemiMajorAxis = np.linspace(0.01,30,200)

# occurence of planets at different semi-major axis lengths according to 
# Fulton model
FultonOccurence = Fulton_model(SemiMajorAxis)

# plot the Fulton Model
plt.figure()
plt.clf()
plt.plot(SemiMajorAxis,FultonOccurence)
plt.xlabel('Semi Major Axis (au)')
plt.ylabel('Occurence')

# cumulative distribution of planet occurence with Fulton Model
OccurenceCumulativeDistn = np.cumsum(FultonOccurence)

# plot cumulative distribution of planet occurence with Fulton Model
plt.figure()
plt.clf()
plt.plot(SemiMajorAxis,OccurenceCumulativeDistn)
plt.xlabel('Semi Major Axis (au)')
plt.ylabel('Cummulative Occurence')
plt.title('Cumulative Distribution')

# normalise the cumulative distribution of planet occurence with Fulton Model
OccurenceCumulativeDistn /= max(OccurenceCumulativeDistn)

# plot normalised cumulative distribution of planet occurence with Fulton Model
plt.figure()
plt.clf()
plt.plot(SemiMajorAxis, OccurenceCumulativeDistn)
plt.xlabel('Semi Major Axis (au)')
plt.ylabel('Normalised Cummulative Occurence')
plt.title('Normalised Cumulative Distribution')


def eccentricity_model(e,a = 0.867,b=3.03):
    '''

    Parameters
    ----------
    e : type - number
        the eccentricity, a value between 0 and 1
    a : type - number, optional
        parameter one. The default is 0.867.
    b : type - number, optional
        parameter two. The default is 3.03.

    Returns
    -------
    type - number
        The occurence of this eccentricity.
        
    Model From - https://academic.oup.com/mnrasl/article/434/1/L51/1166128
    '''
    top = scipygamma(a+b)
    bottom= scipygamma(a)*scipygamma(b)
    return top/bottom * e**(a-1) * (1-e)**(b-1)

# set of eccentricities
Eccentricities = np.linspace(0.001,1,1000)

# occurence of each eccentricity 
EccentricityDistn = eccentricity_model(Eccentricities)

# plotting the eccentrictiy model
plt.figure()
plt.clf()
plt.plot(Eccentricities,EccentricityDistn )
plt.xlabel('Eccentricities')
plt.ylabel('Occurence')

# find the cumulative distribution of eccentricites 
EccentricityCumulative = np.cumsum(EccentricityDistn)

# normalise the cumulative distribution of eccentricities
EccentricityCumulative /= EccentricityCumulative[-1]

# plot the normalised distribution of eccentricities 
plt.figure()
plt.clf()
plt.plot(Eccentricities, EccentricityCumulative)
plt.xlabel('Eccentricities')
plt.ylabel('Cumulative Probability of Eccentricity')

# mike's way of interpolating values
EccentricityInterp = np.interp(np.random.uniform(size = 10000),
                              EccentricityCumulative,EccentricityDistn)

# my way interpolating values
EccentricityInterp2 = np.random.choice(np.linspace(0,1,1000),size = 10000, 
                         p = EccentricityDistn/sum(EccentricityDistn))

def interpolate_eccentricity(low = 0.001,high = 1,num = 10000):
    '''

    Parameters
    ----------
    low : TYPE, optional
        DESCRIPTION. The default is 0.001.
    high : TYPE, optional
        DESCRIPTION. The default is 1.
    num : TYPE, optional
        DESCRIPTION. The default is 10000.

    Returns
    -------
    EccentricityInterp2 : TYPE
        DESCRIPTION.

    '''
    
    # eccentricities to interpolate with random values
    e = np.linspace(low,high,num)
    
    # occurence of each eccentricity 
    EccentricityDistn = eccentricity_model(e)
    
    # find the cumulative distribution of eccentricites 
    EccentricityCumulative = np.cumsum(EccentricityDistn)
    
    # normalise the cumulative distribution of eccentricities
    EccentricityCumulative /= EccentricityCumulative[-1]
    
    
    # interpolating values
    EccentricityInterp2 = np.random.choice(e,size = num, 
                             p = EccentricityDistn/sum(EccentricityDistn))
    
    return EccentricityInterp2
    
def inclination_model(i, edge_angle = 90):
    '''
    Parameters
    ----------
    i : type - number
        inclination angle in degrees.
        
    edge_angle : type - number, optional
        the angle to consider as edge on in degrees, choose 0 or 90. 
        The defult is 90.

    Returns
    -------
    type - number
        occurence of inclination angle.
        
    Model from https://arxiv.org/pdf/1805.08211.pdf and mike

    '''
    if edge_angle == 90:
        i = np.deg2rad(i)
        return np.sin(i)
    elif edge_angle == 0:
        i = np.deg2rad(i)
        return np.cos(i)
    else:
        return 'Not valid edge_angle, should be 0 deg or 90 deg.'

# set of inclination angles

inclinations = np.linspace(0,180,1000)

# inclination distribution
InclinationDistn = inclination_model(inclinations)

# plot of inclination angles
plt.figure()
plt.clf()
plt.plot(inclinations,InclinationDistn)
plt.xlabel('Inclination')
plt.ylabel('Occurence')

# cumulative distribution of inclination angles
InclinationCumulative = np.cumsum(InclinationDistn)

# normalising cumulative distribution of inclination angles
InclinationCumulative /= InclinationCumulative[-1]

# plotting the normalised cumulative distribution of inclination angles
plt.figure()
plt.clf()
plt.plot(inclinations, InclinationCumulative)
plt.xlabel('Inclinations')
plt.ylabel('Cumulative Probability of Inclination')

# mike's way of interpolating values
InclinationInterp = np.interp(np.random.uniform(size = 10000),
                              InclinationCumulative,InclinationDistn)

# my way interpolating values
InclinationInterp2 = np.random.choice(np.linspace(0,180,1000),size = 10000, 
                         p = InclinationDistn/sum(InclinationDistn))

def interpolate_inclination(low = 0,high = 180,num = 10000):
    '''

    Parameters
    ----------
    low : TYPE, optional
        DESCRIPTION. The default is 0.001.
    high : TYPE, optional
        DESCRIPTION. The default is 1.
    num : TYPE, optional
        DESCRIPTION. The default is 10000.

    Returns
    -------
    EccentricityInterp2 : TYPE
        DESCRIPTION.

    '''
    
    # eccentricities to interpolate with random values
    i = np.linspace(low,high,num)
    
    # inclination distribution
    InclinationDistn = inclination_model(i)
    
    # cumulative distribution of inclination angles
    InclinationCumulative = np.cumsum(InclinationDistn)
    
    # normalising cumulative distribution of inclination angles
    InclinationCumulative /= InclinationCumulative[-1]
    
    # my way interpolating values
    InclinationInterp2 = np.random.choice(i, size = num, 
                             p = InclinationDistn/sum(InclinationDistn))
    
    return InclinationInterp2

def interpolate_Omega(low = 0, high = 360, num = 10000):
    '''

    Parameters
    ----------
    low : TYPE, optional
        DESCRIPTION. The default is 0.
    high : TYPE, optional
        DESCRIPTION. The default is 360.
    num : TYPE, optional
        DESCRIPTION. The default is 10000.

    Returns
    -------
    Omega : TYPE
        DESCRIPTION.

    '''
    Omega = np.linspace(low,high,num)
    return Omega

def interpolate_omega(low = 0, high = 180, num = 10000):
    '''

    Parameters
    ----------
    low : TYPE, optional
        DESCRIPTION. The default is 0.
    high : TYPE, optional
        DESCRIPTION. The default is 180.
    num : TYPE, optional
        DESCRIPTION. The default is 10000.

    Returns
    -------
    Omega : TYPE
        DESCRIPTION.

    '''
    Omega = np.linspace(low,high,num)
    return Omega

def interpolate_semimajoraxis(low = 3, high = 30, num = 10000):
    '''

    Parameters
    ----------
    low : TYPE, optional
        DESCRIPTION. The default is 3.
    high : TYPE, optional
        DESCRIPTION. The default is 30.
    num : TYPE, optional
        DESCRIPTION. The default is 10000.

    Returns
    -------
    a : TYPE
        DESCRIPTION.

    '''
    
    a = np.linspace(low,high,num)
    
    return a

def interpolate_period(a, num = 10000):
    '''

    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    num : TYPE, optional
        DESCRIPTION. The default is 10000.

    Returns
    -------
    T : TYPE
        DESCRIPTION.

    '''
    T = []
    
    for val in a:
        T.append((np.sqrt(val**3)*u.yr).to(u.day).value)
    
    return T

# -----------------------------------------------------------------------------
# testing making a cumulative distribution - mike's way
# -----------------------------------------------------------------------------

# some pdf
test = np.arange(10)

# cumulative distribution
cumdist = np.cumsum(test).astype(float)

# normalise the values
cumdist /= cumdist[-1] 

# interpolate 10 000 random values onto this distribution
vals = np.interp(np.random.uniform(size=10000),cumdist,test) 

# verify by plotting (should look like original distribution)
plt.figure()
plt.clf()
plt.hist(vals) 

# -----------------------------------------------------------------------------
# testing making a cumulative distribution - my way
# -----------------------------------------------------------------------------

# some pdf
test = np.arange(100)

# interpolate 10 000 random values onto this distribution
vals2 = np.random.choice(100,size=10000, p = test/sum(test))