#This function is based on Christian Hummel's 2001 Michelson School
#web document. It's designed to be used directly for fitting...
#
#It comes from "binary_position.pro", and IDL program.
#
#NB for differentiation, see mathematica/binary_diff.nb
#Done some timing tests on my laptop: this procedure is dominated
#by calculation, not interpretation when there are more than 30 jds.

#Meaning of negative numbers:
#Negative eccentricity: change \omege by 180 degrees and T by half a
# period.
#Negative inclination: the same as positive inclination!

import numpy as np
import matplotlib.pyplot as plt

def binary_orbit(params, jds, niter_anomaly=5, do_deriv=False):
    """Compute the separation and position angle of a binary given
    a list of epochs.
    
    Parameters
    ----------
    params: numpy array(7)
        T0: The time of periastron passage
        P:  The period
        a:  The semi-major axis
        e:  The eccentricity
        n:  Capital omega ( an orbital angle )
        w:  Little omega
        i:  The inclination
    jds: numpy array
        The list of dates to compute the binary position. Same time units as
        T0 and P
    niter_anomaly: int
        Number of fixed point iterations to compute the eccentric anomaly
    do_deriv: bool
        Do we compute the derivative of the orbital parameters? Needed for 
        explicit chi-squared Hessian computation (which can be super-useful)
        
    Returns
    -------
    (rho,theta,vr,[deriv])
        Separation in units of a, position angle in degrees, 
        velocity in units of semi major per unit per time  
    """
    #jds is a numpy array.
    t = jds-params[0]
    P = params[1]
    a = params[2]
    e = abs(params[3])
    n = params[4]*np.pi/180.
    w = params[5]*np.pi/180.
    i = params[6]*np.pi/180.
    #Allow a negative eccentricity to have the obvious meaning.
    if (params[3] < 0):
        t += P/2.
        w += np.pi

    #The mean anomaly 
    #(p,t) -> M 
    #Tr1 (transformation 1)
    #Sequence is (p,t,e) -> M -> bE -> nu -> alpha. We want: 
    #dAdp = dnudbE.dbEdM.dMdp
    #dAdt = dnudbE.dbEdM.dMdp
    #dAde = dnudbE.dbEde + dnude
    dMdt = -2*np.pi/P     #- sign because time _since_ periastron passage.
    dMdp = -2*np.pi*t/P**2 #NB Can be very large.
    M = 2*np.pi*(t % P)/P

    #The eccentric anomaly, M = E - esinE
    #(M,e) -> (bE,e) ;Tr2
    #1 = dbEdM - e dbEdM cos(bE) 
    #0 = dbEde - e*dbEde*cos(bE) - sin(bE)
    bE = M+e*np.sin(M)+e**2/2*np.sin(2*M)
    for k in range(0,niter_anomaly):
        bE=bE+(M-bE+e*np.sin(bE))/(1-e*np.cos(bE))

    #The `true anomaly'. With a pi ambiguity,
    #nu = 2*atan(sqrt((1+e)/(1-e))*tan(bE/2))
    #(bE,e) -> (nu,e) ;Tr3
    nu=2*np.arctan2(np.sqrt((1+e)/(1-e))*np.sin(bE/2), np.cos(bE/2))

    #Offset for nu
    #(nu,w) -> alpha ;Tr4
    alpha=nu+w

    if do_deriv:
         dbEdM = 1./(1-e*np.cos(bE))
         dbEde =  np.sin(bE)/(1-e*np.cos(bE))
         #Derivatives are now for alpha (just an offset from nu)
         #From mathematica...
         dnude  = np.sin(bE)/(e-1)/np.sqrt((1+e)/(1-e))/(e*np.cos(bE)-1)
         dnudbE = (e-1)*np.sqrt((1+e)/(1-e))/(e*np.cos(bE) - 1)
         dAdp = dnudbE*dbEdM*dMdp
         dAdt = dnudbE*dbEdM*dMdt
         dAde = dnudbE*dbEde + dnude
    
    #Final calculations (square brackets mean trivial):
    #(alpha,e,i) [+a,n] -> (rho,theta) Tr5
    #We have dAd(p,t,e,w), with alpha=A. Only complex for
    #e, where we need:
    #drhode = drhodA.dAde + drhodnu.dAde + drhode
    #dthde = dthdA.dAde + dthde
    #Also, drhodp = drhodA.dAdp etc...

    #For the derivative part of the code
    sqtmp = np.sqrt(np.cos(alpha)**2 + np.sin(alpha)**2*np.cos(i)**2)
    rho=a*(1-e**2)/(1 + e*np.cos(nu))*sqtmp
    theta=np.arctan2(np.sin(alpha)*np.cos(i),np.cos(alpha))+n

    if do_deriv:
        drhodA = a*(e**2-1)/(1+e*np.cos(nu))*np.cos(alpha)*np.sin(alpha)*(np.sin(i))**2/sqtmp
        drhodnu = -a*e*(e**2-1)*np.sin(nu)/(1+e*np.cos(nu))**2*sqtmp
        drhode =  -a*(2*e+(1+e**2)*np.cos(nu))*sqtmp/(1 + e*np.cos(nu))**2
        drhodi = -a*(1-e**2)/(1+e*np.cos(nu))*np.cos(i)*np.sin(i)*(np.sin(alpha))**2/sqtmp
        dthdA = np.cos(i)/(np.cos(alpha))**2/(1+(np.cos(i)*np.tan(alpha))**2)
        dthdi = -np.sin(i)*np.tan(alpha)/(1+(np.cos(i)*np.tan(alpha))**2)
        #[T0,P,a,e,n,w,i]
        drho =  np.array([(drhodA+drhodnu)*dAdt, (drhodA+drhodnu)*dAdp,  rho/a, drhode+(drhodA+drhodnu)*dAde, \
              np.zeros(len(jds)), drhodA*np.pi/180., drhodi*np.pi/180.])
        dth  =  np.array([dthdA*dAdt, dthdA*dAdp, np.zeros(len(jds)), dthdA*dAde, np.ones(len(jds))*np.pi/180., dthdA*np.pi/180., \
              dthdi*np.pi/180.])*180/np.pi
        deriv =  (drho,dth)
    
    #The radial velocity is in units of semi major axis units per time unit.
    #e.g. if a is in km and P is in seconds, then vr is in km/s.
    #if a is in milli-arcsec and P is in days, then vr has to be multiplied by
    #(1 AU in km) / (1 day in s) / (parallax in milli-arcsec)
    # [Note that the old IDL code didn't have the amplitude out the front]
    vr = 2*np.pi*a*np.sin(i)/P/np.sqrt(1 - e**2) * (np.cos(alpha) + e*np.cos(w))

    if do_deriv:
        return rho, theta*180/np.pi, vr, deriv
    else:
        return rho, theta*180/np.pi, vr
    
def plot_orbit(params,jds,rho,rho_sig,theta,theta_sig):
    """Make a pretty orbital plot
    """
    rho_mod,theta_mod,dummy = binary_orbit(params, jds)
    plt.clf()
    w_ao = np.where(theta_sig > 0)[0]
    plt.plot(rho[w_ao]*np.sin(np.radians(theta[w_ao])),rho[w_ao]*np.cos(np.radians(theta[w_ao])),'.')
    w_lunar = np.where(theta_sig < 0)[0]
    for ix in w_lunar:
        midpoint = np.array([rho[ix]*np.sin(np.radians(theta[ix])),rho[ix]*np.cos(np.radians(theta[ix]))])
        segment = np.array([100*np.cos(np.radians(theta[ix])),-100*np.sin(np.radians(theta[ix]))])
        start_pt = midpoint - segment
        end_pt = midpoint + segment
        plt.plot([start_pt[0],end_pt[0]],[start_pt[1],end_pt[1]],'r-')
        plt.plot([midpoint[0],rho_mod[ix]*cos(np.radians(theta_mod[ix]))],[midpoint[1],rho_mod[ix]*sin(np.radians(theta_mod[ix]))])
    mint = np.min([params[0],np.min(jds)])
    maxt = np.max([params[0] + params[1],np.max(jds)])
    times = mint + (maxt-mint)*np.arange(1001)/1e3
    rho_orbit,theta_orbit,dummy = binary_orbit(params, times)
    plt.plot(rho_orbit*np.sin(np.radians(theta_orbit)),rho_orbit*np.cos(np.radians(theta_orbit)))

def leastsq_orbit_fun(params,jds,rho,rho_sig,theta,theta_sig):
    """A function for scipy.optimize.leastsq. Lunar occultation
    measurements can be indicated by theta_sig < 0. These will 
    be placed at the """
    model_rho,model_theta,dummy = binary_orbit(params,jds)
    theta_diff = ((model_theta - theta + 180) % 360)-180
    if (np.sum(model_rho) != np.sum(model_rho)):
        raise UserWarning
    w_ao = np.where(theta_sig > 0)[0]
    retvect = np.append( (model_rho[w_ao] - rho[w_ao])/rho_sig[w_ao], theta_diff[w_ao]/theta_sig[w_ao])
    w_lunar = np.where(theta_sig < 0)[0]
    if len(w_lunar)>0:
        proj = model_rho[w_lunar]*np.cos(np.radians(theta[w_lunar] - model_theta[w_lunar]))
        retvect = np.append(retvect,(proj-rho[w_lunar])/rho_sig[w_lunar])
#    import pdb; pdb.set_trace()
    return retvect
    
def leastsq_orbit_deriv(params,jds,rho,rho_sig,theta,theta_sig):
    """A function returning the derivative for scipy.optimize.leastsq"""
    model_rho,model_theta,dummy,deriv = binary_orbit(params,jds, do_deriv=True)
    w_lunar = np.where(theta_sig < 0)[0]
    if (len(w_lunar)>0):
        import pdb; pdb.set_trace() #Not implemented yet!
    orbit_fun_deriv = np.concatenate((deriv[0]/np.tile(rho_sig,7).reshape(7,len(jds)),
        deriv[1]/np.tile(theta_sig,7).reshape(7,len(jds))),axis=1)
    return orbit_fun_deriv.T
    
def binary_lnprob(params,jds,rho,rho_sig,theta,theta_sig):
    """For use with e.g. emcee. Return chi^2/2.
    """
    if (params[3] > 1):
        return -np.inf
    if (params[3] < 0):
        return -np.inf
    retval = -0.5*np.sum(leastsq_orbit_fun(params,jds,rho,rho_sig,theta,theta_sig)**2)
#    model_rho,model_theta,dummy = binary_orbit(params,jds)
#    #Difference modulo 360
#    theta_diff = ((theta -model_theta + 180) % 360)-180
#    retval = -0.5*(np.sum( ((rho-model_rho)/rho_sig)**2) + np.sum((theta_diff/theta_sig)**2))
    #We really shouldn't have NaNs... but if we do, it would be good to stop here for
    #bugshooting.
    if (retval != retval):
        raise UserWarning
    return retval

if __name__=="__main__":
    t = np.linspace(0,1,100)
    params0 = [0,1,1,0,0,0,0]
    params_out = [0.5, 1, 1, 0.01, 0, 180, 0]
    params_in = [0,1,1,0.01,0,0,0]
    rho_out, theta_out, vr = binary_orbit(params_out, t, niter_anomaly=20)
    rho_in, theta_in, vr = binary_orbit(params_in, t, niter_anomaly=20)
    rho0, theta0, vr = binary_orbit(params0, t, niter_anomaly=20)
    xy_rots = [rho_in*np.cos(np.radians(theta_in)) - rho_out*np.cos(np.radians(theta_out)), \
        rho_in*np.sin(np.radians(theta_in)) - rho_out*np.sin(np.radians(theta_out))]
    xy_rots = np.array(xy_rots)
    xys = np.zeros_like(xy_rots)
    for i in range(len(theta0)):
        mat = [[np.cos(np.radians(theta0[i])), np.sin(np.radians(theta0[i]))],\
            [-np.sin(np.radians(theta0[i])), np.cos(np.radians(theta0[i]))]]
        xys[:,i] = np.dot(mat, xy_rots[:,i])