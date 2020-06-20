import numpy as np
import astropy.units as u

def LISA_Sn(f, Tobs=4 * u.yr):
    """ 
        Calculate LISA sensitivity curve Sn(f) using equations from Robson et al. 2019
    
        Args:
            f    --> [array_like, Hz]   Frequency range
            Tobs --> [float, yr]        LISA mission time (default 4 years)
            
        Returns:
            Sn   --> [array_like, 1/Hz] Sensitivity curve 
    """
    L = 2.5e9
    fstar = 19.09e-3
    f = f.to(u.Hz).value
    
    # single link optical metrology noise (Robson+ Eq. 10)
    def Poms(f):
        return (1.5e-11)**2 * (1 + (2e-3 / f)**4)

    # single test mass acceleration noise (Robson+ Eq. 11)
    def Pacc(f):
        return (3e-15)**2 * (1 + (0.4e-3 / f)**2) * (1 + (f / (8e-3))**4)

    # galactic confusion noise (Robson+ Eq. 14)
    def Sc(f, Tobs):
        
        # work out the nearest year to the given time and use parameters from Robson+ 2019
        nearest_year = int(round(Tobs.to(u.yr).value))
        if (nearest_year == 1):
            alpha  = 0.133
            beta   = 243.
            kappa  = 482.
            gamma  = 917.
            fk = 2.58e-3  
        elif (nearest_year == 2):
            alpha  = 0.171
            beta   = 292.
            kappa  = 1020.
            gamma  = 1680.
            fk = 2.15e-3 
        elif (nearest_year == 3):
            alpha  = 0.165
            beta   = 299.
            kappa  = 611.
            gamma  = 1340.
            fk = 1.73e-3  
        else:
            alpha  = 0.138
            beta   = -221.
            kappa  = 521.
            gamma  = 1680.
            fk = 1.13e-3

        return 9e-45 * f**(-7/3.) * np.exp(-f**(alpha) + beta * f * np.sin(kappa * f)) \
                * (1 + np.tanh(gamma * (fk - f)))
    
    # Calculate sensitivity curve
    Sn = 10 / (3.0 * L**2) * (Poms(f) + 4 * Pacc(f) / (2 * np.pi * f)**4) * (1 + 0.6 * (f / fstar)**2) + Sc(f, Tobs)
    
    return Sn / u.Hz


# sensitivity with WD confusion noise and MBHB confusion noise
def muAres_Sn(f):
    """ 
        Calculate muAres sensitivity curve Sn(f) (code provided by Valeriya Korol)

        Args:
            f    --> [array_like, Hz]     Frequency range

        Returns:
            Sn   --> [array_like, 1 / Hz] Sensitivity curve 
    """
    f = f.to(u.Hz).value

    cmks=299792458. # m/s
    Larm=395e9 # m
    
    #acc noise
    Sacc=1e-30*(2.*np.pi*f)**(-4.)*(1.+5e-6/f) 

    #shot noise
    Ssn=2.96e-23 * 30

    #other meas noise
    Somn=2.65e-23

    #transfer frequency
    f0T=cmks/(2.*Larm)
    TT=(1.+(f/(0.41*f0T))**2.)**(1./2.)

    #instrument curve
    rSh=( 5**(1./2.)*    #sky average
         (2./3**(1./2.))*    #geometric factor, 1 if 60deg
         TT*(Somn+Ssn+4.*Sacc)**(1./2.)/Larm )

    # WD noise parameters (from Cornish et al. 2018)
    Amp = 0.5*1.8*1e-44
    alpha = 0.138
    beta = -221.
    kappa = 521.
    gamma = 1680.
    fk = 0.00113

    if isinstance(f, float):
        if f > 0.1:
            Sgal = 0.0
        else:
            Sgal = Amp*f**(-7./3.)* np.exp(-f**alpha + beta * f * np.sin(kappa*f))*(1 + np.tanh(gamma * (fk - f)))
    else:
        Sgal = np.zeros(len(f))
        ff = np.array(f[f<=0.1])
        Sgal[f<=0.1] = Amp*ff**(-7./3.)* np.exp(-ff**alpha + beta * ff * np.sin(kappa*ff))*(1 + np.tanh(gamma * (fk - ff)))

    ff2 = 1e-5
    ff3 = 1e-4

    # WD curve
    Sgal[f<ff2] = 8e-20**2/f[f<ff2] # asymptotic limit
    Sgal[(f<ff3)] = ( 2.65e-17*np.tanh(10*(3e-4-f[(f<ff3)])) )**2 / f[(f<ff3)] # smoothed with hyperbolic tangent


    # MBHB gwb curve
    f_new = f*3.
    f_cut_gwb = 1e-6
    ff_new2 = f_new[f_new<f_cut_gwb]
    ff_new3 = f_new[f_new>=f_cut_gwb]
    hc_MBHB = np.zeros(len(f))
    hc_MBHB[f_new<f_cut_gwb] = 1e-10*np.tanh(0.1*(2e-6-ff_new2))*(ff_new2/2e-6)**-0.1
    hc_MBHB[f_new>=f_cut_gwb] = 1.35e-20*(ff_new3/1e-4)**-1.45


    return rSh**2 + Sgal + hc_MBHB**2/f