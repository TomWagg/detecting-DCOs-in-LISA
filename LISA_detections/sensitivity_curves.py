import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib

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
        
        # use parameters from Robson+ 2019 by finding closest number of years
        years = Tobs.to(u.yr).value
        if (years < 0.75):
            alpha  = 0.133
            beta   = 243.
            kappa  = 482.
            gamma  = 917.
            fk = 2.58e-3  
        elif (years < 1.5):
            alpha  = 0.171
            beta   = 292.
            kappa  = 1020.
            gamma  = 1680.
            fk = 2.15e-3 
        elif (years < 3.0):
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
def muAres_Sn(f, Tobs):
    """ 
        Calculate muAres sensitivity curve Sn(f) (code provided by Valeriya Korol)

        Args:
            f    --> [array_like, Hz]     Frequency range
            Tobs --> [array_like, yr]     Length of mission

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

def plot_sensitivity_curve(Sn, frequency_range=np.logspace(-5, 0, 1000) * u.Hz, Tobs=4 * u.yr, 
                           asd=False, shade=True, filepath=None, fig=None, ax=None, color="#18068b"):
    """ 
        Plot the sensitivity curve of a gravitational wave detector
        
        Args:
            Sn              --> [function]       function that computes sensitivity curve
            frequency_range --> [array_like, Hz] array of frequencies at which to evalute the curve
            asd             --> [boolean]        whether to plot amplitude spectral density
                                                     instead of characteristic strain
            shade           --> [boolean]        whether the shade in the area under the plot
            filepath        --> [string]         if not None, save plot to this file
            fig             --> [figure]         if not None, plot directly onto this figure
            ax              --> [axis]           if not None, plot directly onto this axis
            
        Returns:
            fig             --> [figure]         matplotlib figure
            ax              --> [axis]           matplotlib axis
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(1, figsize=(12, 9))
    fs = 20
    
    # work out what the noise amplitude should be
    if asd:
        noise_amplitude = np.sqrt(Sn(frequency_range, Tobs))
        ax.set_ylabel(r'ASD $[\rm Hz^{-1/2}]$', fontsize=1.5*fs, labelpad=10)
    else:
        noise_amplitude = np.sqrt(frequency_range * Sn(frequency_range, Tobs))
        ax.set_ylabel(r'Characteristic Strain', fontsize=1.5*fs, labelpad=10)
    
    # plot the curve and shade if needed
    ax.loglog(frequency_range, noise_amplitude, color=color)
    if shade:
        ax.fill_between(frequency_range, 0, noise_amplitude, alpha=0.2, color=color)
    
    # adjust labels, sizes and frequency limits to plot is flush to the edges
    ax.set_xlabel(r'Frequency [Hz]', fontsize=1.5*fs, labelpad=10)
    ax.tick_params(axis='both', which='major', labelsize=1.5*fs)
    ax.set_xlim(min(frequency_range).value, max(frequency_range).value)
    
    if filepath is not None:
        fig.savefig(filepath)
        
    return fig, ax

def plot_binaries_on_sc(Sn, binary_amplitude, frequencies, 
                        weights=None, cvar=None, cbar_label=None, snr=None,
                        asd=False, shade=True, filepath=None, cmap="plasma_r"):
    """
        Plot an array of binaries on a gravitational wave detector sensitivity curve
        
        Args:
            Sn               --> [function]             sensitivity curve function
            binary_amplitude --> [array_like]           characteristic strain (or ASD if asd=True)
            frequencies      --> [array_like, Hz]       binary gravitational wave frequency
            weights          --> [array_like, unitless] weight for each binary
            cvar             --> [array_like]           variable to colour each binary
            cbar_label       --> [string]               colourbar label
            cmap             --> [string]               colourmap to apply to binaries
            norm             --> [function]             colour normalising function
            snr              --> [array_like, unitless] signal-to-noise ratios of binaries
            asd              --> [boolean]              whether to plot amplitude spectral density
                                                        instead of characteristic strain
            shade            --> [boolean]              whether to shade in the area under the plot
            filepath         --> [string]               if not None, save the plot at this file
            
        Returns:
            fig              --> [figure]               matplotlib figure
            ax               --> [axis]                 matplotlib axis
    """
    
    if snr:
        SNR_CUTOFF = 7
        
        # mask extremely undetectable binaries
        mask = snr > 1e-4
        yvals, snr, frequency = yvals[mask], snr[mask], frequency[mask]
        if e is not None:
            e = e[mask]

        # define detectable binaries as above a certain SNR cutoff
        detectable = snr > SNR_CUTOFF
    
    # define the range of frequencies to plot
    minfpow = int(np.floor(np.log10(frequencies.min().value)))
    frequency_range = np.logspace(min(-5, minfpow), -2, 1000) * u.Hz
        
    # create a new figure
    color = matplotlib.cm.get_cmap(cmap)(1.0)
    fig, ax = plot_sensitivity_curve(Sn, frequency_range, asd=asd, shade=shade, color=color)
    plt.tight_layout()
        
    if weights is None:
        weights = np.ones(len(frequencies))
        
    if cvar:
        if snr is None:
            boundaries = np.linspace(0, 1, 11)
            norm = matplotlib.colors.BoundaryNorm(boundaries, plt.cm.plasma.N)
            cmap = 'plasma'
        else:
            boundaries = np.logspace(-4, 2, 13)
            norm = matplotlib.colors.BoundaryNorm(boundaries, plt.cm.plasma_r.N)
            
        scatt = ax.scatter(frequencies.value, binary_amplitude, s=weights*100, c=cvar,
                           cmap=cmap, norm=norm, label="Undetectable Binaries")
    else:
        scatt = ax.scatter(frequencies.value, binary_amplitude, s=weights*100, label="Undetectable Binaries")
        
    
    if snr:
        if cvar:
            ax.scatter(frequencies[detectable], binary_amplitude[detectable], s=150,
                       marker="*", c=cvar[detectable], cmap=cmap, norm=norm, label="Detectable Binaries")
            legend = plt.legend(loc='best', fontsize=fs)
            legend.legendHandles[0].set_sizes([150])
            legend.legendHandles[0].set_facecolor(matplotlib.cm.get_cmap('plasma_r')(0.2))
            legend.legendHandles[1].set_sizes([250])
            legend.legendHandles[1].set_facecolor(matplotlib.cm.get_cmap('plasma_r')(1.0))
        else:
            ax.scatter(frequencies[detectable], binary_amplitude[detectable], s=150,
                       marker="*", label="Detectable Binaries")
            legend = plt.legend(loc='best', fontsize=fs)
            legend.legendHandles[0].set_sizes([150])
            legend.legendHandles[1].set_sizes([250])
    
    if cvar:
        if snr is not None:
            def fmt(x, pos):
                a, b = '{:0e}'.format(x).split('e')
                b = int(b)
                return r'$10^{{{}}}$'.format(b)
            cbar = plt.colorbar(scatt, ticks=boundaries[::2], format=matplotlib.ticker.FuncFormatter(fmt))
            yval = norm(SNR_CUTOFF) / norm(1e2)
            cbar.ax.plot([-0.1, 1.1], [yval, yval], linewidth=3, color='white')
        else:
            cbar = plt.colorbar(scatt, ticks=boundaries[::2])
        cbar.set_label(label=cbar_label, fontsize=1.5*fs, labelpad=15)
        cbar.ax.tick_params(labelsize=fs)

    if filepath is not None:
        fig.savefig(filepath, bbox_inches='tight')

    return fig, ax