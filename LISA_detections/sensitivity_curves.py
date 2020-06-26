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

def plot_sensitivity_curve(Sn, frequency_range=np.logspace(-5, 0, 1000) * u.Hz, Tobs=4 * u.yr,
                           shade=True, filepath=None, fig=None, ax=None, color="#18068b"):
    """ 
        Plot the sensitivity curve of a gravitational wave detector
        
        Args:
            Sn              --> [function]       function that computes sensitivity curve
            frequency_range --> [array_like, Hz] array of frequencies at which to evalute the curve
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
    noise_amplitude = np.sqrt(Sn(frequency_range, Tobs))
    
    # plot the curve and shade if needed
    ax.loglog(frequency_range, noise_amplitude, color=color)
    if shade:
        ax.fill_between(frequency_range, 0, noise_amplitude, alpha=0.2, color=color)
    
    # adjust labels, sizes and frequency limits to plot is flush to the edges
    ax.set_xlabel(r'Frequency [Hz]', fontsize=1.5*fs, labelpad=10)
    ax.set_ylabel(r'ASD $[\rm Hz^{-1/2}]$', fontsize=1.5*fs, labelpad=10)
    ax.tick_params(axis='both', which='major', labelsize=1.5*fs)
    ax.set_xlim(min(frequency_range).value, max(frequency_range).value)
    
    if filepath is not None:
        fig.savefig(filepath)
        
    return fig, ax

def plot_binaries_on_sc(Sn, frequencies, strain,
                        weights=None, snr=None, eccentricity=None, Tobs=4 * u.yr,
                        shade=True, filepath=None, cmap="plasma_r", development=False):
    """
        Plot an array of binaries on a gravitational wave detector sensitivity curve
        
        Args:
            Sn               --> [function]             sensitivity curve function
            frequencies      --> [array_like, Hz]       Gravitational wave frequency of binaries
            strain           --> [array_like, unitless] Strain of binary
            weights          --> [array_like, unitless] Weight for each binary
            snr              --> [array_like, unitless] Signal-to-noise ratios of binaries
            eccentricity     --> [array_like, unitless] Eccentricity of binaries, if present use as colour
            Tobs             --> [float, yr]            Observation time (default=LISA mission length)
            shade            --> [boolean]              whether to shade in the area under the plot
            filepath         --> [string]               if not None, save the plot at this file
            cmap             --> [string]               matplotlib colourmap to use
            
        Returns:
            fig              --> [figure]               matplotlib figure
            ax               --> [axis]                 matplotlib axis
    """
    fs = 20
    SNR_CUTOFF = 7
    
    # if no weights are provided then give them uniform weights
    if weights is None:
        weights = np.ones(len(frequencies))
    
    if snr is not None:
        # mask extremely undetectable binaries
        mask = snr > 1e-4
        strain, snr, weights, frequencies = strain[mask], snr[mask], weights[mask], frequencies[mask]
        if eccentricity is not None:
            eccentricity = eccentricity[mask]

        # define detectable binaries as above a certain SNR cutoff
        detectable = snr > SNR_CUTOFF
    
    # define the range of frequencies to plot (extending to lower frequencies if necessary)
    minfpow = int(np.floor(np.log10(frequencies.min().value)))
    frequency_range = np.logspace(min(-5, minfpow), -2, 1000) * u.Hz
    
    binary_amplitude = np.sqrt(2 * Tobs.to(u.s)) * strain
    
    # plot the sensitivity curve
    color = matplotlib.cm.get_cmap(cmap)(1.0)
    fig, ax = plot_sensitivity_curve(Sn, frequency_range, shade=shade, color=color, Tobs=Tobs)
    plt.tight_layout()
    
    # plot the binaries, colour either by eccentricity, snr or nothing
    if eccentricity is not None:
        cvar = eccentricity
        boundaries = np.linspace(0, 1, 11)
        norm = matplotlib.colors.BoundaryNorm(boundaries, plt.cm.plasma.N)
        scatt = ax.scatter(frequencies.value, binary_amplitude, s=weights*10, c=eccentricity,
                           cmap=cmap, norm=norm, label="Undetectable Binaries")
        cbar = plt.colorbar(scatt, ticks=boundaries[::2])
        cbar.set_label(label="Eccentricity", fontsize=1.5*fs, labelpad=15)
        cbar.ax.tick_params(labelsize=fs)
    elif snr is not None:
        cvar = snr
        boundaries = np.logspace(-4, 2, 13)
        norm = matplotlib.colors.BoundaryNorm(boundaries, plt.cm.plasma_r.N)
            
        scatt = ax.scatter(frequencies.value, binary_amplitude, s=weights*10, c=snr,
                           cmap=cmap, norm=norm, label="Undetectable Binaries")
        def fmt(x, pos):
                a, b = '{:0e}'.format(x).split('e')
                b = int(b)
                return r'$10^{{{}}}$'.format(b)
        cbar = plt.colorbar(scatt, ticks=boundaries[::2], format=matplotlib.ticker.FuncFormatter(fmt))
        yval = norm(SNR_CUTOFF) / norm(1e2)
        cbar.ax.plot([-0.1, 1.1], [yval, yval], linewidth=3, color='white')
        cbar.set_label(label="Signal-to-noise Ratio", fontsize=1.5*fs, labelpad=15)
        cbar.ax.tick_params(labelsize=fs)
    else:
        scatt = ax.scatter(frequencies.value, binary_amplitude, s=weights*100, color=matplotlib.cm.get_cmap(cmap)(0.5))
    
    # pay no attention to the man behind the curtain
    # ignore this section, used for development
    if development:
        # duplicate the y-axis
        ax_right = ax.twinx()
        ax_right.plot(frequencies, binary_amplitude / np.sqrt(Sn(frequencies)), lw=2, color=color)
        ax_right.set_ylabel("SNR Contribution", fontsize=fs)
        ax_right.tick_params(axis="both", which="major", labelsize=0.7 * fs)
        print(np.sum(binary_amplitude / np.sqrt(Sn(frequencies))))
        
        if len(weights) > 1:
            ax.scatter(frequencies[1].value, binary_amplitude[1], s=weights[1]*110, color=matplotlib.cm.get_cmap(cmap)(1.0), 
                       label="GW Freq")

        max_strain = np.where(binary_amplitude == max(binary_amplitude))
        ax.scatter(frequencies[max_strain].value, binary_amplitude[max_strain], s=weights[max_strain]*110,
                   color=matplotlib.cm.get_cmap(cmap)(0.75), label="Max Strain")

        power = binary_amplitude / np.sqrt(Sn(frequencies))
        max_power = np.where(power == max(power))
        print(max_power)
        ax.scatter(frequencies[max_power].value, binary_amplitude[max_power], s=weights[max_power]*110, 
                   color=matplotlib.cm.get_cmap(cmap)(0.25), label="Max SNR")

        ax.legend(loc="lower left", fontsize=0.7 * fs)
    
    if snr is not None:
        ax.scatter(frequencies[detectable], binary_amplitude[detectable], s=150,
                   marker="*", c=cvar[detectable], cmap=cmap, norm=norm, label="Detectable Binaries")
        legend = plt.legend(loc='best', fontsize=fs)
        legend.legendHandles[0].set_sizes([150])
        legend.legendHandles[0].set_facecolor(matplotlib.cm.get_cmap('plasma_r')(0.2))
        legend.legendHandles[1].set_sizes([250])
        legend.legendHandles[1].set_facecolor(matplotlib.cm.get_cmap('plasma_r')(1.0))

    if filepath is not None:
        fig.savefig(filepath, bbox_inches='tight')

    return fig, ax