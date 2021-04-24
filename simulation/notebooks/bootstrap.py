import numpy as np
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d


def bootstrapped_kde(variable, weights, seeds, ax, bw_adjust=None,
                     bootstraps=200, x_count=500, log_scale=(False, False),
                     color="tab:blue", label=None, **kwargs):
    """Create a bootstrapped weighted KDE plot.

    Parameters
    ----------
    variable : `float/array`
        Variable that you want to make a KDE of.
    weights : 'float/array'
        Weights associated with each variable (see all to 1 for unweighted)
    seeds : `int/array`
        Seeds that make the binaries in COMPAS
    ax : `matplotlib Axis`
        Axis on which to plot
    bw_adjust : `float`, optional
        Factor by which to adjust the bandwidth, by default None
    bootstraps : `int`, optional
        How many bootstraps to do, by default 200
    x_count : `int`, optional
        How many x values to evaluate at, by default 500
    log_scale : `tuple`, optional
        Whether each axis should be log scaled, by default (False, False)
    color : `str`, optional
        Colour for the KDE, by default "tab:blue"
    label : `str`, optional
        Label for the plotted KDE, by default None

    Returns
    -------
    ax : `matplotlib Axis`
        Axis on which KDE is plotted
    """

    # store the KDE values for each bootstrap
    kde_vals = np.zeros((bootstraps, x_count))

    # make everything unique and adjust weights
    _, uni_index, uni_count = np.unique(seeds,
                                        return_index=True, return_counts=True)
    variable = variable[uni_index]
    weights = weights[uni_index] * uni_count

    # record indices to sample from
    indices = np.arange(len(variable))

    # decide on x values to evaluate at (based on log scaling)
    if log_scale[0]:
        print("WARNING: I think this doesn't work")
        x_vals = np.logspace(np.log10(np.min(variable)),
                             np.log10(np.max(variable)), x_count)
    else:
        x_vals = np.linspace(np.min(variable), np.max(variable), x_count)

    # perform bootstrapping
    for i in range(bootstraps):
        # sample indices
        boot_index = np.random.choice(indices, size=len(indices), replace=True)

        # find kde of sample and adjust bandwidth if need be
        kde = gaussian_kde(variable[boot_index], weights=weights[boot_index])
        if bw_adjust is not None:
            kde.set_bandwidth(kde.factor * bw_adjust)

        # evaluate kde
        kde_vals[i] = kde.evaluate(x_vals)

    # calculate 1- and 2- sigma percentiles
    percentiles = np.percentile(kde_vals, [15.89, 84.1, 2.27, 97.725], axis=0)

    # plot uncertainties as filled areas
    ax.fill_between(x_vals, percentiles[2], percentiles[3],
                    alpha=0.15, color=color, **kwargs)
    ax.fill_between(x_vals, percentiles[0], percentiles[1],
                    alpha=0.3, color=color, **kwargs)

    # plot the regular kde
    kde = gaussian_kde(variable, weights=weights)
    if bw_adjust is not None:
        kde.set_bandwidth(kde.factor * bw_adjust)
    ax.plot(x_vals, kde(x_vals), color=color, label=label, **kwargs)

    # adjust scales if needed
    if log_scale[0]:
        ax.set_xscale("log")
    if log_scale[1]:
        ax.set_yscale("log")

    return ax


def bootstrapped_ecdf(variable, weights, seeds, ax,
                      bootstraps=200, normalisation=None, x_count=10000,
                      log_scale=(False, False), color="tab:blue", label=None,
                      **kwargs):
    """Create a bootstrapped weighted ECDF plot.

    Parameters
    ----------
    variable : `float/array`
        Variable that you want to make a ECDF of.
    weights : 'float/array'
        Weights associated with each variable (see all to 1 for unweighted)
    seeds : `int/array`
        Seeds that make the binaries in COMPAS
    ax : `matplotlib Axis`
        Axis on which to plot
    bootstraps : `int`, optional
        How many bootstraps to do, by default 200
    normalisation : `float`, optional
        A value to normalise the CDF to
    x_count : `int`, optional
        How many x values to evaluate at, by default 500
    log_scale : `tuple`, optional
        Whether each axis should be log scaled, by default (False, False)
    color : `str`, optional
        Colour for the ECDF, by default "tab:blue"
    label : `str`, optional
        Label for the plotted ECDF, by default None

    Returns
    -------
    ax : `matplotlib Axis`
        Axis on which ECDF is plotted
    """
    # store the ECDF values for each bootstrap
    ecdf_vals = np.zeros((bootstraps, x_count))

    # make everything unique and adjust weights
    _, uni_index, uni_count = np.unique(seeds,
                                        return_index=True, return_counts=True)
    variable = variable[uni_index]
    weights = weights[uni_index] * uni_count

    # record indices to sample from
    indices = np.arange(len(variable))

    # decide on x values to evaluate at (based on log scaling)
    if log_scale[0]:
        print("WARNING: I think this doesn't work")
        x_vals = np.logspace(np.log10(np.min(variable)),
                             np.log10(np.max(variable)), x_count)
    else:
        x_vals = np.linspace(np.min(variable), np.max(variable), x_count)

    # perform bootstrapping
    for i in range(bootstraps):
        # sample indices
        boot_index = np.random.choice(indices, size=len(indices), replace=True)

        boot_var = variable[boot_index]
        boot_weight = weights[boot_index]

        # create a CDF
        sorted_index = np.argsort(boot_var)
        y_vals = np.cumsum(boot_weight[sorted_index])
        if normalisation is not None:
            y_vals = y_vals / np.sum(boot_weight) * normalisation

        # interpolate the CDF
        func = interp1d(boot_var[sorted_index], y_vals, bounds_error=False,
                        fill_value=(0.0, np.max(y_vals)))

        # evaluate the interpolation
        ecdf_vals[i] = func(x_vals)

    # calculate 1- and 2- sigma percentiles
    percentiles = np.percentile(ecdf_vals, [15.89, 84.1, 2.27, 97.725], axis=0)

    # plot uncertainties as filled areas
    ax.fill_between(x_vals, percentiles[2], percentiles[3],
                    alpha=0.15, color=color, **kwargs)
    ax.fill_between(x_vals, percentiles[0], percentiles[1],
                    alpha=0.3, color=color, **kwargs)

    sorted_index = np.argsort(variable)
    y_vals = np.cumsum(weights[sorted_index])
    if normalisation is not None:
        y_vals = y_vals / np.sum(weights) * normalisation
    ax.plot(variable[sorted_index], y_vals, zorder=10,
            color=color, label=label)

    if log_scale[0]:
        ax.set_xscale("log")
    if log_scale[1]:
        ax.set_yscale("log")

    return ax


def bootstrapped_Z(Z, weights, seeds, ax, bootstraps=200,
                   color="tab:blue", label=None, **kwargs):
    """Create a bootstrapped weighted metallicity histogram plot.

    Parameters
    ----------
    Z : `float/array`
        Metallicity that you want to make a plot of.
    weights : 'float/array'
        Weights associated with each variable (see all to 1 for unweighted)
    seeds : `int/array`
        Seeds that make the binaries in COMPAS
    ax : `matplotlib Axis`
        Axis on which to plot
    bootstraps : `int`, optional
        How many bootstraps to do, by default 200
    color : `str`, optional
        Colour for the plot, by default "tab:blue"
    label : `str`, optional
        Label for the plotted line, by default None

    Returns
    -------
    ax : `matplotlib Axis`
        Axis on which stuff is plotted
    """
    # work metallicity bins
    Z_vals = np.concatenate((np.logspace(-4, np.log10(0.022), 50).round(5),
                             [0.0244,  0.02705, 0.03]))
    inner_bins = np.array([Z_vals[i]
                           + (Z_vals[i+1] - Z_vals[i])
                           / 2 for i in range(len(Z_vals) - 1)])
    Z_bins = np.concatenate(([Z_vals[0]], inner_bins,
                             [Z_vals[-1]]))

    # make everything unique and adjust weights
    _, uni_index, uni_count = np.unique(seeds,
                                        return_index=True, return_counts=True)
    Z = Z[uni_index]
    weights = weights[uni_index] * uni_count

    indices = np.arange(len(Z))
    hist_vals = np.zeros((bootstraps, len(Z_vals)))
    for i in range(len(hist_vals)):

        boot_index = np.random.choice(indices, size=len(indices), replace=True)

        hist_vals[i], _ = np.histogram(Z[boot_index], bins=Z_bins,
                                       weights=weights[boot_index],
                                       density=True)

    hist, _ = np.histogram(Z, bins=Z_bins, weights=weights, density=True)
    nonzero = hist > 0

    # calculate 1- and 2- sigma percentiles
    percentiles = np.percentile(hist_vals, [15.89, 84.1, 2.27, 97.725], axis=0)

    # plot uncertainties as filled areas
    ax.fill_between(Z_vals[nonzero], percentiles[2][nonzero],
                    percentiles[3][nonzero], alpha=0.15, color=color, **kwargs)
    ax.fill_between(Z_vals[nonzero], percentiles[0][nonzero],
                    percentiles[1][nonzero], alpha=0.3, color=color, **kwargs)

    ax.plot(Z_vals[nonzero], hist[nonzero], color=color,
            zorder=-1, label=label, **kwargs)
    ax.scatter(Z_vals[nonzero], hist[nonzero], s=25, color=color, **kwargs)

    return ax
