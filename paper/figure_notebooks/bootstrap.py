import numpy as np
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d


def bootstrapped_kde(variable, weights, seeds, ax, bw_adjust=None,
                     bootstraps=200, x_count=500, x_min=None, x_max=None, log_scale=(False, False),
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

    if x_min is None:
        x_min = np.min(variable)
    if x_max is None:
        x_max = np.max(variable)

    # decide on x values to evaluate at (based on log scaling)
    if log_scale[0]:
        print("WARNING: I think this doesn't work", variable)
        x_vals = np.logspace(np.log10(x_min), np.log10(x_max), x_count)
    else:
        x_vals = np.linspace(x_min, x_max, x_count)

    sorted_order = np.argsort(seeds)
    sorted_seeds = seeds[sorted_order]

    # perform bootstrapping
    for i in range(bootstraps):
        _, starts, counts = np.unique(sorted_seeds, return_counts=True, return_index=True)
        res = np.split(sorted_order, starts[1:])
        inds = np.array([np.random.choice(r) if len(r) > 1 else r[0] for r in res])

        loop_variable = variable[inds]
        loop_weights = weights[inds] * counts

        # record indices to sample from
        indices = np.arange(len(loop_variable))

        # sample indices
        boot_index = np.random.choice(indices, size=len(indices), replace=True)

        # find kde of sample and adjust bandwidth if need be
        kde = gaussian_kde(loop_variable[boot_index], weights=loop_weights[boot_index])
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
    ax.plot(x_vals, np.median(kde_vals, axis=0), color=color, label=label, **kwargs)

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

    # record indices to sample from
    indices = np.arange(len(variable))

    # decide on x values to evaluate at (based on log scaling)
    if log_scale[0]:
        x_vals = np.logspace(np.log10(np.min(variable)),
                             np.log10(np.max(variable)), x_count)
    else:
        x_vals = np.linspace(np.min(variable), np.max(variable), x_count)

    sorted_order = np.argsort(seeds)
    sorted_seeds = seeds[sorted_order]

    # perform bootstrapping
    for i in range(bootstraps):
        _, starts, counts = np.unique(sorted_seeds, return_counts=True, return_index=True)
        res = np.split(sorted_order, starts[1:])
        inds = np.array([np.random.choice(r) if len(r) > 1 else r[0] for r in res])

        loop_variable = variable[inds]
        loop_weights = weights[inds] * counts

        # record indices to sample from
        indices = np.arange(len(loop_variable))

        # sample indices
        boot_index = np.random.choice(indices, size=len(indices), replace=True)

        boot_var = loop_variable[boot_index]
        boot_weight = loop_weights[boot_index]

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
    ax.fill_between(x_vals, percentiles[2], percentiles[3], alpha=0.15, color=color, **kwargs)
    ax.fill_between(x_vals, percentiles[0], percentiles[1], alpha=0.3, color=color, **kwargs)

    ax.plot(x_vals, np.median(ecdf_vals, axis=0), color=color, label=label, zorder=10)

    if log_scale[0]:
        ax.set_xscale("log")
    if log_scale[1]:
        ax.set_yscale("log")

    return ax
