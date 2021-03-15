from matplotlib import pyplot as plt

__all__ = ["variations"]

# object for all physics variations with file name, short label, long label
# and the colour used in plotting with matplotlib
variations = [{
        "file": "fiducial",
        "short": "A",
        "long": "Fiducial",
        "colour": "grey"
    },
    {
        "file": "optimistic",
        "short": "B",
        "long": "Optimistic CE",
        "colour": plt.get_cmap("tab20")(0)
    },
    {
        "file": "alpha0_5",
        "short": "C",
        "long": r"CE efficiency $\alpha = 0.5$",
        "colour": plt.get_cmap("tab20")(3 / 20)
    },
    {
        "file": "alpha2_0",
        "short": "D",
        "long": r"CE efficiency $\alpha = 2.0$",
        "colour": plt.get_cmap("tab20")(2 / 20)
    },
    {
        "file": "unstableCaseBB",
        "short": "E",
        "long": "Unstable case BB mass transfer",
        "colour": plt.get_cmap("tab20")(4 / 20)
    },
    {
        "file": "massTransferEfficiencyFixed_0_25",
        "short": "F",
        "long": r"Mass transfer efficiency fixed $\beta=0.25$",
        "colour": plt.get_cmap("tab20")(16 / 20)
    },
    {
        "file": "massTransferEfficiencyFixed_0_5",
        "short": "G",
        "long": r"Mass transfer efficiency fixed $\beta=0.5$",
        "colour": plt.get_cmap("tab20")(17 / 20)
    },
    {
        "file": "massTransferEfficiencyFixed_0_75",
        "short": "H",
        "long": r"Mass transfer efficiency fixed $\beta=0.75$",
        "colour": plt.get_cmap("tab20")(17 / 20)
    },
    {
        "file": "ccSNkick_100km_s",
        "short": "I",
        "long": "SN kick dispersion " +
                r"$\sigma_{\rm RMS}^{\rm 1D} = 100 \ {\rm km\ s^{-1}}$",
        "colour": plt.get_cmap("tab20")(18 / 20)
    },
    {
        "file": "ccSNkick_30km_s",
        "short": "J",
        "long": "SN kick dispersion " +
                r"$\sigma_{\rm RMS}^{\rm 1D} = 30 \ {\rm km\ s^{-1}}$",
        "colour": plt.get_cmap("tab20")(19 / 20)
    },
    {
        "file": "noBHkick",
        "short": "J",
        "long": "No BH kicks",
        "colour": plt.get_cmap("tab20")(8 / 20)
    },
    {
        "file": "rapid",
        "short": "L",
        "long": "Fryer radid prescription",
        "colour": plt.get_cmap("tab20")(6 / 20)
    },
    {
        "file": "maxNSmass2_0",
        "short": "M",
        "long": r"Maximum neutron star mass = 2.0 ${\rm M_{\odot}}$",
        "colour": plt.get_cmap("tab20")(6 / 20)
    },
    {
        "file": "maxNSmass3_0",
        "short": "N",
        "long": r"Maximum neutron star mass = 3.0 ${\rm M_{\odot}}$",
        "colour": plt.get_cmap("tab20")(6 / 20)
    },
    {
        "file": "noPISN",
        "short": "O",
        "long": "No pair instability supernova",
        "colour": plt.get_cmap("tab20")(6 / 20)
    },
]
