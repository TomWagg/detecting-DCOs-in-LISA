from matplotlib import pyplot as plt
import numpy as np

__all__ = ["variations"]

# object for all physics variations with file name, short label, long label
# and the colour used in plotting with matplotlib
variations = np.array(
    [
        {
            "file": "fiducial",
            "short": "A",
            "med": "fiducial",
            "long": "Fiducial",
            "colour": "grey"
        },
        {
            "file": "massTransferEfficiencyFixed_0_25",
            "short": "B",
            "med": r"$\beta=0.25$",
            "long": r"Mass transfer efficiency fixed $\beta=0.25$",
            "colour": plt.get_cmap("tab20b")(15)
        },
        {
            "file": "massTransferEfficiencyFixed_0_5",
            "short": "C",
            "med": r"$\beta=0.5$",
            "long": r"Mass transfer efficiency fixed $\beta=0.50$",
            "colour": plt.get_cmap("tab20b")(14)
        },
        {
            "file": "massTransferEfficiencyFixed_0_75",
            "short": "D",
            "med": r"$\beta=0.75$",
            "long": r"Mass transfer efficiency fixed $\beta=0.75$",
            "colour": plt.get_cmap("tab20b")(13)
        },
        {
            "file": "unstableCaseBB",
            "short": "E",
            "med": "unstable\ncase BB",
            "long": "Unstable case BB mass transfer",
            "colour": plt.get_cmap("tab20")(2)
        },
        {
            "file": "alpha0_5",
            "short": "F",
            "med": r"$\alpha_{\rm CE}=0.5$",
            "long": r"CE efficiency $\alpha = 0.5$",
            "colour": plt.get_cmap("tab20b")(11)
        },
        {
            "file": "alpha2_0",
            "short": "G",
            "med": r"$\alpha_{\rm CE}=2.0$",
            "long": r"CE efficiency $\alpha = 2.0$",
            "colour": plt.get_cmap("tab20b")(10)
        },
        {
            "file": "optimistic",
            "short": "H",
            "med": "optimistic\nCE",
            "long": "Optimistic CE",
            "colour": plt.get_cmap("tab20")(4)
        },
        {
            "file": "rapid",
            "short": "I",
            "med": "rapid SN",
            "long": "Fryer rapid prescription",
            "colour": plt.get_cmap("tab20")(18)
        },
        {
            "file": "maxNSmass2_0",
            "short": "J",
            "med": r"max $m_{\rm NS}$" + "\n" + r"$2.0 \, {\rm M_{\odot}}$",
            "long": r"Maximum neutron star mass = 2.0 ${\rm M_{\odot}}$",
            "colour": plt.get_cmap("tab20")(1)
        },
        {
            "file": "maxNSmass3_0",
            "short": "K",
            "med": r"max $m_{\rm NS}$" + "\n" + r"$3.0 \, {\rm M_{\odot}}$",
            "long": r"Maximum neutron star mass = 3.0 ${\rm M_{\odot}}$",
            "colour": plt.get_cmap("tab20")(0)
        },
        {
            "file": "noPISN",
            "short": "L",
            "med": "no PISN",
            "long": "No pair instability supernova",
            "colour": plt.get_cmap("tab10")(4)
        },
        {
            "file": "ccSNkick_100km_s",
            "short": "M",
            "med": r"$\sigma_{\rm cc}$" + "\n" + r"$100 \, {\rm km s^{-1}}$",
            "long": "SN kick dispersion "
                     + r"$\sigma_{\rm RMS}^{\rm 1D} = 100 \ {\rm km\ s^{-1}}$",
            "colour": plt.get_cmap("tab20b")(2)
        },
        {
            "file": "ccSNkick_30km_s",
            "short": "N",
            "med": r"$\sigma_{\rm cc}$" + "\n" + r"$30 \, {\rm km s^{-1}}$",
            "long": "SN kick dispersion "
                     + r"$\sigma_{\rm RMS}^{\rm 1D} = 30 \ {\rm km\ s^{-1}}$",
            "colour": plt.get_cmap("tab20b")(1)
        },
        {
            "file": "noBHkick",
            "short": "O",
            "med": "no BH\nkicks",
            "long": "No BH kicks",
            "colour": plt.get_cmap("tab20b")(18)
        },
    ]
)
