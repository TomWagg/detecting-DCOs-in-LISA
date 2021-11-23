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
            "long": r"Fixed mass transfer efficiency of $\beta=0.25$",
            "colour": plt.get_cmap("tab20b")(15)
        },
        {
            "file": "massTransferEfficiencyFixed_0_5",
            "short": "C",
            "med": r"$\beta=0.5$",
            "long": r"Fixed mass transfer efficiency of $\beta=0.5$",
            "colour": plt.get_cmap("tab20b")(14)
        },
        {
            "file": "massTransferEfficiencyFixed_0_75",
            "short": "D",
            "med": r"$\beta=0.75$",
            "long": r"Fixed mass transfer efficiency of $\beta=0.75$",
            "colour": plt.get_cmap("tab20b")(13)
        },
        {
            "file": "unstableCaseBB",
            "short": "E",
            "med": "Case BB\nunstable",
            "long": "Case BB mass transfer is always unstable",
            "colour": plt.get_cmap("tab20")(3)
        },
        {
            "file": "unstableCaseBB_opt",
            "short": "F",
            "med": "E + K",
            "long": "Case BB Unstable + Optimistic CE",
            "colour": plt.get_cmap("tab20")(2)
        },
        {
            "file": "alpha0_1",
            "short": "G",
            "med": r"$\alpha_{\rm CE}=0.1$",
            "long": r"CE efficiency $\alpha = 0.1$",
            "colour": plt.get_cmap("tab20b")(11)
        },
        {
            "file": "alpha0_5",
            "short": "H",
            "med": r"$\alpha_{\rm CE}=0.5$",
            "long": r"CE efficiency $\alpha = 0.5$",
            "colour": plt.get_cmap("tab20b")(10)
        },
        {
            "file": "alpha2_0",
            "short": "I",
            "med": r"$\alpha_{\rm CE}=2.0$",
            "long": r"CE efficiency $\alpha = 2.0$",
            "colour": plt.get_cmap("tab20b")(9)
        },
        {
            "file": "alpha10",
            "short": "J",
            "med": r"$\alpha_{\rm CE}=10.0$",
            "long": r"CE efficiency $\alpha = 10.0$",
            "colour": plt.get_cmap("tab20b")(8)
        },
        {
            "file": "optimistic",
            "short": "K",
            "med": "optimistic\nCE",
            "long": "HG donor stars initiating a CE survive CE",
            "colour": plt.get_cmap("tab20")(11)
        },
        {
            "file": "rapid",
            "short": "L",
            "med": "rapid SN",
            "long": "Fryer rapid SN remnant mass prescription",
            "colour": "lightseagreen"
        },
        {
            "file": "maxNSmass2_0",
            "short": "M",
            "med": r"max $m_{\rm NS}$" + "\n" + r"$2.0 \, {\rm M_{\odot}}$",
            "long": r"Maximum NS mass = 2.0 ${\rm M_{\odot}}$",
            "colour": plt.get_cmap("tab20")(1/20)
        },
        {
            "file": "maxNSmass3_0",
            "short": "N",
            "med": r"max $m_{\rm NS}$" + "\n" + r"$3.0 \, {\rm M_{\odot}}$",
            "long": r"Maximum NS mass = 3.0 ${\rm M_{\odot}}$",
            "colour": plt.get_cmap("tab20")(0/20)
        },
        {
            "file": "noPISN",
            "short": "O",
            "med": "no PISN",
            "long": "PISN and pulsational-PISN not implemented",
            "colour": "deepskyblue"
        },
        {
            "file": "ccSNkick_100km_s",
            "short": "P",
            "med": r"$\sigma_{\rm cc}$" + "\n" + r"$100 \, {\rm km s^{-1}}$",
            "long": r"$\sigma_{\rm RMS}^{\rm 1D} = 100 \ {\rm km\ s^{-1}}$ for core-collapse supernova",
            "colour": plt.get_cmap("tab20c")(11/20)
        },
        {
            "file": "ccSNkick_30km_s",
            "short": "Q",
            "med": r"$\sigma_{\rm cc}$" + "\n" + r"$30 \, {\rm km s^{-1}}$",
            "long": r"$\sigma_{\rm RMS}^{\rm 1D} = 30 \ {\rm km\ s^{-1}}$ for core-collapse supernova",
            "colour": plt.get_cmap("tab20c")(8/20)
        },
        {
            "file": "noBHkick",
            "short": "R",
            "med": "no BH\nkicks",
            "long": "Black holes receive not natal kick",
            "colour": "darkgreen"
        },
        {
            "file": "wolf_rayet_multiplier_0_1",
            "short": "S",
            "med": r"$f_{\rm WR} = 0.1$",
            "long": r"Wolf-Rayet wind factor $f_{\rm WR} = 0.1$",
            "colour": "thistle"
        },
        {
            "file": "wolf_rayet_multiplier_5",
            "short": "T",
            "med": r"$f_{\rm WR} = 5$",
            "long": r"Wolf-Rayet wind factor $f_{\rm WR} = 5.0$",
            "colour": "purple"
        }
    ]
)
