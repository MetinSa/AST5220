import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# Path for figures to be saved
savepath = "../milestone2/figures/"

# Reading in files
x_rec, a_rec, z_rec, X_e, n_e, tau, tau2, tau22, g, g2, g22 = np.loadtxt("milestone2.dat", unpack = True)

# Latex plotting format
plt.rc("text", usetex = True)
plt.rc("font", family = "serif", size = 15)

# Plotting functions
def X_ePlot():
    
    """
    Function that plots X_e for comparison with Callin

    """

    plt.title("Electron Density")
    plt.semilogy(z_rec[np.argwhere(z_rec<1800)[0,0]:], X_e[np.argwhere(z_rec<1800)[0,0]:], color = "royalblue")
    plt.grid(linestyle = "--")
    plt.xlabel("$z$")
    plt.ylabel("$X_e$")
    plt.savefig(savepath + "xe.pdf")
    plt.gca().invert_xaxis()
    plt.show()

#X_ePlot()

def tauPlot():

    """
    Function which plots the optical depth: tau, tau' and tau''

    """

    plt.title("Optical Depth")
    plt.semilogy(x_rec[200:-3], tau22[200:-3], color = "crimson", label = r"$\tau''$")
    plt.semilogy(x_rec[200:], np.abs(tau2[200:]), color = "mediumseagreen", label = r"$\mid{\tau'}\mid$")
    plt.semilogy(x_rec[200:], tau[200:], color = "royalblue", label = r"$\tau$")
    plt.xlabel("$x$")
    plt.ylabel(r"$\tau$")
    plt.legend()
    plt.grid(linestyle = "--")
    plt.savefig(savepath + "tau.pdf")
    plt.show()

#tauPlot()

def gPlot():

    """
    Function which plots the visibility function: g, g' and g''

    """

    plt.title("Visibility Function")
    plt.plot(x_rec, g22/300, color = "crimson", label = "$g''$")
    plt.plot(x_rec, g2/10, color = "mediumseagreen", label = "$g'$")
    plt.plot(x_rec, g, color = "royalblue", label = "$g$")
    plt.xlabel("$x$")
    plt.ylabel(r"$g$")
    plt.xlim([-7.5,-6])
    plt.legend()
    plt.grid(linestyle = "--")
    plt.savefig(savepath + "g.pdf")
    plt.show()

#gPlot()
