import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

savepath = "../milestone1/figures/"

MPc= 3.085e22

# Reading in grids and eta
x_eta, a_eta, z_eta, eta  = np.loadtxt("grids.dat", unpack = True)

# Reading in densities
Omega_m, Omega_b, Omega_r, Omega_lambda  = np.loadtxt("densities.dat", unpack = True)

# Reading hubble parameter
H_x  = np.loadtxt("hubble.dat", unpack = True)

# Reading hubble parameter
x_t, eta_splint  = np.loadtxt("splined.dat", unpack = True)

# General plotting
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = 15)

def etaplot():

    # Plotting eta
    plt.figure(1)
    plt.title("Conformal time")
    plt.semilogy(x_eta, eta, color = "royalblue", label = r"$\eta(x)$")
    plt.semilogy(x_t, eta_splint, color = "darkorange", label = r"splined $\eta(x)$")
    plt.grid(linestyle = "--")
    plt.xlabel("x = log $a$")
    plt.ylabel(r"$\eta$")
    plt.legend()
    plt.savefig(savepath + "eta.pdf")
    plt.show()

#etaplot()

def hubbleplot():

    # Plotting the hubble parameter as a function of x and z
    plt.figure(1)
    plt.title("Hubble parameter as a function of x")
    plt.semilogy(x_eta, H_x, color = "royalblue")
    plt.grid(linestyle = "--")
    plt.xlabel("x = log $a$")
    plt.ylabel(r"H$(x)$")
    plt.savefig(savepath + "H_x.pdf")

    plt.figure(2)
    plt.title("Hubble parameter as a function of z")
    plt.loglog(z_eta, H_x, color = "royalblue")
    plt.grid(linestyle = "--")
    plt.xlabel("z")
    plt.ylabel(r"$H(z) $")
    plt.gca().invert_xaxis()
    plt.savefig(savepath + "H_z.pdf")


    plt.show()

hubbleplot()


def densityplot():

    # Plotting the evolution of the density distribution
    plt.title("Relative Densities")
    plt.grid(linestyle = "--")
    plt.xlabel("x = log $a$")
    plt.ylabel(r"$\Omega$")
    plt.plot(x_eta, Omega_m, label = r"$\Omega_m$", color = "royalblue")
    plt.plot(x_eta, Omega_b, label = r"$\Omega_b$", color = "darkorange")
    plt.plot(x_eta, Omega_r, label = r"$\Omega_r$", color = "mediumseagreen")
    plt.plot(x_eta, Omega_lambda, label = r"$\Omega_\lambda$", color = "crimson")
    plt.legend()
    plt.savefig(savepath + "densities.pdf")
    plt.show()

#densityplot()
