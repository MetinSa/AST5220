import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from cycler import cycler

# Constants
c = 2.99792458e8
Mpc = 3.08568025e22
H0 = (0.7*100* (1e3))/Mpc

# Computing kc/H0 for the chosen k values
k_min = (0.1*H0)/c
k_max = (1000*H0)/c
nk = 100
k_plot = [1,10,30,50,80,100]
ks = np.zeros(nk+1)
k_val = []
for i in range(1,nk+1):
    ks[i] = k_min + (k_max - k_min)*((i-1)/(nk-1))**2
    if i in k_plot:
        k_val.append(ks[i])

# Path for figures to be saved
savepath = '../milestone3/figures/'

# Reading in files
quantities = np.loadtxt('milestone3.dat') 

# Distributing the quantities into their respective arrays
N_quantities = 11
N_k = 6
len_quantities = []
delta = []
delta_b = []
v = []
v_b = []
Phi = []
Psi = []
dPhi = []
dPsi = []
Theta0 = []
Theta1 = []

# Making a list with the start and end point of all variables
for i in range(0,N_k+1):
    len_quantities.append(i*1500)

# Distributing x_t
x_t = quantities[len_quantities[0]:len_quantities[1]-1,0]

# Distributing rest of the variables
for i in range(N_k):
    delta.append(quantities[len_quantities[i]:(len_quantities[i+1]-1),1])
    delta_b.append(quantities[len_quantities[i]:(len_quantities[i+1]-1),2])
    v.append(quantities[len_quantities[i]:(len_quantities[i+1]-1),3])
    v_b.append(quantities[len_quantities[i]:(len_quantities[i+1]-1),4])
    Phi.append(quantities[len_quantities[i]:(len_quantities[i+1]-1),5])
    Psi.append(quantities[len_quantities[i]:(len_quantities[i+1]-1),6])
    dPhi.append(quantities[len_quantities[i]:(len_quantities[i+1]-1),7])
    dPsi.append(quantities[len_quantities[i]:(len_quantities[i+1]-1),8])
    Theta0.append(quantities[len_quantities[i]:(len_quantities[i+1]-1),9])
    Theta1.append(quantities[len_quantities[i]:(len_quantities[i+1]-1),10])

# Latex plotting format
plt.rc('text', usetex = True)
plt.rc('font', family = 'serif', size = 14)
default_cycler  = (cycler(color=['sienna', 'darkorchid', 'darkorange', 'mediumseagreen', 'crimson', 'royalblue']
))
plt.rc('axes', prop_cycle=default_cycler)


# Plotting functions
def deltaplot():
    """
    Plotting delta

    """
    for k in range(N_k):
        plt.semilogy(x_t, delta[k], label = r'$kc$/$H_0$ = %.2f' %((k_val[k]*c)/H0) )
 
    plt.xlabel('$x$')
    plt.ylabel(r'$\delta$')
    plt.legend()
    plt.grid(linestyle = '--')
    plt.savefig(savepath + "delta.pdf")
    plt.show()

#deltaplot()

def delta_bplot():
    """
    Plotting delta_b

    """
    for k in range(N_k):
        plt.semilogy(x_t, delta_b[k], label = r'$kc$/$H_0$ = %.2f' %((k_val[k]*c)/H0) )
 
    plt.xlabel('$x$')
    plt.ylabel(r'$\delta_b$')
    plt.legend()
    plt.grid(linestyle = '--')
    plt.savefig(savepath + "delta_b.pdf")
    plt.show()

#delta_bplot()

def vplot():
    """
    Plotting v

    """
    for k in range(N_k):
        plt.plot(x_t, v[k], label = r'$kc$/$H_0$ = %.2f' %((k_val[k]*c)/H0) )
 
    plt.xlabel('$x$')
    plt.ylabel(r'$v$')
    plt.legend()
    plt.grid(linestyle = '--')
    plt.savefig(savepath + "v.pdf")
    plt.show()

#vplot()

def v_bplot():
    """
    Plotting v_b

    """
    for k in range(N_k):
        plt.plot(x_t, v_b[k], label = r'$kc$/$H_0$ = %.2f' %((k_val[k]*c)/H0) )
 
    plt.xlabel('$x$')
    plt.ylabel(r'$v_b$')
    plt.legend()
    plt.grid(linestyle = '--')
    plt.savefig(savepath + "v_b.pdf")
    plt.show()

#v_bplot()

def Phiplot():
    """
    Plotting Phi

    """
    for k in range(N_k):
        plt.plot(x_t, Phi[k], label = r'$kc$/$H_0$ = %.2f' %((k_val[k]*c)/H0) )
 
    plt.xlabel('$x$')
    plt.ylabel(r'$\Phi$')
    plt.legend()
    plt.grid(linestyle = '--')
    plt.savefig(savepath + "Phi.pdf")
    plt.show()

#Phiplot()

def Psiplot():
    """
    Plotting Psi

    """
    for k in range(N_k):
        plt.plot(x_t, Psi[k], label = r'$kc$/$H_0$ = %.2f' %((k_val[k]*c)/H0) )
 
    plt.xlabel('$x$')
    plt.ylabel(r'$\Psi$')
    plt.legend()
    plt.grid(linestyle = '--')
    plt.savefig(savepath + "Psi.pdf")
    plt.show()

#Psiplot()

def Theta0plot():
    """
    Plotting mononpole Theta0

    """
    for k in range(N_k):            
        plt.plot(x_t, Theta0[k], label = r'$kc$/$H_0$ = %.2f' %((k_val[k]*c)/H0) )
 
    plt.xlabel('$x$')
    plt.ylabel(r'$\Theta_0$')
    plt.legend()
    plt.grid(linestyle = '--')
    plt.savefig(savepath + "Theta0.pdf")
    plt.show()

#Theta0plot()

def Theta1plot():
    """
    Plotting dipole Theta1

    """
    for k in range(N_k):        
        plt.plot(x_t, Theta1[k], label = r'$kc$/$H_0$ = %.2f' %((k_val[k]*c)/H0) )
 
    plt.xlabel('$x$')
    plt.ylabel(r'$\Theta_1$')
    plt.legend(loc = "upper left", fontsize = 12)
    plt.grid(linestyle = '--')
    plt.savefig(savepath + "Theta1.pdf")
    plt.show()

#Theta1plot()
