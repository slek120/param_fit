import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def plot(params):
    # Extract parameters
    t1 = 1.
    t2, mu, x1, x2, epsf, V = params

    # Format filename
    filename = 't1_%.2f_t2_%.2f_mu_%.2f_x1_%.2f_x2_%.2f_epsf_%.2f_V_%.2f.png'%(t1,t2,mu,x1,x2,epsf,V)
    
    # Define function to plot
    def energy(k):
        epsk = -2*t1*(np.cos(k[0])+np.cos(k[1]))-4*t2*np.cos(k[0])*np.cos(k[1])-mu
        x0   = -2*x1*(np.cos(k[0])+np.cos(k[1]))-4*x2*np.cos(k[0])*np.cos(k[1])-epsf
        det  = np.sqrt(0.25*(epsk-x0)*(epsk-x0)+V*V)
        return 0.5*(epsk+x0)-det, 0.5*(epsk+x0)+det

    # Set steps
    delta = 0.025
    # Set bounds
    x = np.arange(-3.14, 3.14, delta)
    y = np.arange(-3.14, 3.14, delta)
    # Set values
    X, Y = np.meshgrid(x, y)
    Z1, Z2 = energy([X,Y])

    # Textbox
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    textstr = (
        '$t=%.2e$\n'
        '$t\'=%.2e$\n'
        '$\mu=%.2e$\n'
        '$X_0=%.2e$\n'
        '$X_0\'=%.2e$\n'
        '$\epsilon_f=%.2e$\n'
        '$V=%.2e$'
        )%(t1,t2,mu,x1,x2,epsf,V)


    # Heatmap
    extent = [-3.14, 3.14, -3.14, 3.14]
    plt.clf()
    plt.text(-2.9, 2.9 , textstr, fontsize=14,
        verticalalignment='top', bbox=props)

    plt.imshow(Z1, extent=extent)
    plt.savefig("Ea/heatmaps/"+filename)

    plt.imshow(Z2, extent=extent)
    plt.savefig("Eb/heatmaps/"+filename)

    # Contour plots
    plt.figure()
    CS1 = plt.contour(X, Y, Z1)
    plt.clabel(CS1, inline=1, fontsize=10)
    plt.text(-2.9, 2.9 , textstr, fontsize=14,
        verticalalignment='top', bbox=props)
    plt.savefig("Ea/contours/"+filename)

    plt.figure()
    CS2 = plt.contour(X, Y, Z2)
    plt.clabel(CS2, inline=1, fontsize=10)
    plt.text(-2.9, 2.9 , textstr, fontsize=14,
        verticalalignment='top', bbox=props)
    plt.savefig("Eb/contours/"+filename)

# plot([2.29663137949, 0.00702885711482, 1.80269022673e-07, -7.64627539075e-08, -2.39035812789e-07, -0.000417686911034])