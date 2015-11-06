import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def heatmap(params):
    t1 = 1
    t2, mu, x1, x2, epsf, V = params

    def energy(k):
        epsk = -2*t1*(np.cos(k[0])+np.cos(k[1]))-4*t2*np.cos(k[0])*np.cos(k[1])-mu
        x0   = -2*x1*(np.cos(k[0])+np.cos(k[1]))-4*x2*np.cos(k[0])*np.cos(k[1])-epsf
        det  = np.sqrt(0.25*(epsk-x0)*(epsk-x0)+V*V)
        return 0.5*(epsk+x0)-det, 0.5*(epsk+x0)+det

    delta = 0.025
    x = np.arange(-3.14, 3.14, delta)
    y = np.arange(-3.14, 3.14, delta)
    X, Y = np.meshgrid(x, y)
    Z1, Z2 = energy([X,Y])

    extent = [-3.14, 3.14, -3.14, 3.14]
    plt.clf()
    plt.title(str(params))

    plt.imshow(Z1, extent=extent)
    plt.savefig("Ea/heatmaps/"+str(params)+".png")

    plt.imshow(Z2, extent=extent)
    plt.savefig("Eb/heatmaps/"+str(params)+".png")

    plt.figure()
    plt.title(str(params))
    CS1 = plt.contour(X, Y, Z1)
    plt.clabel(CS1, inline=1, fontsize=10)
    plt.savefig("Ea/contours/"+str(params)+".png")

    plt.figure()
    plt.title(str(params))
    CS2 = plt.contour(X, Y, Z2)
    plt.clabel(CS2, inline=1, fontsize=10)
    plt.savefig("Eb/contours/"+str(params)+".png")