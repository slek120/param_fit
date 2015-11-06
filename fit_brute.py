import numpy as np
from scipy.optimize import brute
from scipy.optimize import fmin
from plot import heatmap

# Array of (x,y) points
xdata = []
with open("fit01.dat") as f:
    for line in f:
        line = line.split()
        xdata.append((float(line[0]), float(line[1])))
# Array of values at (x,y)
ydata = [0.0 for x in range(len(xdata))]
# Array of experimental error at (x,y)
sigma = [1.0 for x in range(len(xdata))]
# Ranges for fit parameters
# x0 = [2.29663137951, 0.00702885711719, 2.09847731881e-07, -8.84759981594e-08, -2.99660935836e-07, -0.000417687931814]
# ranges = (slice(x0[i]*0.5, x0[i]*1.5, x0[i]*0.5) for i in range(len(x0)))

ranges = (
    slice( 1.0   , 3.0   , 1.0   ),  # t'_min,   t'_max,   t'_step
    slice( 0.0   , 1.e-2 , 0.2e-3),  # mu_min,   mu_max,   mu_step
    slice( 0.0   , 1.e-6 , 0.2e-7),  # x0_min,   x0_max,   x0_step
    slice(-1.e-7 , 0.0   , 0.2e-8),  # x0'_min,  x0'_max,  x0'_step
    slice(-1.e-6 , 0.0   , 0.2e-7),  # epsf_min, epsf_max, epsf_step
    slice(-1.e-3 , 0.0   , 0.2e-4)   # V_min,    V_max,    V_step
)

# Define the objective function to be minimised
# params ... array holding the values of the fit parameters
# X      ... array holding inputs of observed data
# Y      ... array holding values of observed data
# Err    ... array holding errors of observed data
def func(params, X, Y, Err):
    # Extract current values of fit parameters
    t1   = 1
    t2, mu, x1, x2, epsf, V = params

    # Compute sum of chi-squared values
    chi2a = 0.0
    chi2b = 0.0
    for n in range(len(X)):
        k = X[n]

        epsk = -2*t1*(np.cos(k[0])+np.cos(k[1]))-4*t2*np.cos(k[0])*np.cos(k[1])-mu
        x0   = -2*x1*(np.cos(k[0])+np.cos(k[1]))-4*x2*np.cos(k[0])*np.cos(k[1])-epsf
        det  = np.sqrt(0.25*(epsk-x0)*(epsk-x0)+V*V)
        Ea   = 0.5*(epsk+x0)+det
        Eb   = 0.5*(epsk+x0)-det

        chi2a = chi2a + (Y[n] - Ea)*(Y[n] - Ea)#/(Err[n]*Err[n]) # Uncomment to include error
        chi2b = chi2b + (Y[n] - Eb)*(Y[n] - Eb)#/(Err[n]*Err[n]) # Uncomment to include error
    return min(chi2a,chi2b)

# First do a brute force grid optimization and then find local minimum around grid point
params, chi2, grid, Jout = brute(func, ranges, args=(xdata, ydata, sigma), finish=fmin, full_output=True, disp=True)

output = "[" + str(params[0]) + ", " \
             + str(params[1]) + ", " \
             + str(params[2]) + ", " \
             + str(params[3]) + ", " \
             + str(params[4]) + ", " \
             + str(params[5]) + "]\n"
with open("data.txt", "a") as f:
    f.write(output)

output = "fitt1=" + str(1) + ";\n" + \
         "fitt2=" + str(params[0]) + ";\n" + \
         "fitmu=" + str(params[1]) + ";\n" + \
         "fitx1=" + str(params[2]) + ";\n" + \
         "fitx2=" + str(params[3]) + ";\n" + \
         "fitf="  + str(params[4]) + ";\n" + \
         "fitV="  + str(params[5]) + ";"
output = output.replace("e","*10^")
print output

# Plot energy dispersion
heatmap(params)