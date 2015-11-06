import numpy as np
from scipy.optimize import fmin as simplex # Rename fmin to "simplex"
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
x0   = [1.947, 0.2913, 4.5e-9, -7.87e-9, -6.969e-9, -3.9e-5]
xMin = [x0[i]*0.50 for i in range(len(x0))]
xMax = [x0[i]*1.50 for i in range(len(x0))]

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

# Initial params and chi2
params = list(xMin)
chi2   = func(params, xdata, ydata, sigma)

# Recursively iterate
def iterate(index):
    global chi2, params
    if index >= len(xMin):
        return
    steps = 4
    if xMin[index] == xMax[index]:
        steps = 0
    for j in range(1,steps+1):
        # iterate through each index
        iterate(index+1)
        # Linear interpolate from xMin to xMax
        x0[index] = j/float(steps)*(xMax[index]-xMin[index])+xMin[index]
        # Test for improvement in chi2
        testChi2 = func(x0, xdata, ydata, sigma)
        if (testChi2 < chi2):
            params = list(x0)
            chi2 = testChi2

# Start iteration
iterate(0)
print params, chi2
# Find local minimum with downhill simplex algorithm
params = simplex(func, params, args=(xdata, ydata, sigma))
print params, chi2

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