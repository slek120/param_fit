# Parameter Fit
Find parameters for a given fermi surface

# Setup
If not installed, run

    pip install numpy scipy matplotlib


# Extract kx and ky values from fermi surface

- Install GraphClicker
- Open image
- Set coordinates
- Click on the fermi surface
- Export to a dat file

# Find a fit

- Change the data file in one of the fit_TYPE.py files
- Set the initial guess
- Run

```
    python fit_TYPE.py
```

# Output

The parameters will be printed to the console and also saved to data.txt

Contour and heatmap plots will also be generated

A mathematica file is also included where the parameters can be manipulated

Copy and paste the output into the mathematica file