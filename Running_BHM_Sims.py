# Example of how to implement shear calculations of BAHAMAS Simulations
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck18
from astropy.cosmology import WMAP9
import BHMsim_py

## Instantiate a BAHAMAS object
baham0 = BHMsim_py.BHM_cluster(name = 'GrNm_004.npz', horizontal_axis = 'x', vertical_axis = 'y', zlens = 0.7)

## Plot the magnification for the object
baham0.plot_magnification0(hfov_hor = 5*u.arcmin, hfov_ver = 5*u.arcmin, vmin = -0.2, vmax = 1.5, cmap = 'plasma')


# Plotting many convergence and shear maps for varioius BAHAMAS Simulations
## generate a list of simulation paths
def num_converter(num):
    if num % 10 == num:
        return str(f"0{num}")
    return str(num)

sim_paths = ["GrNm_0{}.npz".format(num_converter(i)) for i in range(1, 21)]
print(sim_paths)

## for each simulation, plot the convergence and shear
for sim_path in sim_paths:
    baham = BHMsim_py.BHM_cluster(name = sim_path, horizontal_axis = 'x', vertical_axis = 'y', zlens = 0.7)
    baham.plot_convergence0_and_shear0(hfov_hor = 20*u.arcmin, hfov_ver = 20*u.arcmin, shearlines = False, convergence0_kwargs={}, shear0_kwargs={'vmin' : -3, 'cmap' : 'inferno'})