import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck18
from astropy.cosmology import WMAP9

def Sigma(name, BHMpath, horizontal_axis, vertical_axis, bhm_hfov):
    ''' 
    Calculates the Sigma (projected surface mass density) using matter distributions from the BAHAMAS simulations.
    
    Parameters
    ----------
    name : str
        name of the BAHAMAS simulation file
    
    horizontal_axis : str
        'x', 'y', or 'z'. Horizontal axis coordinate
    
    vertical_axis : str
        'x', 'y', or 'z'. Vertical axis coordinate
        
    bhm_hfov : int, float
        half of the field-of-view of the BAHAMAS simulation
    
    Returns
    -------
    Sigma : 2D array
       the mass density projected onto the inputted horizontal and vertical axes
    
    '''
    
    # parameters
    bins = 1200
    h=0.700 # some sort of astronomical factor
    
    # Mapping between axis labels and numerical indices
    axis_mapping = {'x': 0, 'y': 1, 'z': 2}
    
    # Check if horizontal and vertical axes are the same
    if horizontal_axis == vertical_axis:
        raise ValueError("Horizontal and vertical axes cannot be the same.")
      
    # Assign axes indecies through axis_mapping
    horizontal_index = axis_mapping[horizontal_axis]
    vertical_index = axis_mapping[vertical_axis]

    
    # Load the BAHAMAS cluster
    # data = np.load("/Users/amischel/Downloads/BAHAMAS_cutouts" + name)
    data = np.load(BHMpath + '/BAHAMAS_cutouts/' + name)

    dm_mass = (data['dm_mass'])*u.solMass
    bh_mass = data['bh_mass']*u.solMass
    gas_mass = data['gas_mass']*u.solMass
    star_mass = data['star_mass']*u.solMass

    # account for the fact that the BAHAMAS files were intitially loaded in units of 1/h
    CoP = (data['CoP']*h)*u.Mpc

    dm_pos = (data['dm_pos']*h)*u.Mpc
    bh_pos = (data['bh_pos']*h)*u.Mpc
    gas_pos = (data['gas_pos']*h)*u.Mpc
    star_pos = (data['star_pos']*h)*u.Mpc
    
    hfov_unitless = bhm_hfov.to_value('Mpc')

    # generate 2D histograms of positions with masses for the weights
    STARS_weighted, xedges_stars, yedges_stars = np.histogram2d(star_pos[:,horizontal_index]-CoP[horizontal_index], star_pos[:,vertical_index]-CoP[vertical_index],range=[[-hfov_unitless,hfov_unitless],[-hfov_unitless,hfov_unitless]],bins=bins, density=False, weights=star_mass[:])
    GAS_weighted, xedges_gas, yedges_gas = np.histogram2d(gas_pos[:,horizontal_index]-CoP[horizontal_index], gas_pos[:,vertical_index]-CoP[vertical_index],range=[[-hfov_unitless,hfov_unitless],[-hfov_unitless,hfov_unitless]],bins=bins, density=False, weights=gas_mass[:])
    BH_weighted, xedges_bh, yedges_bh = np.histogram2d(bh_pos[:,horizontal_index]-CoP[horizontal_index], bh_pos[:,vertical_index]-CoP[vertical_index],range=[[-hfov_unitless,hfov_unitless],[-hfov_unitless,hfov_unitless]],bins=bins, density=False, weights=bh_mass[:])
    DM_weighted, xedges_dm, yedges_dm = np.histogram2d(dm_pos[:,horizontal_index]-CoP[horizontal_index], dm_pos[:,vertical_index]-CoP[vertical_index],range=[[-hfov_unitless,hfov_unitless],[-hfov_unitless,hfov_unitless]],bins=bins, density=False, weights=dm_mass[:])

    bin_area = (xedges_stars[1] - xedges_stars[0])**2

    # calculate convergence
    Sigma = (STARS_weighted + GAS_weighted + BH_weighted + DM_weighted)/bin_area 
    
    return Sigma