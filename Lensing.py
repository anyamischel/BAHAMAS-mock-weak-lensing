# Creating a function that calculates the intensity of the source at any requested location using the Gaussian distribution for brightness.
import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck18
from astropy.cosmology import WMAP9
from scipy.fft import fft, ifft, fftfreq, fftshift, fft2, ifft2, ifftshift

def gaussian_source_intensity(beta_x, beta_y, rs, e=0, I_0=1*u.W/(u.arcsec)**2, source_x=0*u.arcsec, source_y=0*u.arcsec):
    ''' 
    Given a location (or locations) on the lens plane and parameters about the source, calculate the light intensity using the
    Gaussian distribution.
    This function is only applicable when the ellipse has vertical and horizontal symmetry.
    
    Parameters
    ----------
    beta_x : integer, float, array, or list
        x-coordinate at which the intensity is calculated
    
    beta_y : integer, float, array, or list
        y-coordinate at which the intensity is calculated
    
    I_0 : integer, float
        Maximum intensity of light
    
    rs : integer, float
        Radius of the source, defined for the elliptical case as half the length of the ellipse in the y-direction
        
    e : integer or float between 0 and 1
        Eccentricity of source
        
    source_x : integer, float
        x-coordinate of the center of the source
    
    source_y : integer, float
        y-coordinate of the center of the source
    
    
    Returns
    -----------
    intensity : integer, float, array, or list
        intensity evaluated at the given point (or points)

    '''

    rx = rs*np.sqrt(1 - e**2)
    ry = rs
    
    intensity = I_0 * np.exp(-(beta_x - source_x)**2 / (2*rx**2) - (beta_y - source_y)**2 / (2*ry**2))
    return intensity

def gaussian_source_intensity_r_theta(r, theta, rs, e=0, I_0=1*u.W/(u.arcsec)**2, source_x=0*u.arcsec, source_y=0*u.arcsec):
    ''' Takes polar arugments for the location the intensity will be evaluated instead of cartesian and calculates the intensity
    
     Parameters
    ----------
    r : float, integer
        the radius or distance from the center (0,0) to the location of intensity evaluation
    
    theta : float, integer
        the angle beteween the positive x-axis and the vector pointing from the origin to the location of intensity evaluation. Defined as CCW positive.
        
    Returns
    -------
    intensity : integer, float, array, or list
        intensity evaluated at the given point (or points)
    '''
    
    beta_x = r*np.cos(theta)
    beta_y = r*np.sin(theta)
    
    return gaussian_source_intensity(beta_x, beta_y, I_0=I_0, rs=rs, e=e, source_x=source_x, source_y=source_y)

def gaussian_source_intensity_rot(beta_x, beta_y, rs, phi=0*u.rad, e=0, I_0=1*u.W/(u.arcsec)**2, source_x=0*u.arcsec, source_y=0*u.arcsec):
    ''' Calculates the light intensity at any given location (or locations) from an elliptical source which does not have vertical or horizontal symmetry
    
    Parameters
    ----------
    phi : float, integer
        the angle between the "horizontal" axis of symmetry and the positive x-axis (maybe I should fix this wording)
        
    Returns
    -------
    intensity : integer, float, array, or list
        intensity evaluated at the given point (or points)
        
    '''
    
    ## initial location of evaluation
    theta0 = np.arctan2(beta_y, beta_x)
    r = np.sqrt(beta_x**2 + beta_y**2)
    
    ## transform location by rotating back by phi  
    beta_x_tf = r*np.cos(theta0 - phi)
    beta_y_tf = r*np.sin(theta0 - phi)
    
    ## initial location of the center
    theta_source = np.arctan2(source_y, source_x)
    r_source = np.sqrt(source_x**2 + source_y**2)
    
    ## transform the location of the center by phi as well
    source_x_tf = r_source*np.cos(theta_source - phi)
    source_y_tf = r_source*np.sin(theta_source - phi)
    
    return gaussian_source_intensity(beta_x = beta_x_tf, beta_y = beta_y_tf, rs=rs, I_0 = I_0, e = e, source_x = source_x_tf, source_y = source_y_tf)


def SIS_deflection_angle(lens_center_x, lens_center_y, ee_r, theta_x, theta_y):
    '''
    Calulates a deflection angle field for a given SIS and location of evaluation on the lens plane (theta). 
    Can be called several times and summed to calculate the total deflection angle field due to multiple lens sources.
    
    Parameters
    ----------
    lens_center_x : int
        x-coordinate of the center of the SIS lens
    
    lens_center_y : int
        y-coordinate of the center of the SIS lens
    
    ee_r : float, int
        the einstein ring radius of the particular lens
    
    theta_x : nxn array
        contains the x coordinates of theta, or the point of evaluation of intensity on the lens plane
    
    theta_y : nxn array
        contains the x coordinates of theta, or the point of evaluation of intensity on the lens plane

    Returns
    -------
    [alpha_x, alpha_y] : 1x2 array, where alpha_x and alpha_y are nxn arrays
        contains the x component of alpha, or the deflection angle, corresponding to each theta on the lens plane
        
    '''
    xi_x = theta_x - lens_center_x
    xi_y = theta_y - lens_center_y
    
    # calculating the components of the deflection angle field
    alpha_x = ee_r*xi_x/np.sqrt(xi_x**2 + xi_y**2)
    alpha_y = ee_r*xi_y/np.sqrt(xi_x**2 + xi_y**2)
    
    return [alpha_x, alpha_y]

def ee_r(sigma_v, dL, dS):
    '''
    Calculates the Einstein ring radius
    
    Paramters
    ---------
    sigma_v : float
        the velocity dispersion, around 1400 km/s for a standard galaxy
    
    dL : float
        the distance from the viewer to the lens
    
    dS : float
        the distance from the viewer to the source
    
    Returns
    -------
    ee_r : float
        the einstein ring radius
    
    '''
    
    ee_r = ((4*np.pi*sigma_v**2/const.c**2 * dL/dS).to_value(''))*u.rad
    return ee_r

def Sigma_crit(dS, dLS, dL):
    '''
    Calculates Sigma crit, which is a constant for which >1 indicates strong lensing
    
    Parameters
    ----------
    dS : float
        distance from the viewer to the source
        
    dLS : float
        distance from the lens to the source
        
    dL : float
        distance from the viewer to the lens
        
    Returns
    Sigma_crit : float
        
    '''
    Sigma_crit = const.c**2/(4*np.pi*const.G) * dS/(dL*dLS)
    return Sigma_crit


# Weak lensing
def convergence(Sigma, Sigma_crit):
    '''
    Calculates the convergence (Kappa)
    
    Parameters
    ----------
    Sigma : 2D array
        projected mass density
    
    Sigma_crit : float
        
    '''
    return (Sigma/Sigma_crit).to_value('')

def shear(convergence, padding=1):
    '''
    Calculates the shear given convergence using FFT and IFFT
    
    Parameters
    ----------
    convergence : 2D array
        the convergence (kappa)
        
    padding : float
        indicates how much zero padding to include when computing the fourier transform
    
    Returns
    -------
    shear1 : 2D array
        the first component of the shear, indicating diagonal stretching/squishing
    
    shear2 : 2D array
        the second component of the shear, indicating vertical/horizontal stretching/squishing
    
    shear : 2D array
        the magnitude of the shear
        
    '''
    
    n = len(convergence)
    edgex = len(convergence)
    edgey = len(convergence[0])
    
    # evaluate the 2D FFT on the convergence
    convergence_ft = fft2(convergence, s = [n*padding, n*padding])
    
    kx = fftfreq(n*padding, 2*edgex/(n*padding))
    ky = fftfreq(n*padding, 2*edgey/(n*padding))

    kx,ky = np.meshgrid(kx, ky, indexing='ij')
    
    kx[0, 0] = 1
    ky[0, 0] = 1

    # solving for the components of the fourier transorm of the shear 
    shear1_ft = (kx**2 - ky**2)/(kx**2 + ky**2) * convergence_ft
    shear2_ft = (2*kx*ky)/(kx**2 + ky**2) * convergence_ft

    # setting the origin point of the shears to be 0 to avoid errors
    shear1_ft[0, 0] = 0
    shear2_ft[0, 0] = 0
    
    # performing the inverse fourier transform to get the components of the shear
    shear1 = (np.real(ifft2(shear1_ft)))[:n, :n]
    shear2 = (np.real(ifft2(shear2_ft)))[:n, :n]
    shear = np.sqrt(shear1**2 + shear2**2)
    
    return shear1, shear2, shear

def draw_shearlines(center_x, center_y, shear1, shear2, a, ax, length, color='white'):
    '''
    Draw lines indicating the direction of the stretching of background sources due to the shear on a plot
    
    Parameters
    ----------
    center_x : 2D array
        the x-coordinate of the center point of these lines
    
    center_y : 2D array
        the y-coordinate of the center point of these lines
    
    shear1 : 2D array
        the fist component of the shear
    
    shear2 : 2D array
        the second component of the shear
    
    a : int
        an "inverse density" that controls the quantity of shear lines. the function plots a shear line at every ath element of the center point array
        
    '''
    
    center_samp_x = center_x[::a, ::a].flatten()
    center_samp_y = center_y[::a, ::a].flatten()

    angle = (0.5*np.arctan2((shear2[::a, ::a]).flatten(), (shear1[::a, ::a]).flatten()))*u.rad

    # Calculate the starting and ending points of the lines
    start_x = center_samp_x - 0.5 * length * np.cos(angle)
    start_y = center_samp_y - 0.5 * length * np.sin(angle)
    end_x = center_samp_x + 0.5 * length * np.cos(angle)
    end_y = center_samp_y + 0.5 * length * np.sin(angle)
    
    # Draw the line
    ax.plot([start_x, end_x], [start_y, end_y], color)