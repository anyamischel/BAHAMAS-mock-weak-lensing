import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
import Lensing as l
import Simulation as s
from astropy.cosmology import Planck18
from astropy.cosmology import WMAP9
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline

class BHM_cluster:
    def __init__(self, name, BHMpath, horizontal_axis, vertical_axis, zlens, bhm_hfov=10*u.Mpc, cosmo=Planck18):
        
        self.name = name
        self.cosmo = cosmo # Version used for astropy
        self.BHMpath = BHMpath # path of the saved BAHAMAS-mock-weak-lensing repository
        
        # Chosen axes for the BAHAMAS simualtion to look at
        self.horizontal_axis = horizontal_axis
        self.vertical_axis = vertical_axis
        
        # Default redshifts for the lens and source
        self.zlens0 = 0.5
        self.zsource0 = 1.0
        
        # Actual redshift of the lens
        self.zlens = zlens
        
        # Convergence and shear assuming default lens and source redshifts - initialize as None initially
        self.convergence0 = None
        self.shear0 = None ## magnitude of shear
        self.shear1_0 = None ## first component of shear
        self.shear2_0 = None ## second component of shear
        
        # The half of the field-of-view of the BAHAMAS simulation (in distance units, not angular distance)
        self.bhm_hfov = bhm_hfov
        
        # The default half of the field-of-view of the final image, in this case the size of the SuperBIT field-of-view.
        self.SuperBIT_hfov_hor = 15*u.arcmin
        self.SuperBIT_hfov_ver = 10*u.arcmin
      
    
## Calculating convergence and shear functions =========================================

    def calculate_Sigma_crit(self, zsource1, zlens1):
        '''
        Calculates Sigma_crit, which is greater than 1 for strong lensing and less than 1 for weak lensing
        
        Parameters
        ----------
        zsource1 : float
            redshift of a background source
        
        zlens1 : float
            redshift of a foreground lens
            
        Returns
        -------
        Sigma_crit1 : float
            Sigma_crit for the given redshifts
            
        '''
        
        dL1 = self.cosmo.angular_diameter_distance(zsource1) # distance from the viewer to the lens plane
        dS1 = self.cosmo.angular_diameter_distance(zlens1) # distance from the viewer to the source plane
        dLS1 = self.cosmo.angular_diameter_distance_z1z2(zlens1, zsource1) # distance from the lens plane to the source
        
        Sigma_crit1 = l.Sigma_crit(dS1, dLS1, dL1)
        
        return Sigma_crit1

    def calculate_convergence0_map(self):
        '''
        Calculates the baseline convergence (convergence0) everywhere over the size of the BAHAMAS simulation, which assumes a default lens and source redshift as seen in the def __init__ function
        
        Returns
        -------
        convergence0 : 2d array
            the baseline convergence everywhere over the simulation
        
        '''
        
        # calculating Sigma
        Sigma = s.Sigma(self.name, self.BHMpath, self.horizontal_axis, self.vertical_axis, self.bhm_hfov)
        
        # default redshift
        Sigma_crit0 = self.calculate_Sigma_crit(self.zsource0, self.zlens0)
        
        # since this is the default convergence, it should only need to be calculated once and then saved
        convergence0 = l.convergence(Sigma, Sigma_crit0)
        self.convergence0 = convergence0 # save this value

        return convergence0
    
    def get_convergence0_map(self):
        '''
        Checks if the convergence0 is already saved and if not, calculates it
        
        Returns
        -------
        convergence0 : 2d array
            the baseline convergence everywhere over the simulation
        
        '''
        
        # checks if convergence0 is already calculated
        if self.convergence0 is None:
            self.calculate_convergence0_map()
            return self.convergence0
        else:
            return self.convergence0
    
    def calculate_shear0_map(self):
        '''
        Calculates the default shear (shear0) everywhere over the size of the BAHAMAS simulation, which assumes a default lens and source redshift as seen in the def __init__ function
        
        Returns
        -------
        shear1_0 : 2d array
            first component of the shear, representing a diagonal stretching or squashing of the original image
        
        shear2_0 : 2d array
            second component of the shear, representing a horizontal/vertical stretching or squashing of the original image
        
        shear0 : 2d array
            the magnitude of the shear, which is a measure of how strong the stretching/squashing of the original image is
        
        '''
        
        self.shear1_0, self.shear2_0, self.shear0 = l.shear(self.get_convergence0_map())
        return self.shear1_0, self.shear2_0, self.shear0
    
    
    def get_shear0_map(self):
        '''
        Checks if the shear0 is already saved and if not, calculates it
        
        Returns
        -------
        shear1_0 : 2d array
            first component of the shear, representing a diagonal stretching or squashing of the original image
        
        shear2_0 : 2d array
            second component of the shear, representing a horizontal/vertical stretching or squashing of the original image
        
        shear0 : 2d array
            the magnitude of the shear, which is a measure of how strong the stretching/squashing of the original image is
        
        '''
        
        # checks if shear0 is already calculated
        if self.shear0 is None:
            self.calculate_shear0_map()
            return self.shear1_0, self.shear2_0, self.shear0
        else:
            return self.shear1_0, self.shear2_0, self.shear0
        
        
    
    def evaluate_convergence(self, loc_hor, loc_vert, zsource, zlens=None):
        '''
        Evaluates the actual convergence for a given source and lens redshift at a given location(s) by performing a 2D interpolation of convergence0
        
        Parameters
        ----------
        loc_hor : float, 1D array
            the location (or locations) on the horizontal axis at the point of evaluaton of the convergence
        
        loc_vert : float, 1D array
            the location (or locations) on the vertical axis at the point of evaluaton of the convergence
        
        zsource : float
            the true redshift of the background source
        
        zlens : float
            the true redshift of the foreground lens which defaults to the object's zlens if no input is given
        
        Returns
        -------
        convergence_eval : float, 2D array
            the true value of the convergence at the inputted location (or locations)
        
        '''
        
        if zlens is None:
            zlens = self.zlens
                
        # default redshifts
        Sigma_crit0 = self.calculate_Sigma_crit(self.zsource0, self.zlens0)
        
        # actual redshift
        Sigma_crit = self.calculate_Sigma_crit(zsource, zlens)
        
        # adjust actual value of convergence to account for non-default redshifts
        redshift_adjustment = Sigma_crit0/Sigma_crit
        
        convergence = self.get_convergence0_map() * redshift_adjustment
                
        # extrapolating convergence to any value
        
        # makes the "x" and "y" coordinates be consistent with the size of the map of convergence
        extent_img = (((self.bhm_hfov/(self.cosmo.angular_diameter_distance(zlens))).to_value(''))*u.rad).to_value('arcmin')
        x = np.arange(-extent_img, extent_img, 2*extent_img/len(convergence))
        y = np.arange(-extent_img, extent_img, 2*extent_img/len(convergence))
        
        convergence_interp = interpolate.RectBivariateSpline(x, y, convergence, kx=3, ky=3)
        
        # evaluate the convergence at the given locations using interpolation
        convergence_eval = convergence_interp(loc_hor.to_value('arcmin'), loc_vert.to_value('arcmin'))
        
        return convergence_eval
    
        
        
    def evaluate_shear(self, loc_hor, loc_vert, zsource, zlens=None):
        '''
        Evaluates the actual shear for a given source and lens redshift at a given location(s) by performing a 2D interpolation of shear0
        
        Parameters
        ----------
        loc_hor : float, 1D array
            the location (or locations) on the horizontal axis at the point of evaluaton of the shear
        
        loc_vert : float, 1D array
            the location (or locations) on the vertical axis at the point of evaluaton of the shear
        
        zsource : float
            the true redshift of the background source
        
        zlens : float
            the true redshift of the foreground lens which defaults to the object's zlens if no input is given
        
        Returns
        -------
        shear_eval : float, 2D array
            the true value of the magnitude of the shear at the inputted location (or locations)
            
        '''
        
        
        if zlens is None:
            zlens = self.zlens
                
        # default redshifts
        Sigma_crit0 = self.calculate_Sigma_crit(self.zsource0, self.zlens0)
        
        # actual redshift
        Sigma_crit = self.calculate_Sigma_crit(zsource, zlens)
        
        # adjust actual value of shear to account for non-default redshifts
        redshift_adjustment = Sigma_crit0/Sigma_crit
        
        shear = (self.get_shear0_map())[2] * redshift_adjustment
                
        
        # extrapolating shear to any value
        
        # makes the "x" and "y" coordinates be consistent with the size of the map of convergence
        extent_img = (((self.bhm_hfov/(self.cosmo.angular_diameter_distance(zlens))).to_value(''))*u.rad).to_value('arcmin')
        
        ## x and y should represent point on grid 
        x = np.arange(-extent_img, extent_img, 2*extent_img/len(shear))
        y = np.arange(-extent_img, extent_img, 2*extent_img/len(shear))
        
        shear_interp = interpolate.RectBivariateSpline(x, y, shear, kx=3, ky=3)
        
        # evaluate the shear at the given locations using interpolation
        shear_eval = shear_interp(loc_hor.to_value('arcmin'), loc_vert.to_value('arcmin'))
        
        return shear_eval
        
        
## Plotting functions =========================================================
    
    def plot_convergence0(self, fig=None, ax=None, hfov_hor=None, hfov_ver=None, **kwargs):
        ''' 
        Plots the log (base 10) of the convergence0 map
        
        Parameters
        ----------
        ax : axis
            the axis which the convergence0 will be plotted on, which is only important when both convergence0 and shear0 are plotted side by side
        
        hfov_hor : int, float
            half of the horizontal field-of-view for viewing convergence0, which defaults to the SuperBIT field-of-view if no argument is given
        
        hfov_ver : int, float
            half of the vertical field-of-view for viewing convergence0, which defaults to the SuperBIT field-of-view if no argument is given
        
        '''
        
        if hfov_hor is None:
            hfov_hor = self.SuperBIT_hfov_hor
            
        if hfov_ver is None:
            hfov_ver = self.SuperBIT_hfov_ver

        if ax is None:
            fig, ax = plt.subplots()
            ax.set_xlabel(self.horizontal_axis + '-coord (arcmin)')
            ax.set_ylabel(self.vertical_axis + '-coord (arcmin)')
            ax.set_title('Convergence$_0$')
            
        # Since BAHAMAS simulations are performed independent of the distance to the clusters, the scales are in Mpc
        # To convert to an extent which is dependent on the distance to the lens (as used by SuperBIT), the scale will instead be in arcmin
        extent_img = (((self.bhm_hfov/(self.cosmo.angular_diameter_distance(self.zlens))).to_value(''))*u.rad).to_value('arcmin')

        # plots log (base 10) of convergence0
        convergence_plt = ax.imshow((np.log10(self.get_convergence0_map()+1e-5)).T, origin = 'lower', extent = [-extent_img, extent_img, -extent_img, extent_img], **kwargs)
        
        # sets the field of view from the input accordingly
        ax.set_xlim(-hfov_hor.to_value('arcmin'), hfov_hor.to_value('arcmin'))
        ax.set_ylim(-hfov_ver.to_value('arcmin'), hfov_ver.to_value('arcmin'))
        #plt.colorbar(convergence_plt, fraction=0.046, pad=0.04)
        
        # Add colorbar, make sure to specify tick locations to match desired ticklabels
        cbar = fig.colorbar(convergence_plt, ticks=[0, -1, -2, -3, -4, -5], fraction=0.046, pad=0.04)
        cbar.ax.set_yticklabels(['10$^0$', '10$^{-1}$', '10$^{-2}$', '10$^{-3}$', '10$^{-4}$', '10$^{-5}$'])  # vertically oriented colorbar                    
        
    def plot_shear0(self, fig=None, ax=None, hfov_hor=None, hfov_ver=None, shearlines=False, shearlines_count=100, **kwargs):
        '''
        Plots the log (base 10) of the shear0 map
        
        Parameters
        ----------
        ax : axis
            the axis which the shear0 will be plotted on, which is only important when both convergence0 and shear0 are plotted side by side
        
        hfov_hor : int, float
            half of the horizontal field-of-view for viewing shear0, which defaults to the SuperBIT field-of-view if no argument is given
        
        hfov_ver : int, float
            half of the vertical field-of-view for viewing shear0, which defaults to the SuperBIT field-of-view if no argument is given
        
        shearlines : boolean
            if True, plots lines on the shear0 plot which represent the direction in which a background source would be stretched, calculated using the components of shear0. the magnitude of the lines is left constant since the shear0 plot itself describes the magnitude
        
        shearlines_count : int
            determines the "inverse density" of shearlines plotted on shear0. For a higher value of shearlines_count, fewer shearlines are plotted.
        
        '''
        
        if hfov_hor is None:
            hfov_hor = self.SuperBIT_hfov_hor
            
        if hfov_ver is None:
            hfov_ver = self.SuperBIT_hfov_ver

        if ax is None:
            fig, ax = plt.subplots()
            ax.set_xlabel(self.horizontal_axis + '-coord (arcmin)')
            ax.set_ylabel(self.vertical_axis + '-coord (arcmin)')
            ax.set_title('Shear$_0$')
            
        # Since BAHAMAS simulations are performed independent of the distance to the clusters, the scales are in Mpc
        # To convert to an extent which is dependent on the distance to the lens (as used by SuperBIT), the scale will instead be in arcmin
        
        extent_img = (((self.bhm_hfov/(self.cosmo.angular_diameter_distance(self.zlens))).to_value(''))*u.rad).to_value('arcmin')
        
        # plots log (base 10) of shear0
        shear_plt = ax.imshow((np.log10((self.get_shear0_map())[2]+1e-5)).T, origin = 'lower', extent = [-extent_img, extent_img, -extent_img, extent_img], **kwargs)
        
        # sets the field of view from the input accordingly
        ax.set_xlim(-hfov_hor.to_value('arcmin'), hfov_hor.to_value('arcmin'))
        ax.set_ylim(-hfov_ver.to_value('arcmin'), hfov_ver.to_value('arcmin'))
        #plt.colorbar(shear_plt, fraction=0.046, pad=0.04)
       
        # Add colorbar, make sure to specify tick locations to match desired ticklabels
        cbar = fig.colorbar(shear_plt, ticks=[0, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0], fraction=0.046, pad=0.04)
        cbar.ax.set_yticklabels(['10$^0$', '10$^{-0.5}$', '10$^{-1.0}$', '10$^{-1.5}$', '10$^{-2.0}$', '10$^{-2.5}$', '10$^{-3.0}$'])  # vertically oriented colorbar                    
        
    
        # checks if shearlines are turned on and if so, plot them
        if shearlines is True:
            ## x0 and y0 represent a grid for the horizontal and vertical axes over the shear map 
            x0 = np.arange(-extent_img, extent_img, 2*extent_img/len(self.shear0))
            y0 = np.arange(-extent_img, extent_img, 2*extent_img/len(self.shear0))
        
            x,y = (np.meshgrid(x0, y0, indexing = 'ij'))*u.arcmin
            
            l.draw_shearlines(center_x=x, center_y=y, shear1=(self.get_shear0_map())[0], shear2=(self.get_shear0_map())[1], length=1*u.arcmin, a=shearlines_count, ax=ax)


    def plot_convergence0_and_shear0(self, hfov_hor=None, hfov_ver=None, shearlines=False, shearlines_count=100, convergence0_kwargs = {}, shear0_kwargs = {}):
        '''
        Plots the log (base 10) of the convergence side by side
        
        Parameters
        ----------
        hfov_hor : int, float
            half of the horizontal field-of-view for viewing convergence0 and shear0, which defaults to the SuperBIT field-of-view if no argument is given
        
        hfov_ver : int, float
            half of the vertical field-of-view for viewing convergence0 and shear0, which defaults to the SuperBIT field-of-view if no argument is given
        
        shearlines : boolean
            if True, plots lines on the shear0 plot which represent the direction in which a background source would be stretched, calculated using the components of shear0. the magnitude of the lines is left constant since the shear0 plot itself describes the magnitude
        
        shearlines_count : int
            determines the "inverse density" of shearlines plotted on shear0. For a higher value of shearlines_count, fewer shearlines are plotted.
        
        convergence0_kwargs : list of kwargs
            keyword arguments for plotting convergence0
        
        shear0_kwargs : list of kwargs
            keyword arguments for plotting shear0
            
        '''
        
        if hfov_hor is None:
            hfov_hor = self.SuperBIT_hfov_hor
            
        if hfov_ver is None:
            hfov_ver = self.SuperBIT_hfov_ver
            
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,10))
        
        # call the earlier functions to plot convergence0 and shear0 on their respective axes
        self.plot_convergence0(fig=fig, ax=ax1, hfov_hor=hfov_hor, hfov_ver = hfov_ver, **convergence0_kwargs)
        self.plot_shear0(fig=fig, ax=ax2, hfov_hor=hfov_hor, hfov_ver = hfov_ver, shearlines=shearlines, shearlines_count=shearlines_count, **shear0_kwargs)

        # set corresponding title and axis labels
        ax1.set_title('Convergence$_0$ ($\kappa$)')
        ax2.set_title('Shear$_0$ ($\gamma$)')
        ax1.set_xlabel(self.horizontal_axis + '-coord (arcmin)')
        ax2.set_xlabel(self.horizontal_axis + '-coord (arcmin)')
        
        ax1.set_ylabel(self.vertical_axis + '-coord (arcmin)')
        ax2.tick_params(left=False)
        ax2.set(yticklabels=[])
            
        plt.tight_layout()
        
        
    def plot_magnification0(self, ax=None, hfov_hor=None, hfov_ver=None, **kwargs):
        '''
        Plots the log (base 10) of the magnitude of the magnification
        
        Parameters
        ----------
        ax : axis
            the axis magnification0 will be plotted on, only important if it is plotted next to something else
        
        hfov_hor : int, float
            half of the horizontal field-of-view for viewing magnification0, which defaults to the SuperBIT field-of-view if no argument is given
        
        hfov_ver : int, float
            half of the vertical field-of-view for viewing magnification0, which defaults to the SuperBIT field-of-view if no argument is given
        
        '''
        
        if hfov_hor is None:
            hfov_hor = self.SuperBIT_hfov_hor
            
        if hfov_ver is None:
            hfov_ver = self.SuperBIT_hfov_ver

        if ax is None:
            fig, ax = plt.subplots()
            ax.set_xlabel(self.horizontal_axis + '-coord (arcmin)')
            ax.set_ylabel(self.vertical_axis + '-coord (arcmin)')
            ax.set_title('Magnification$_0$')
        
        # checks to make sure if convergence0 and shear0 are already calculated and if not, calculate and save them
        if self.convergence0 is None:
            self.calculate_convergence0_map()
            
        if self.shear0 is None:
            self.calculate_shear0_map()
        
        extent_img = (((self.bhm_hfov/(self.cosmo.angular_diameter_distance(self.zlens))).to_value(''))*u.rad).to_value('arcmin')
        
        # calculate magnification0
        magnification0 = 1/((1-self.get_convergence0_map())**2 - (self.get_shear0_map())[2]**2)
        
        # plot magniication0
        magnification_plt = ax.imshow((np.log10(np.abs(magnification0)+1e-5)).T, origin = 'lower', extent = [-extent_img, extent_img, -extent_img, extent_img], **kwargs)
        
        ax.set_xlim(-hfov_hor.to_value('arcmin'), hfov_hor.to_value('arcmin'))
        ax.set_ylim(-hfov_ver.to_value('arcmin'), hfov_ver.to_value('arcmin'))
        
        # Add colorbar, make sure to specify tick locations to match desired ticklabels
        cbar = fig.colorbar(magnification_plt, ticks=[1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -0.2], fraction=0.046, pad=0.04)
        cbar.ax.set_yticklabels(['10$^{1.4}$', '10$^{1.2}$', '10$^{1.0}$', '10$^{0.8}$', '10$^{0.6}$', '10$^{0.4}$', '10$^{0.2}$', '10$^{0}$', '10$^{-0.2}$'])  # vertically oriented colorbar 
