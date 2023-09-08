# weak-lensing-mock-generation

## Calculating Shear from BAHAMAS Simulations for Mock Weak Gravitational Lensing Image Generation

## Project Description:
This repository loads pre-run BAHAMAS Simulations, which are large-scale N-body simulations of the universe's baryonic and dark matter content, and evaluates and plots the convergence, shear, and magnification fields.

## Installation
Git clone the repository and ensure that your local machine has all the python packages listed in requirements.txt.

## File Overview & Description
This repository consists of the following main files:
- `Lensing.py`
- `Simulation.py`
- `BHMSim.py`
- `Running_BHM_Sims.py`
- `Report.pdf`
- `requirements.txt`


Along with a couple jupyter notebook versions (.ipynb) of these files:
BHMSim_nb.ipynb
Running_BHM_Sims_nb.ipynb

As well as a directory of the pre-run BAHAMAS Simulations:
BAHAMAS_cutouts

---

Below are descriptions for the function and contents of each of these files:

`Lensing.py`  
This is a python file containing a multitude of functions related to lensing calculations. This includes functions related to both strong lensing and weak lensing, though the only functions used in this project are weak lensing functions. The strong lensing functions calculate quantities like the Einstein ring radius, deflection angle field, etc., while the weak lensing functions calculate quantities such as $\Sigma_{\mathrm{crit}}$, the convergence, and the shear.

`Simulation.py`  
This is a python file containing a function which processes and analyzes the pre-loaded BAHAMAS simulations.

`BHMsim.py`  
This is a python file containing the class BHMsim, which is a class for loading the pre-run BAHAMAS simulations and calculating and plotting the convergence, shear, and magnification fields.

`Running_BHM_Sims.py`  
This is a python file with an example of how to instantiate an object of BHMsim and plot its convergence and shear. More details will be included in the 'Running the Project' section.

`Report.pdf`  
This is a full report of the project, including the motivation for the project, an explanation of lensing theory, a derivation of the main equations used to calculate the shear from the convergence, and an overview of limitation and future improvements for the project.

`requirements.txt`  
This is a text file containing a list of all the python packages needed in order to run the project.

---
## Important functions

###### Lensing.py
The main functions used for this project within Lensing.py are:

Sigma_crit - calculates the critical surface density, $\Sigma_{\mathrm{crit}}$.

convergence - takes in a Sigma, and Sigma crit, and calculates the convergence $\kappa$.

shear - uses 2D discrete fourier transforms to convert from the convergence to the shear. This fourier relationship is derived in Report.pdf.

draw_shearlines - draws lines representing the direction that background sources would be stretched that can be placed on top of an image of the shear plot.



###### Simulation.py

Sigma - calculates the projected surface density of a BAHAMAS simulation by projecting all the particles onto one coordinate plane and making a 2D histogram to turn the particle's weighted positions into a 2D array.


###### BHMsim.py
calculate_Sigma_crit - calls the Sigma_crit function from Lensing.py and inputs the redshift parameters used to instantiate the BAHAMAS object.

calculate_convergence0_map - calls the Sigma function from Lensing.py and its own calculate_Sigma_crit function to calculate a 2D array of the convergence using default redshifts and stores it. This way, the convergence does not have to be re-calculated every time a new lens or source redshift is chosen.

get_convergence0_map - checks if the convergence0 is already saved and if not, calculates it

calculate_shear0_map - calls shear from Lensing.py and its own saved convergence to calculate the shear components and magnitude.

get_shear0_map - checks if the shear is already saved and if not, calculates it

evaluate_convergence - scales the saved convergence by a $\Sigma_{\mathrm{crit}}$ factor to reflect the correct inputted redshifts and then uses linear intepolation to find the convergence at any point.

evaluate_shear - scales the saved shear by a $\Sigma_{\mathrm{crit}}$ factor to reflect the correct inputted redshifts and then uses linear intepolation to find the shear at any point.

plot_convergence0 - plots the log base 10 of the convergence0 map (plotting the true convergence would only be a difference of scaling the colormap, which will not change the shape or characteristics of the image).

plot_shear0 - plots the log base 10 of the shear0 map (magnitude of the shear vector). Includes the option to plot the shearlines, which represent the direction background sources would be stretched.

plot_convergence0_and_shear0 - plots both the log base 10 of the convergence and shear side by side using separate color maps.

plot_magnification0 - plots the log base 10 of the magnification.





## Running the Project

## Credits
Thank you to Dr. Andrew Robertson, JPL and Dr. Eric Huff, JPL for their enormous support as co-mentors of this project. Thank you as well to the JPL Education Office for organizing the internship program at JPL, and to Dr. Edward Stone for funding my project.