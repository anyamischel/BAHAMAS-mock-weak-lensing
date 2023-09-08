# weak-lensing-mock-generation

## Calculating Shear from BAHAMAS Simulations for Mock Weak Gravitational Lensing Image Generation

### Project Description:
This repository loads pre-run BAHAMAS Simulations, which are large-scale N-body simulations of the universe's baryonic and dark matter content, and evaluates and plots the convergence, shear, and magnification fields.

### Installation
Git clone the repository and ensure that your local machine has all the python packages listed in requirements.txt.

### File Overview & Description
This repository consists of the following main files:
Lensing.py
Simulation.py
BHMSim.py
Running_BHM_Sims.py
Report.pdf
requirements.txt


Along with a couple jupyter notebook versions (.ipynb) of these files:
BHMSim_nb.ipynb
Running_BHM_Sims_nb.ipynb

As well as a directory of the pre-run BAHAMAS Simulations:
BAHAMAS_cutouts

-------------------------

Below are descriptions for the function and contents of each of these files:

Lensing.py
----------
This is a python file containing a multitude of functions related to lensing calculations. Lensing.py includes both some functions for strong lensing and weak lensing. The strong lensing functions include a function to calculate the Einstein ring radius, the deflection angle field due to a Singular Isothermal Sphere (SIS), and additional functions to generate mock galaxy clusters at a chosen location, size, and ellipticity using a Gaussian light intensity distribution. The weak lensing functions include a function to calculate the convergence, a function to calculate $\Sigma_{\mathrm{crit}}$, and a function to calculate the shear components. In this project, only the weak lensing functions were used.

Simulation.py
-------------
This is a python file containing a function related to loading and analyzing the BAHAMAS simulations. This function, called Sigma, calculates the projected mass density of one of the BAHAMAS simulations by putting all of the particles from a given simulation into a 2D histogram along a chosen coordinate axis (XY, YZ, XZ).

BHMsim.py
---------
This is a python file containing the class BHMsim, which is a class for loading the pre-run BAHAMAS simulations and calculating and plotting the convergence, shear, and magnification fields. A more detailed description of each function is listed in section ENTER LATER BECAUSE I DONT KNOW WHAT SECTION YET!!!!!!!!!!!!!!

Running_BHM_Sims.py
-------------------
This is a python file with an example of how to instantiate an object of BHMsim and plot its convergence and shear. More details will be included in the 'Running the Project' section.

Report.pdf
----------
This is a full report of the project, including the motivation for the project, an explanation of lensing theory, a derivation of the main equations used to calculate the shear from the convergence, and an overview of limitation and future improvements for the project.

requirements.txt
----------------
This is a text file containing a list of all the python packages needed in order to run the project.


### Running the Project

### Credits
Thank you to Dr. Andrew Robertson, JPL and Dr. Eric Huff, JPL for their enormous support as co-mentors of this project. Thank you as well to the JPL Education Office for organizing the internship program at JPL, and to Dr. Edward Stone for funding my project.