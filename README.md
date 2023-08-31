Plan for the README

# weak-lensing-mock-generation
Takes cosmic simulations from the BAHAMAS project and evaluates the shear at any given point.

To do: if I want to be able to do the reshift thing, where I change the extent of the. map, I have to recognize that the extent or sizing of the BAHAMAS simulation is basically fixed nby hfov. I think hfov is controlled by somethign called r200, and the sizing is like 1.2-2 or something times r200, so I can figure out what teh scale fo teh BAHAMAS maps (in terms of Mpc) based on that. I willkeep this scale for the sigma plot fo bahamas. Then, when I actually go to caluclate shear, THAT's wehn i switch over to using arcmin or whatever for the sizing of the plots, since superbit measures in arcmin (i think 30 arcmin across for one axis?). i do that by scaling with the redshift information. However, what might be a good idea is to expand the plot of the bahamas simulation by just changing the range to be like 3 times hfov instead of just hfov because that should make some empty space around the bahamas plot which will be better for having not rough edges when computing the fourier transform.

For tomorrow, make sure to push my changes more frequently and to work in branches if necessary. practice smaller commits on my practice repository first. meeting with andrew either right before the darksector meeting or right after journal club.