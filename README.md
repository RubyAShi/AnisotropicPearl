# Example simulations in example.png
This script calcualtes the z-component of an anisotropic Pearl vortex. It takes in the Pearl length in the x-, y-direction, and the simulation height as inputs. 
Figure 346 plots the z-field and figure 138 plots the z-flux of the vortex when using a 3 um IBM SUQID. 
The simulation range is hardcoded as xRange and yRange, both defaulted at 400 microns. Consider going up when simulating Pearl vortices with larger Pearl lengths.
The simulation range should be at least greater than the shorter Pearl length. 
The size of the effective SQUID pickup loop is hardcoded as 3.53 um, which is the effective loop size of an IBM SQUID used by the Moler group at Stanford. 
One could try calculating with a smaller SQUID, but typically the simulated flux is below the electrical noise floor of the SQUID. This is why we used our biggest SQUID to image Pearl vortices. 
Do not change the definition of the Fourier grid. Doing so will yield incorrect results. 
** example.png is a figure in an unpublished manuscript to be submitted soon
** The theory was developed in Kogan's paper Phys. Rev. B 104, 144512
