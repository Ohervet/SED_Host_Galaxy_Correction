# SED_Host_Galaxy_Correction

Code to remove the elliptical host galaxy emission in IR-to-UV SEDs of AGN. 
This code makes use of PEGASE 2 spectral templates of galaxies [Fioc & Rocca-Volmerange 1999](https://ui.adsabs.harvard.edu/abs/1999astro.ph.12179F/abstract). 

Author: Olivier Hervet (ohervet@ucsc.edu)

Acknowledgments and citation
-------
If you use this code in a scientific publication or communication, please cite the paper [Fioc & Rocca-Volmerange 1999](https://ui.adsabs.harvard.edu/abs/1999astro.ph.12179F/abstract). When this GitHub page goes public, please also mention it (e.g. with a footnote URL).


- Compile the code with:
`make gal`

- Run the code with a parameter file such as:
`./gal galaxy.par`

- This should create a host galaxy spectrum such as:
`data/HGS_13.dat`

- The columns are:
`nu[Hz], F(nu)[erg/s/Hz], nuF(nu)[erg/s], log10(nu)[Hz], log10(nuF(nu))[erg/s]`

- Open `host_gal.py`
This script is not really optimized, so there is some manual work here

- Create a flag with your source name at the beginning of your script
Following the example of other sources, enter the following values associated with your observation details
`Reff = xx` effective radius in arcsec 
`N = xx`    sersic profile index
optical-UV fluxes measured by your instrument in which you want to remove the host

- Then you need to calibrate your host galaxy template to a measured flux of the host in a given filter and given aperture
There is a selection of possible filters already available in the script. If you have the host in another filter you will need to create a new flag with your filter characteristics
You will need to iteratively run ./gal galaxy.par and host_gal.py playing with the parameter galaxy_mass until your template perfectly match the measured flux.
We assume here that you already know the redshift and age of the galaxy. For the age, if you don't have any idea, use 13Gyr as default.
If you know the color (e.g. B-V) of the host galaxy you can also contraint the age and redshift of your host
There is no fitting method implemented here, everthing is done manually.

- Once you calibrated the host, you can retrieve the host-removed fluxes by turning on the flag of your instrument.


## Extra notes
- Be careful on propagating the errors from the measured host flux to the host-removed fluxes.
- None of the galaxy templates are perfectly accurate, systematics uncertainty can be quite large for strongly host-contaminated SED. I would be careful on not overinterpreting results for fluxes being >50% host contaminated. Fluxes >90% contaminated are not safe to remove the host. These fluxes should be considered as upper limits for the non-thermal emission.

