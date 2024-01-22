# SED_Host_Galaxy_Correction

Code to remove the elliptical host galaxy emission in IR-to-UV SEDs of AGN. 
This code makes use of PEGASE 2 spectral templates of galaxies [Fioc & Rocca-Volmerange 1999](https://ui.adsabs.harvard.edu/abs/1999astro.ph.12179F/abstract). 

Author: Olivier Hervet (ohervet@ucsc.edu)


Acknowledgments and citation
-------
If you use this code in a scientific publication or communication, please cite the paper [Fioc & Rocca-Volmerange 1999](https://ui.adsabs.harvard.edu/abs/1999astro.ph.12179F/abstract). When this GitHub page goes public, please also mention it (e.g. with a footnote URL).


## Quick tutorial

- You need to first find a good reference of the Sersic profile of your host galaxy [Wikipedia](https://en.wikipedia.org/wiki/S%C3%A9rsic_profile)  and its flux measurement in a given filter that we will call a "calibrator" in the following. An example of such a reference can be [Nilsson et al. 2003](https://ui.adsabs.harvard.edu/abs/2003A%26A...400...95N/abstract).

- Convert the calibrator flux in [erg cm-2 s-1]. Most of calibrator fluxes are given in magnitude. You can use the function `MagnToFlux` in `scripts/host_utils.py` for some help.

- Find the frequency range of the calibrator. This task may be tricky because, for some reason, the exact filter and passband are only sometimes mentioned in photometry papers. You may have to investigate the instrument technical document or other publications using the same instrument. If you find the filter transmission curve, you can estimate the frequency range as the FWHM of this curve.

- Create a .dat file in `SED_data` where you will set your calibrator SED point and all other SED points that need to be corrected from the host emission.

- You need to use the exact same format as the example file `B3_2250.dat`, the columns are:
`!nu(Hz)         F(ergcm-2s-1)   delta_nu(-)     delta_nu(+)     delta_F(-)      delta_F(+)	aperture(arcsec)	instrument`

- The row where you put the calibrator does not matter. For the calibrator, set the aperture at 0. Be careful that each SED data point may have a different aperture. You may have to look to the technical documentation of your instrument or photometry analysis chain to retrieve the aperture. Note that what I call "aperture" here is the radius of PSF extraction, which is different from the original meaning of a telescope aperture.

- Now you need to create a configuration file `host_config.txt` in the main folder. Here is an example of a configuration file:
```
  # Configuration file for running host galaxy correction

# Data file:
data_file=SED_data/B3_2250.dat

output_file=B3_2250.out


# Host Calibrator:
# it needs to match the name given in the "instrument" column of the SED data file
calibrator=R_NOT_HIRAC


# Host effective radius [arcsec]
Reff=5.7


# Sersic index
# use 4 for the standard Devaucouleur profile
N=4


# Redshift
redshift=0.1187


# Host template file
template_file=templates/EG_13Gy.dat
```
  

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

