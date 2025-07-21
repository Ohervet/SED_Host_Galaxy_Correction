# SED_Host_Galaxy_Correction

Code to remove the elliptical host galaxy emission in IR-to-UV SEDs of AGN. 
This code makes use of PEGASE 2 spectral templates of galaxies [Fioc & Rocca-Volmerange 1999](https://ui.adsabs.harvard.edu/abs/1999astro.ph.12179F/abstract). 

Author: Olivier Hervet (ohervet@ucsc.edu)

License
-------
The code is licensed under a [BSD-3-Clause License](LICENSE).


Acknowledgments and citation
-------
If you use this code in a scientific publication or communication, please cite the Zenodo version [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10800919.svg)](https://doi.org/10.5281/zenodo.10800919), and the PEGASE 2 paper [Fioc & Rocca-Volmerange 1999](https://ui.adsabs.harvard.edu/abs/1999astro.ph.12179F/abstract).


## Quick tutorial

- You need to first find a good reference of the Sersic profile of your host galaxy ([Wikipedia](https://en.wikipedia.org/wiki/S%C3%A9rsic_profile))  and its flux measurement in a given filter that we will call a "calibrator" in the following. An example of such a reference can be [Nilsson et al. 2003](https://ui.adsabs.harvard.edu/abs/2003A%26A...400...95N/abstract).

- Convert the calibrator flux in [erg cm-2 s-1]. Most of calibrator fluxes are given in magnitude. You can use the function `MagnToFlux` in `scripts/host_utils.py` for some help (works only for AB magnitude system).

- Find the frequency range of the calibrator. This task may be tricky because, for some reason, the exact filter and passband are only sometimes mentioned in photometry papers. You may have to investigate the instrument technical document or other publications using the same instrument. If you find the filter transmission curve, you can estimate the frequency range as the FWHM of this curve. Multiple filters' passbands are available through the [2SVO database](https://svo2.cab.inta-csic.es/theory/fps/index.php?&mode=browse&gname=FLWO&gname2=KeplerCam&zoom=1&all=0). 

- Create a .dat file in `SED_data` where you will set your calibrator SED point and all other SED points that need to be corrected from the host emission.

- You need to use the exact same format as the example file `B3_2250.dat`, the columns are:
`!nu(Hz)         F(ergcm-2s-1)   delta_nu(-)     delta_nu(+)     delta_F(-)      delta_F(+)	aperture(arcsec)	instrument`

- The row where you put the calibrator does not matter. For the calibrator, set the aperture at 0. Be careful that each SED data point may have a different aperture. You may have to look to the technical documentation of your instrument or photometry analysis chain to retrieve the aperture. Note that what I call "aperture" here is the radius of PSF extraction, which is different from the original meaning of a telescope aperture.

- Now you need to create a configuration file `host_config.txt` in the main folder. Here is an example of a configuration file:
```
  # Configuration file for running host galaxy correction

# Data file:
data_file=SED_data/B3_2250.dat

# Output name:
output_name=B3_2250


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
```

- Launch `scripts/host_corr.py` to apply the host correction. Corrected SED and associated SED plot should be written in `results/`


## Extra notes
- You cannot use this tool if your host galaxy is not an elliptical galaxy.
- None of the galaxy templates are perfectly accurate, and systematic uncertainty can be quite large for strongly host-contaminated SED. You should be mindful of not overinterpreting the results, especially for host corrections >50%. Fluxes >95% contaminated are not safe to be host-removed. These fluxes should be considered as upper limits for the non-thermal emission.

