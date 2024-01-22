#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 13:09:02 2024

@author: olivier hervet, ohervet@ucsc.edu
"""

import numpy as np
import host_utils


#cosmological constant
H0 = 69.6



#load configs
configs = host_utils.read_configs()
redshift = configs["redshift"]
calibrator = configs["calibrator"]
Reff = configs["Reff"]
N = configs["N"]

#read data from SED data file
data = host_utils.read_data(configs["data_file"])

#find the calibrator measurements
cal_index = list(data["instrument"]).index(calibrator)
cal_data = list(data[cal_index])

#keep only uncorrected SED data
data.remove_row(cal_index)

#read the host template
template = host_utils.read_template(configs["template_file"])

#distance luminosity of the host [cm]
DL = host_utils.CosmoCalc(redshift, H0)


#host  calibration process
def calibration(M_gal_log, template=template):
    host_spectrum = np.asarray(template) 
    host_spectrum[1] = 10**M_gal_log *  2.997924e+10/(template[0]*1e-08)  * template[1] / (4.0*np.pi*DL**2)
    #flux of the host in the calibrator energy band
    simulated_calib = host_utils.host_average(cal_data[0]-cal_data[2], cal_data[0]+cal_data[3], host_spectrum)
    root = simulated_calib - cal_data[1]
    return root
#fit the host galaxy spectrum to the calibrator, we consider the mass of the galaxy as the only free parameter
M_gal_log = host_utils.bisection(calibration, 7, 13, tol=1e-5, dbg=False)


#observed calibrated host spectrum
host_spectrum = np.asarray(template) 
host_spectrum[1] = 10**M_gal_log *  2.997924e+10/(template[0]*1e-08)  * template[1] / (4.0*np.pi*DL**2)



#host fraction in a all apertures
aperture = np.array(data["aperture(arcsec)"])
host_frac = host_utils.YOUNG_Sersic(aperture,Reff, N)

#calculate host flux in each SED spectral points [erg cm-2 s-1]
data_freq = np.array(data['!nu(Hz)'])
data_flux = np.array(data['F(ergcm-2s-1)'])
data_errfreq_low = np.array(data['delta_nu(-)'])
data_errfreq_high = np.array(data['delta_nu(+)'])
data_errflux_low = np.array(data['delta_F(-)'])
data_errflux_high = np.array(data['delta_F(+)'])
data_instruments = np.array(data["instrument"])
input_SED = data_freq, data_flux, data_errfreq_low, data_errfreq_high, data_errflux_low, data_errflux_high, data_instruments

host_in_SED = np.zeros(len(data_freq))
for i in range(len(data_freq)):
    host_in_SED[i] = host_utils.host_average(data_freq[i]-data_errfreq_low[i], data_freq[i]+data_errfreq_high[i], host_spectrum) * host_frac[i]


#host removed SED
clean_flux = data_flux - host_in_SED
    
#Error propagation from the calibrator and SED flux errors
A_errn = np.ones(len(data_freq))*cal_data[4]
clean_errflux_low = np.sqrt(A_errn**2 + data_errflux_low**2)
A_errp = np.ones(len(data_freq))*cal_data[5]
clean_errflux_high = np.sqrt(A_errp**2 + data_errflux_high**2)

clean_SED = data_freq, clean_flux, data_errfreq_low, data_errfreq_high, clean_errflux_low, clean_errflux_high, data_instruments


#write output file
host_utils.make_output(configs["output_file"], clean_SED)

print(f"Fitted mass of the galaxy: {10**M_gal_log:.4e} Mo\n")
print(f"Host galaxy flux in SED without aperture correction [erg.cm-2.s-1]: \n", host_in_SED/host_frac,"\n")
print(f"Host galaxy flux in SED within instrumental apertures [erg.cm-2.s-1]: \n", host_in_SED ,"\n")
print(f"Percentage of host contamination in SED\n",  (1-(data_flux-host_in_SED)/data_flux)*100)


#plot and save results
host_utils.plotting(configs["output_file"],host_spectrum, data,host_frac, cal_data, input_SED, clean_SED)



