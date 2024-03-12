#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 13:04:57 2024

@author: olivier
"""

import numpy as np
import matplotlib.pyplot as plt
import pathlib
import tempfile
import re
from astropy.io import ascii


#-----physical constants-----#
c = 2.997924e+10 #[cm/s]
h = 4.135667662E-15  # eV * s
pctocm = 3.086e+18 # [cm]




PROGRAM_NAME = "/SED_Host_Galaxy_Correction"
TMP = False
INITIALIZE = False
MAIN_FOLDER = "scripts"
OUPTUT_FOLDER = "results"

def _get_path():
    base_path = str(pathlib.Path().resolve())
    stop = base_path.find(PROGRAM_NAME)
    if stop == -1:
        raise Exception(PROGRAM_NAME + " is not in file path")
    return base_path[:-len(MAIN_FOLDER)] #+ "/"


if TMP:
    TEMP_DIR = tempfile.mkdtemp()
    BASE_PATH = TEMP_DIR + "/"
    FOLDER_PATH = _get_path()
    print(BASE_PATH)
else:
    TEMP_DIR = None
    BASE_PATH = _get_path()
    FOLDER_PATH = _get_path()
    
    
    
    


def read_data(data_file):
    table = ascii.read(FOLDER_PATH + data_file, format='csv',data_start = 1,delimiter='\s')
    return table
    
    
    
    
def read_configs(config_file="host_config.txt", config_string=None, verbose=False):
    """
    Given a relative path to a file, read the host configs.

    See README for format of the configs file.

    This can also be used to parse a string of dictionary values of configs.

    Args:
        config_file (optional): str
            Relative path to file with configurations; default is "mcmc_config.txt"
        config_string (optional): str
            This is a string of the form {'key': value, ...}. If a config_string
            is given, it will be parsed instead of reading from the file. This is
            useful when the values from a previous configuration dictionary are
            read from a file.
        verbose (optional): bool
            Controls whether values are shown; defaults to False

    Returns:
        (dictionary)
        configurations

    Dictionary format:
    "data_file"             (str) relative path to data
    """

    attributes = ["data_file", "output_name", "calibrator", "Reff", "N", "redshift"]
    configurations = {}  # dictionary of parameters
    if config_string is None:
        # read configurations
        with open(FOLDER_PATH + config_file, 'r') as file:
            if verbose:
                print("Reading configuration options from:", config_file)
            for line in file:
                elements = re.split('=|# ', line)
                if elements[0].strip() in attributes:  # determine if line with an attribute
                    configurations[elements[0].strip()] = elements[1].strip()
    else:
        configurations = eval(config_string)

    # check if all configs present
    for att in attributes:
        if att not in configurations:
            raise Exception("No " + att + " provided!")
            
    if config_string is None:       
         # change int params to floats
         configurations["redshift"] = float(configurations["redshift"])
         configurations["Reff"] = float(configurations["Reff"])
         configurations["N"] = float(configurations["N"])

    if verbose:
        # show parameters
        print("Configurations:")
        for att in attributes:
            print("  ", att, "=", configurations[att])

    return configurations



def read_template(template_file):
    
    #first retrieve the frequency table
    lambd = ascii.read(FOLDER_PATH + "templates/freq.dat", format='csv', data_start = 0, delimiter=' ')
    flux = ascii.read(FOLDER_PATH + template_file, format='csv', data_start = 0, delimiter=' ')
    lambd2 = []
    flux2 = []
    for i in range(len(lambd)):
        lambd2 += list(lambd[i])
        flux2 += list(flux[i])
        
    freq = np.flip(c/(np.array(lambd2[:-2])*1.0e-8))
    flux2 = np.flip(flux2[:-2])
    #the last two values of the file are empty
    #need to reverse the frequency array
    return freq,flux2


def MagnToFlux(magn, freq):
    '''
    Convert magnitude in AB system into a flux in erg.cm-2.s-2

    Parameters
    ----------
    magn : float
        AB magnitude
    freq : float
        center frequency of the considered filter [Hz]

    Returns
    -------
    float
        flux in erg.cm-2.s-2
        
    '''
    return 10**((magn + 48.6)/(-2.5)) * freq


def CosmoCalc(z, H0, WM=0.286, WV=0.714):
    '''
    Cosmology calculator developped by Ned Wright
    https://www.astro.ucla.edu/~wright/CosmoCalc.html    
    
    Parameters
    ----------
    z : float
        redshift
    H0 : Hubble constant at z=0
        km s-1 Mpc-1
    WM : float, optional
        ratio density of matter energy. The default is 0.286.
    WV : TYPE, optional
        ratio density of vaccum energy. The default is 0.714.

    Returns
    -------
    float
        distance luminosity [cm]
        age at redshift z   [Gyr]

    '''
    
    # initialize constants
    
    WR = 0.        # Omega(radiation)
    WK = 0.        # Omega curvaturve = 1-Omega(total)
    c = 299792.458 # velocity of light in km/sec
    Tyr = 977.8    # coefficent for converting 1/H into Gyr
    DTT = 0.5      # time from z to now in units of 1/H0
    age = 0.5      # age of Universe in units of 1/H0
    zage = 0.1     # age of Universe at redshift z in units of 1/H0
    zage_Gyr = 0.0 # value of zage in Gyr
    DCMR = 0.0     # comoving radial distance in units of c/H0
    DA = 0.0       # angular size distance
    DL = 0.0       # luminosity distance
    DL_Mpc = 0.0
    a = 1.0        # 1/(1+z), the scale factor of the Universe
    az = 0.5       # 1/(1+z(object))
    
    h = H0/100.
    WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1-WM-WR-WV
    az = 1.0/(1+1.0*z)
    age = 0.
    n=1000         # number of points in integrals
    for i in range(n):
      a = az*(i+0.5)/n
      adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
      age = age + 1./adot
    
    zage = az*age/n
    DTT = 0.0
    DCMR = 0.0

    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)

    DTT = (1.-az)*DTT/n
    DCMR = (1.-az)*DCMR/n
    age = DTT+zage
    zage_Gyr = (Tyr/H0)*zage

    # tangential comoving distance

    ratio = 1.00
    x = np.sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
          ratio =  0.5*(np.exp(x)-np.exp(-x))/x 
        else:
          ratio = np.sin(x)/x
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/6. + y*y/120.
    DCMT = ratio*DCMR
    DA = az*DCMT
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL

    # comoving volume computation

    ratio = 1.00
    x = np.sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
          ratio = (0.125*(np.exp(2.*x)-np.exp(-2.*x))-x/2.)/(x*x*x/3.)
        else:
          ratio = (x/2. - np.sin(2.*x)/4.)/(x*x*x/3.)
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/5. + (2./105.)*y*y

    return DL_Mpc * 1.0e6 * pctocm, zage_Gyr



def host_average(nu_min, nu_max, host_spectrum):
    '''
    Calculate the flux of the host galaxy template in a given filter
    
    Parameters
    ----------
    nu_min : float
        lower bound of the filter bandpass FWHM [Hz]
    nu_max : float
        higher bound of the filter bandpass FWHM [Hz]
    Host_spectrum : 2D np.array
        spectrum from galaxy host template

    Returns
    -------
    float
        average flux over the spectral band

    '''
    flux = []
    for i in range(len(host_spectrum[0])):
        if (nu_max >= host_spectrum[0][i] >= nu_min):
            flux.append(host_spectrum[1][i])
    
    return np.mean(flux)


def bisection(f, a, b, tol=1e-5, dbg=False):
    """Simple root-finding bisection algorithm."""
    if f(a) == 0:
        return a
    if f(b) == 0:
        return b
    if (np.sign(f(a)) == np.sign(f(b))):
        print("interval must bracket sign change")
        return np.nan

    fa = f(a)
    ii = 0 # keep track of iteration count
    while abs(b - a) > tol*abs(b):
        ii = ii + 1
        c = (a + b)/2
        fc = f(c)
        if np.sign(fc) == np.sign(fa):
            a = c
            fa = fc
        else:
            b = c
        if dbg: # optionally look at intermediate values
            print(f"{ii:5.0f} {c:19.15f} {fc:23.15e}")

    return c


def YOUNG_Sersic(a, ae, N=4):
    #Integrated relative luminosity within a circle (YOUNG 1976)
    #galaxy Sersic profile, luminosity evolves as r^(1/N)
    #Simple DeVaucouleur profile if N = 4
    # a = float(a)
    # ae = float(ae)
    b = 7.66924944
    S = 0
    T = a/ae
    for i in range (0,8):
        S += b**i * T**(i/N)/(np.math.factorial(i))
    return 1-np.exp(-b*T**(1./N))*S



def make_output(output_name, clean_SED):
    #write SED output file
    with open(BASE_PATH + OUPTUT_FOLDER +"/" + output_name+".out", 'w') as f:
        f.write("!nu(Hz)         F(ergcm-2s-1)   delta_nu(-)     delta_nu(+)     delta_F(-)      delta_F(+)	instrument\n")
        for i in range(len(clean_SED[0])):
            f.write(f'{clean_SED[0][i]:.4e}         {clean_SED[1][i]:.4e}   {clean_SED[2][i]:.4e}     {clean_SED[3][i]:.4e}     {clean_SED[4][i]:.4e}      {clean_SED[5][i]:.4e}	{clean_SED[6][i]+"_host_corrected"}\n')
        
    
            


def plotting(output_name,host_spectrum, data,host_frac, cal_data, input_SED, clean_SED):
    plt.close('all')
    plt.plot(host_spectrum[0],host_spectrum[1],color="0.7",label = "Full host")
    instrument0 =data["instrument"][0]
    plt.plot(host_spectrum[0],host_spectrum[1] * host_frac[0],color="0.1",label = f"Host within {instrument0} aperture")
    plt.errorbar(cal_data[0],cal_data[1], xerr=np.array(cal_data[2],cal_data[3]), yerr = np.array(cal_data[4],cal_data[5]), fmt='o', label = "Calibrator")
    plt.errorbar(input_SED[0], input_SED[1], xerr=[input_SED[2],input_SED[3]], 
                 yerr = [input_SED[4],input_SED[5]], fmt='o', label = "Input SED")
    plt.errorbar(clean_SED[0], clean_SED[1], xerr=[clean_SED[2],clean_SED[3]], 
                 yerr = [clean_SED[4],clean_SED[5]], fmt='o', label = "SED corrected from host")
    plt.xlim([1e13,5e15])
    plt.ylim([2e-13,1.5e-11])
    plt.loglog()
    plt.xlabel(r"$\nu$ [Hz]", fontsize=14)
    plt.ylabel(r"$\nu F_\nu [\mathrm{erg}~\mathrm{cm}^{-2}~\mathrm{s}^{-1}]$", fontsize=14)
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(BASE_PATH + OUPTUT_FOLDER +"/" + output_name+".pdf", format="pdf")