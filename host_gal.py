#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

from pylab import*
import csv
from scipy.interpolate import UnivariateSpline

#-----Flags-----#
#Given filter for the host imaging

Filter = "R_NOT_HIRAC"
#available options are:
# Abraham
# Rc
# R_NOT_HIRAC
# R_Johnson
# Ks_SOFI

#Instrument where we want to remove the host
Instrument = "UVOT"
#available options are:
# UVOT
# ATOM
# WISE
# TMASS

#Name of the host galaxy
Host = "B3_2250"
#available options are:
# Aplib
# PKS0625
# Bllac
# Wcomae
# HESSJ1943
# NGC1275


#speed of light [m s-1]
c = 2.99792458e8

#----------Functions----------#

def convfreq(lambd):
    #angström
    return c/(lambd*1e-10)

def factorielle(arg):
  if arg==0:
    return 1
  return arg*factorielle(arg-1)

# def YOUNG(a,ae):
#     #luminosité relative intégrée dans un cercle (YOUNG 1976)
#     #Loi de Vaucouleur en r 1/4
#     a = float(a)
#     ae = float(ae)
#     b = 7.66924944
#     S = 0
#     T = a/ae
#     for n in range (0,8):
#         S += b**n * T**(n/4.)/(factorielle(n))
#         #print n
#     return 1-exp(-b*T**(1./4.))*S
    
def YOUNG_Sersic(a, ae, N=4):
    #luminosité relative intégrée dans un cercle (YOUNG 1976)
    #Profil de Sersic, generalisation de Vaucouleur en r 1/N
    #égal a de Vaucouleur pour N = 4
    a = float(a)
    ae = float(ae)
    b = 7.66924944
    S = 0
    T = a/ae
    for i in range (0,8):
        S += b**i * T**(i/N)/(factorielle(i))
    return 1-exp(-b*T**(1./N))*S

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
    
def ConvMagnFlux_Ks2MASS(magn):
    #conv en erg.cm-2.s-2
    return 1.387e14*10**((magn +1.85 + 48.6)/(-2.5))

def Host_average(nu_min, nu_max, host_spectrum, filtername):
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
    filtername : str
        name of the filter

    Returns
    -------
    float
        average flux over the spectral band

    '''
    flux = []
    for i in range(len(host_spectrum[0])):
        if (nu_max >= host_spectrum[0][i] >= nu_min):
            flux.append(host_spectrum[1][i])
    
    return mean(flux)


#----------Reading PEGASE spectral file of the host----------#

# fichier = csv.reader(open ("data/HGS_13.dat","r"), delimiter=" ")

# matrix = []
# flux_pond= []
# flux_pond_uvot= []
# flux_moy= []
# flux_pond1= []
# flux_moy1= []
# flux_pond2= []
# flux_moy2= []
# 	
# for line in fichier:
#         matrix.append(line)
        
host_file = loadtxt("data/HGS_13.dat")
freq = host_file[:-6,0]
flux = host_file[:-6,2]
#freq are in decreasing order, reverse the arrays
freq = freq[::-1]
flux = flux[::-1]
host_spectrum = array([freq,flux])






#----------Calibration of the host galaxy in a given filter----------#      

if Filter == "Abraham":
    # R kitt peak (abraham 1991)
    nu_min = 4.212e14
    nu_max = 5.102e14

elif Filter == "Rc":
    #R Cousin (scarpa 2000)
    nu_min = 3.997e14
    nu_max = 5.451e14

elif Filter == "R_Johnson":
    #R Johnson (Pursimo 2002)
    nu_min = convfreq(9200)
    nu_max = convfreq(4800)
    
elif Filter == "R_NOT_HIRAC":
    #https://www.not.iac.es/instruments/filters/filters.php
    nu_min = convfreq(6915)
    nu_max = convfreq(5705)
    
elif Filter == "Ks_SOFI":
    #https://www.eso.org/sci/facilities/lasilla/instruments/sofi/inst/Imaging.html
    nu_min = convfreq(22995)
    nu_max = convfreq(20245) 


mean_flux = Host_average(nu_min, nu_max, host_spectrum, Filter)
print(f"Host flux average (from template) in the {Filter} filter: {mean_flux:2.3e} [erg.cm-2.s-1]\n")

#-----------------------------#
#          Ap Lib
#-----------------------------#
if Host == "Aplib":

    #rayon effectif d'Ap Lib en U ["]
    reffU = 2.55
    
    #rayon effectif d'Ap Lib en B ["]
    reffB = 2.9
    
    #rayon effectif d'Ap Lib en V ["]
    reffV = 5.3 # Valeur assez improbable (Visvanatanh 1977)
    
    #rayon effectif d'Ap Lib en R ["]
    reffRkp = 5.66 #old Abraham 1991 (R kitt peak)
    reffRc = 3.7 #HST scarpa 2000
    reffR = 6.72 #Pursimo 2002
    
    reff = [reffRc, reffB, reffU]
    
    #Flux UVOT non corrigé:
    uvot = [2.6093E-11, 1.9526E-11, 1.3217E-11, 1.0630E-11, 1.1317E-11, 1.1070E-11]
    
    #Flux wise non corrigé
    wise = [2.12e-11, 2.13e-11, 2.51e-11, 2.76e-11] #Vizier
    
    #Flux 2MASS non corrigé
    mass = [2.708e-11, 2.652e-11, 2.372e-11] #5 arcsec Vizier


#-----------------------------#
#          PKS 0625-354
#-----------------------------#

if Host == "PKS0625":
    print("PKS 0625-354")

    #rayon effectif en Rc ["] (Govoni et al. 2000)
    #Reff = 13.08 #avec table 3
        
    #rayon effectif Ks ["] (Inskip et al. 2010)
    Reff = 18.77
    
    #flux de la galaxie hôte en Rc [erg.cm-2.s-1], déduit de (Govoni et al. 2000)
    #F_Rc = 7.366e-11 #table 3
    #F_Rc = 4.332e-11 #table 2
    #F_Rc = 6.354e-11 #table 3 Mhost tot en corrigeant l'absorption par une diff de 0.2 magn (entre Schlafly 2011 et leur article 2000)
    
    #flux de la galaxie hôte en Ks [erg.cm-2.s-1]
    #F_Ks = 5.428e-11   (Inskip et al. 2010) 
    
    #Flux UVOT non corrigés V,B,U,UVW1,UVM2,UVW2 :
    uvot = [1.404722e-11, 8.13053906379e-12, 5.05785149231e-12, 3.99987733815e-12, 4.95892865495e-12, 4.75420812174e-12]
    
    
    #flux high uvot
#    uvot = [1.88751485e-11 + 9.67956332e-13, 1.14083968e-11 + 5.45856305e-13, 6.78772635e-12 + 4.53522919e-13,
#            4.39498354e-12 + 4.04559155e-13, 4.82402005e-12+ 2.96862772e-13, 4.08855756e-12 + 2.66645058e-13]
    
    #Flux ATOM non corrigé R:
    atom = 2.00021599e-11
    

#-----------------------------#
#          Bl Lac
#-----------------------------#

if Host == "Bllac":

    #rayon effectif en Rc ["] (HST, Scarpa et al. 2000)
    Reff = 4.8
    
    #densité de flux de la galaxie hôte en Rc [erg.cm-2.s-1], déduit de (Scarpa et al. 2000)
    F_Rc = 9.09e-12
    
    #Flux UVOT non corrigés V,B,U,UVW1,UVM2,UVW2 (from ASDC) :
    uvot = [3.4277E-11, 2.5942E-11, 1.9953E-11, 1.2823E-11, 1.5668E-11, 1.2388E-11]
    

#-----------------------------#
#          W Comae
#-----------------------------#

if Host == "Wcomae":

    #rayon effectif en Rc ["] (HST, Scarpa et al. 2000)
    Reff = 2.1 
    DReff = 0.4
    
    Reff = Reff - DReff
    
    #densité de flux de la galaxie hôte en Rc [erg.cm-2.s-1], déduit de (Nilsson et al. 2003)
    F_Rc = 3.319e-12
    DF_Rc = 0.305e-12
    
    F_Rc =F_Rc + DF_Rc
    
    #Flux UVOT non corrigés V,B,U,UVW1,UVM2,UVW2 (from ASDC) :  (UVM2 non dispo)
    uvot = [2.06063E-011, 2.07014E-011, 1.93197E-011, 1.77828E-011, 0.0, 1.65577E-011]

 
#-----------------------------#
#          HESS J1943
#-----------------------------#

if Host == "HESSJ1943":

    #rayon effectif max en 2MASS K["] (HST, Scarpa et al. 2000)
    Reff = 2.5 
    N = 8.0
    ouv = 4.0
    
    print("Fraction de luminosite de la galaxy hôte pour une ouverture de 5 arcsec:", YOUNG_Sersic(ouv,Reff, N) , "\n")
    

#-----------------------------#
#          NGC 1275
#-----------------------------#

if Host == "NGC1275":

    #rayon effectif de Spitzer (3.6um) [kpc] (Sani et al. 2018)
    Reff_kpc = 42 #+-15
    
    #Mathews 2006 state an effective radius of 6.41 kpc
    Reff_kpc = 6.41
    
    #using the cosmology of Sani 2018 (H0 =70), we have 0.357kpc/"
    Reff = Reff_kpc/0.357 #["]    
    
    #most reliable source "Third reference catalogue of bright galaxies, De Vaucouleur 1991"
    Reff = 16.9#+2.5 -2.2 ["]
    
    N = 4.0 # mean it follows a de Vaucouleur profile (as stated by Sani 2018)
    #Flux UVOT non corrigés V,B,U,UVW1,UVM2,UVW2 [erg cm-2 s-1]:
    #from NGC1275_Flare2017_SED.dat
    F_uvot = [1.2233E-10, 1.0898E-10, 7.3274E-11, 5.8938E-11, 6.4023E-11, 5.1152E-11]
    print('Non corrected UVOT fluxes [erg cm-2 s-1]:',F_uvot)
    

if Host == "B3_2250":
    #Nilsson et al 2003 (https://ui.adsabs.harvard.edu/abs/2003A%26A...400...95N/abstract)
    
    #this is for the Devaucouleur solution
    Reff = 5.7 #± 0.2 ["]
    N = 4
    m_Host = 16.04 #in R_NOT_HIRAC filter
    m_Host_err = 0.03
    
    # #this is for the Sersic (free beta) solution
    # Reff = 6.8 #± 0.4 ["]
    # Beta = 0.21 #± 0.01
    # N = 1./Beta
    # m_Host = 15.92 #in R_NOT_HIRAC filter
    # m_Host_err = 0.05
    
    #Host contaminated UVOT fluxes (as Given by Atreya on Nov 22, 2023) [erg cm-2 s-1] 
    #in V,B,U,UVW1,UVM2,UVW2 bands
    F = array([5.49234345e-12, 3.55839062e-12, 2.73971192e-12, 2.35877958e-12, 2.53642188e-12, 2.08358131e-12])
    F_err = array([3.74659399e-13, 2.25931732e-13, 1.83323751e-13, 1.81832172e-13, 2.01758742e-13, 1.42501809e-13])

nu_center = (nu_min + nu_max)/2.
flux = MagnToFlux(m_Host, nu_center)
errp = MagnToFlux(m_Host-m_Host_err, nu_center) - flux
errn = flux - MagnToFlux(m_Host+m_Host_err, nu_center)
err_mean = (errp + errn)/2.
print(f"Host flux average (measured) in the {Filter} filter: {flux:2.3e} +- {err_mean:2.3e} [erg.cm-2.s-1]\n")

print(f"{Instrument} fluxes before host correction [erg.cm-2.s-1]:",F,"\n")

#----------Instrument properties---------#

if Instrument == "UVOT":
    #UVOT aperture ["]
    aperture = 5.0
    
    #Bandwidth (FWHM) V,B,U,UVW1,UVM2,UVW2 [A]
    lambd_FWHM = array([769, 975, 785, 693, 498, 657])
    #central wavelength V,B,U,UVW1,UVM2,UVW2 [A]
    lambd_center = array([5468, 4392, 3465, 2600, 2246, 1928])

    freq_min = convfreq(lambd_center + lambd_FWHM/2.)
    freq_max = convfreq(lambd_center - lambd_FWHM/2.)
    freq_center = convfreq(lambd_center)
    


if Instrument == "WISE":

    #ouverture wise ["]
    ouv = array([12.0, 6.5, 6.4, 6.1])
    #Pouv = ouv/reffR
    #ponderation1 = [8.51e-1, 7.14e-1, 7.14e-1, 7.03e-1] #tables de young #Luminosite dans ls champ de WISE/ luminosite totale de la galaxie (reffB)
    #ponderation1 = [7.14e-1, 5.51e-1, 5.27e-1, 5.27e-1] #tables de young #Luminosite dans ls champ de WISE/ luminosite totale de la galaxie (reffR)
    #ponderation1 = [6.215e-1, 4.55e-1, 4.4e-1, 4.3e-1] #extrapolation linéaire de reff (voir plot d'en bas pour la freq max de 2Mass) Reff = 7.73"
    
    
    freq_min1= []
    freq_max1= []
    freq_centr =[]
       
    #----------Filtes de WISE----------#
    
    #longueurs d'onde centrales des filtres de WISE [A]
    lambd_centr1=[2.2e5, 1.2e5, 4.6e4, 3.4e4]
    nu_centr1 = c/(array(lambd_centr1)*1e-10)
    reff_W= []
    ponderation1 = zeros(len(ouv))
    for i in range(len(lambd_centr1)):
        #reff_W.append(s(nu_centr1[i])) #extrapolation des rayons effectifs
        reff_W.append(reffR) #limite inférieure du rayon effectif
        ponderation1[i]= YOUNG(ouv[i],reff_W[i]) #fraction de luminosité ds le champ de wise (limite sup)
    
    #largeur a mi-hauteur des filtres W4, W3, W2, W1 [A]
    D_lambd1= [3.5e4, 8.8e4, 3.1e4, 3.5e4]
    
    
    for i in range (len(lambd_centr1)):
        freq_min1.append(convfreq(lambd_centr1[i] + D_lambd1[i]/2.))
        freq_max1.append(convfreq(lambd_centr1[i] - D_lambd1[i]/2.))
        freq_centr.append(convfreq(lambd_centr1[i]))
    
    #lecture des donnees PEGASE
    for i in range(len(freq_min1)):
        freq1=[]
        flux1=[]
        for j in range (len(matrix)):          
            if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= freq_min1[i]) and (float(matrix [j][0]) <= freq_max1[i]) :
                freq1.append(float(matrix [j][0]))
                flux1.append(float(matrix [j][2]))
        
        flux_pond1.append(moy(flux1)*ponderation1[i])
        flux_moy1.append(moy(flux1))
    
    print("flux de la galaxie hote dans les filtres W4, W3, W2, W1 de WISE [erg.cm-2.s-1]: \n",flux_pond1)
    print("flux moyen \n", flux_moy1)
    print("flux corrigé WISE \n", array(wise)- array(flux_pond1), "\n")
    #print "flux  WISE intermediare (entre valeur sans et avec corr)\n", array(wise)- array(flux_pond1)/2, "\n"

   

if Instrument == "TMASS":
    #ouverture 2MASS ["]
    ouv1 = 5.0
    Pouv1 = ouv1/reffRc
    #ponderation2 = 8.67e-1 #reffB 12.9sec
    #ponderation2 = 7.26e-1 #reffR 12.9sec
    #ponderation2 = 4.71e-1 #reffR 5 sec
    ponderation2 = 3.84e-1 #extrapolation linéaire de reff 5 sec (voir plot d'en bas pour la freq max de 2Mass) Reff = 7.73"
    
    freq_min2= []
    freq_max2= []
       
    #----------Filtes de 2MASS----------#
    
    #longueurs d'onde centrales des filtres de 2MASS [A]
    lambd_centr2=[2.16e4, 1.65e4, 1.25e4]
    
    #largeur a mi-hauteur des filtres Ks, H, J [A]
    D_lambd2= [2.2e3, 2.8e3, 3.0e3]
    
    
    
    for i in range (len(lambd_centr2)):
        freq_min2.append(convfreq(lambd_centr2[i] + D_lambd2[i]/2.))
        freq_max2.append(convfreq(lambd_centr2[i] - D_lambd2[i]/2.))


    for i in range(len(freq_min2)):
        freq2=[]
        flux2=[]
        for j in range (len(matrix)):          
            if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= freq_min2[i]) and (float(matrix [j][0]) <= freq_max2[i]) :
                freq2.append(float(matrix [j][0]))
                flux2.append(float(matrix [j][2]))
        
        flux_pond2.append(moy(flux2)*ponderation2)
        flux_moy2.append(moy(flux2))
        
    print("flux de la galaxie hote dans les filtres Ks, H, J de 2MASS[erg.cm-2.s-1]: \n",flux_pond2)
    print("flux moyen \n", flux_moy2)
    print("flux corrigé 2MASS \n", array(mass)- array(flux_pond2), "\n")
    


if Instrument == "ATOM":
    #ouverture ATOM ["]
    ouvA = 4.0
    #Pouv1 = ouv1/reffRc
    
#    freq_min2= []
#    freq_max2= []
       
    #----------Filtes de ATOM----------#
    
    #longueurs d'onde centrales du filtre Rc d'ATOM (these marcus hauser) [A]
    lambd_centr =  6250
    #nu_centr = c/(array(lambd_centr)*1e-10)
    nu_centr = convfreq(lambd_centr)
    ponderation= YOUNG(ouvA,Reff)
    
    
    #largeur a mi-hauteur du filtre R (these marcus hauser)[A]
    D_lambd= 1150
    
    
    freq_min = convfreq(lambd_centr + D_lambd/2.)
    freq_max = convfreq(lambd_centr - D_lambd/2.)
        
    #lecture des donnees PEGASE

    freq=[]
    flux=[]
    for j in range (len(matrix)):
        if (matrix [j][0][0])!= "#" and (float(matrix [j][0]) >= freq_min) and (float(matrix [j][0]) <= freq_max) :
            freq.append(float(matrix [j][0]))
            flux.append(float(matrix [j][2]))
            
    flux_pond = moy(flux)*ponderation
    flux_moy = moy(flux)
        
   
    print("flux moyen sur la bande R de ATOM sans prendre en compte l'ouverture\n", flux_moy,"\n")
    print("flux de la galaxie hote dans le filtre R de ATOM [erg.cm-2.s-1]: \n",flux_pond ,"\n")
    print("flux corrigé ATOM \n", atom - flux_pond, "\n")
        



host_frac = YOUNG_Sersic(aperture,Reff, N)
print(f"Fraction of the host galaxy luminosity through an aperture of {aperture:.1f} arsec:", host_frac , "\n")

    
mean_flux_instr = zeros(len(freq_min))
for i in range(len(freq_min)):
    mean_flux_instr[i] = Host_average(freq_min[i], freq_max[i], host_spectrum, Instrument) * host_frac
    
#Error propagation from the host flux error and instrumental error
A_errp = ones(len(F))*errp
F_corr_errp = sqrt(A_errp**2 + F_err**2)
A_errn = ones(len(F))*errn
F_corr_errn = sqrt(A_errn**2 + F_err**2)

print(f"Host galaxy flux in {Instrument} filters without aperture correction [erg.cm-2.s-1]: \n", mean_flux_instr/host_frac,"\n")
print(f"Host galaxy flux in {Instrument} filters within 5 arcsec aperture [erg.cm-2.s-1]: \n", mean_flux_instr ,"\n")
print(f"{Instrument} flux removed from the host galaxy emission \n", F - mean_flux_instr, "\n")
print(f"{Instrument} flux error after correction")
print(f"delta_F(-):", F_corr_errn)
print(f"delta_F(+):", F_corr_errp, "\n")
print(f"Percentage of host contamination {Instrument} \n",  (1-(F-mean_flux_instr)/F)*100)





plot(host_spectrum[0],host_spectrum[1],color="0.7",label = "Full host")
plot(host_spectrum[0],host_spectrum[1] * host_frac,color="0.1",label = f"Host within {Instrument} aperture")
errorbar(nu_center, flux, xerr = nu_center-nu_min, yerr=err_mean, fmt='o', label='Calibration flux')
errorbar(freq_center, F, xerr = freq_center-freq_min, yerr=F_err, fmt='o', label=f"{Instrument}")
errorbar(freq_center, F - mean_flux_instr, xerr = freq_center-freq_min, yerr=[F_corr_errn, F_corr_errp], fmt='o', label=f"{Instrument} corrected from host")
xlim([1e13,5e15])
ylim([2e-13,1.5e-11])
loglog()
xlabel(r"$\nu$ [Hz]", fontsize=14)
ylabel(r"$\nu F_\nu [\mathrm{erg}~\mathrm{cm}^{-2}~\mathrm{s}^{-1}]$", fontsize=14)
legend(loc="upper left")
tight_layout()




