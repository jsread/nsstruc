# 2006 CODATA recommended values
import pylab

M_b_g = 1.66e-24 # Average mass of baryon in grams

Fm3_per_cm3= 1.0e39 

M_b = M_b_g * Fm3_per_cm3 # Convert from number density to g cm^-3
	
G = 6.6728e-11 # in m^3 kg^-1 s^-2
C = 299792458 # in m s^-1
KM_per_kg = G/C**2 / 1000.0 # Convert kg to km for geometrized units
Gram_per_KM = 1.34659496e33 # Convert km to g for geometrized units
Erg_per_Gram = 8.98744179e20 # convert g to erg for geometrized units
LOGepg = pylab.log10(Erg_per_Gram) 

Gram_per_MeV = 1.782661758e-27 #convert from MeV to gram
#particle physics convention to GR convention for energy

Solar_mass_kg = 1.9891e30 #traditional convention
Solarmass_kg = 1.98892e30 #modern value
Solarmass_km = Solarmass_kg * G / C**2 /1000.0
	
Gram_per_kg = 1.0e3
Cm3_per_km3 = 1.0e15
Cm3_per_m3 = 1.0e6

C_cgs = 29979245800.0
