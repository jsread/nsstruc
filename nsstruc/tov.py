import numpy as np
import eos
from scipy import integrate
from scipy import optimize
from scipy import interpolate

from conversions import *

# G/c^2 in m / kg 
# gram per kg * cm per m3 is m^-2 per g/cm^3
# 1000**2 to convert to km^-2 per g/cm^3
Gcc = G/C**2  / Gram_per_kg * Cm3_per_m3  * 1000**2
            ## Check in Google: G/C^2 * (1 g/cm^3) in (km^-2)
            ## about 7.4247 x 10^-19

def deriv(rm_vector, enthalpy, eos):
    ''' The TOV equations in terms of pseudoenthalpy. Derivatives of vector
    (radius, mass) as a function of enthalpy for given eos.'''
    (rad, mass) = rm_vector
    drdenth =  - rad * (rad - 2 * mass) \
              / (mass + 4 * np.pi * rad **3  * eos.pressure(enthalpy)*Gcc) 
    dmdenth = 4 * np.pi * eos.energy(enthalpy)*Gcc * rad **2 * drdenth
    return np.array([drdenth,dmdenth])

def profile(eos, energy_center , number=200, initr=0.001, energy_0=1e7,
        rad_0=False):
    '''Calculate the enthalpy profile of a neutron star with equation of
    state eos and central energy density energy_center in g/cm^3

    optional arguments:
        
        number: number of points

        initr: size of central seed

        energy_0: energy density at which to stop resolving crust. This doesn't
            need to be tiny (and may cause problems if it is). Integration still
            goes to surface.

        rad_0: Set true to calculate outer radius at energy_0 instead of at
        formal zero enthalpy limit of EOS
    
    '''

    # set up a core with small initial radius and mass
    initm = 4.0/3.0*np.pi*(initr)**3 * energy_center * Gcc
    init_rm = np.array([initr,initm])

    # central enthalpy
    enth_c = eos.enthalpy(energy_center)

    # end the core spacing here 
    enth_0 = eos.enthalpy(energy_center/10) 

    # resolve core with linear spacing in enthalpy
    linenths = np.linspace(enth_c, enth_0, number/2)[:-1]

    # end the crust spacing here and jump to 0
    enth_0 = eos.enthalpy(energy_0) 
    # resolve crust with log spacing in enthalpy
    logenths = np.logspace(np.log10(linenths[-1]),np.log10(enth_0),
                           number/2+1)[1:]

    # set of enthalpies to calculate energy density, radius, and mass at
    if rad_0:
        enths = np.hstack((np.array(linenths),np.array(logenths)))
    else:
        enths = np.hstack((np.array(linenths),np.array(logenths),np.array(0)))
    # enths = np.linspace(enth_c, 0, number) # boring way for debugging

    # integrate the TOV equations specified in deriv above
    rms = integrate.odeint(deriv, init_rm, enths, args=(eos,))

    return np.vstack((enths,np.transpose(rms)))


def rotderiv(rm_vector, enthalpy, eos):
    ''' The TOV equations augmented for a linearly perturbed, slowly
    rotating neutron star.  Also includes rest mass of star.'''
    (rad, mass, restmass, alpha, omega) = rm_vector
  
    pr = eos.pressure(enthalpy) * Gcc
    en = eos.energy(enthalpy) * Gcc
    rho = eos.density(enthalpy) ** Gcc
  
    factor = rad - 2 * mass
    coef = 4.0 * np.pi

    drdenth =  - rad * factor / (mass + coef * rad**3 * pr) 
    dmdenth = coef * en * rad**2 * drdenth
    dmadenth = coef * rho * rad**2 * drdenth / ( factor / rad )**0.5
    domega = alpha * drdenth
    dalpha = drdenth * (- 4.0 * alpha / rad \
      + coef * rad * (pr + en) * (4.0 * omega + rad * alpha ) / factor )

    return np.array([drdenth,dmdenth, dmadenth, dalpha, domega])


def rotprofile(eos, energy_center , number=100, initr=0.001, energy_0=1e7):
    '''Calculate the enthalpy profile of a slowly rotating neutron star with
    equation of state eos and central energy density energy_center in g/cm^3.

    optional arguments:
        
        number: number of points

        initr: size of central seed

        energy_0: energy density at which to stop resolving crust. This doesn't
            need to be tiny (and may cause problems if it is). Integration still
            goes to surface.
    
    '''
    #TODO: consolodate the profile functions; this mostly just repeats
    #profile but with different initial conditions

    # central enthalpy and density
    enth_c = eos.enthalpy(energy_center)
    dens_c = eos.density(enth_c)

 
    # set up a core with small initial radius and mass
    initm = 4.0/3.0*np.pi*(initr)**3 * energy_center * Gcc
    initma = 4.0/3.0*np.pi*(initr)**3 * dens_c * Gcc

    init_rm = np.array([initr,initm, initma, 0.0, 1.0])

    # resolve core with linear spacing
    linenths = np.linspace(enth_c, enth_c/10, number)[:-1]

    # resolve crust with log enth spacing
    enth_0 = eos.enthalpy(energy_0) # end the crust spacing here and jump to 0
    logenths = np.logspace(np.log10(linenths[-1]),np.log10(enth_0),
                           number/2+1)

    # set of enthalpies to calculate energy density, radius, and mass at
    enths = np.hstack((np.array(linenths),np.array(logenths[1:]),np.array(0)))
    # enths = np.linspace(enth_c, 0, number) # boring way for debugging

    # integrate the TOV equations specified in deriv above
    rms = integrate.odeint(rotderiv, init_rm, enths, args=(eos,))

    return np.vstack((enths,np.transpose(rms)))

def tidederiv(rm_vector, enthalpy, eos):
    (rad, mass, restmass, beta, h) = rm_vector

    pr = eos.pressure(enthalpy) * Gcc
    en = eos.energy(enthalpy) * Gcc
    func = 1.0/eos.dprden(enthalpy)

    factor = rad - 2 * mass
    coef = 4.0*np.pi

    drdenth =  - rad * factor / (mass + coef * rad**3 * pr)
    dmdenth = coef * en * rad**2 * drdenth
    dbeta = 2 * drdenth *\
        ( h * ( \
          (- 2. * np.pi * func * rad * (pr + en))/ factor + \
          (2. * mass**2 + \
            rad**2 * (3.0 + 2.0 * np.pi * rad**2 * \
              ( pr * (16.0 * np.pi * pr * rad**2 - 9.0 ) - 5. * en ))\
            + mass * ( - 6 * rad + 4 * np.pi * rad**3 * (13 *pr + 5*en) ) )\
          / rad**2 / factor**2 )
        + beta / rad / factor * ( mass - rad + 2 * np.pi * rad**3 * (en -
          pr)))
    dh = beta * drdenth
    return np.array([drdenth,dmdenth, dbeta, dh])


def tideprofile(eos, energy_center , number=200, initr=0.001, energy_0=1e7):
    '''Calculate the enthalpy profile of a tidally perturbed neutron star
    with equation of state eos and central energy density energy_center in g/cm^3.

    optional arguments:
        
        number: number of points

        initr: size of central seed

        energy_0: energy density at which to stop resolving crust. This doesn't
            need to be tiny (and may cause problems if it is). Integration still
            goes to surface.
    
    '''
    #TODO: consolodate the profile functions; this mostly just repeats
    #profile but with different initial conditions

    # central enthalpy and density
    enth_c = eos.enthalpy(energy_center)
    dens_c = eos.density(enth_c)


    # set up a core with small initial radius and mass
    initm = 4.0/3.0*np.pi*(initr)**3 * energy_center * Gcc
    # seed perturbation amplitude of 0.1
    init_rm = np.array([initr,initm,initm, 0.1* 2 * initr,0.1 *initr**2])

    # resolve core with linear spacing
    linenths = np.linspace(enth_c, enth_c/10, number/2)[:-1]

    # resolve crust with log enth spacing
    enth_0 = eos.enthalpy(energy_0) # end the crust spacing here and jump to 0
    logenths = np.logspace(np.log10(linenths[-1]),np.log10(enth_0),
                           number/2 + 1)

    # set of enthalpies to calculate energy density, radius, and mass at
    enths = np.hstack((np.array(linenths),np.array(logenths[1:]),np.array(0)))
    # enths = np.linspace(enth_c, 0, number) # boring way for debugging

    # integrate the TOV equations specified in deriv above
    rms = integrate.odeint(tidederiv, init_rm, enths, args=(eos,))

    return np.vstack((enths,np.transpose(rms)))



def energy_of_mass(eos, mass, lowenergy, highenergy=None):
    ''' for a given eos and a stellar mass in solar masses, calculate the
    central energy density in g/cm^3. Include a lower bound on the energy
    to avoid convergence problems especially with crusts.  About 10**14
    g/cm^3 works well for most neutron stars.

    Optional params: 

      highenergy = upper bound for search.  If this is not specified, use
      the limiting enthalpy range of the EOS. 
      '''
    if highenergy == None:
        highenergy = eos.energy(eos.enthrange[1])
    def massdiff(x, eos, mass):
        return profile(eos, x)[2,-1] - mass*Solarmass_km
    try:
        return optimize.brenth(massdiff, lowenergy, highenergy, (eos,mass))
    except RuntimeError:
        warnings.warn("Problem getting convergence in optimze.brenth")
        return 0


def energy_maxmass(eos, lowenergy, highenergy=None):
    ''' Determine properties of the maximum mass star supported by the EOS.
    Returns (central energy density in g/cm^3, mass in Solar masses, and
    radius in km  for the star). 
    
    Include a lower bound on the energy to
    avoid convergence problems especially with crusts.  About 10**14 g/cm^3
    works well for most neutron stars.

    Optional params: 

      highenergy = upper bound for search. If this is not specified, use
      the limiting enthalpy range of the EOS. 
      '''
    if highenergy == None:
        highenergy = eos.energy(eos.enthrange[1])
    # Calculate the maximum mass energy density
    massen = lambda x: 10.0 - profile(eos, x)[2,-1]
    maxen = optimize.brent(massen, brack=(lowenergy,highenergy))
    # Calculate the maximum mass model
    maxprofile = profile(eos,maxen)[:,-1]
    return maxen, maxprofile[2]/Solarmass_km, maxprofile[1]

def mrarray(eos, lowenergy, highenergy, number=20):
    eps = np.logspace(np.log10(lowenergy),np.log10(highenergy),number)
    table = []
    for ep_c in eps:
        prof = profile(eos, ep_c)
        (en, rad, mass) = prof[:,-1]
        table.append((ep_c, mass/Solarmass_km,rad))
    return np.array(table)


def moment_of_inertia( rotprof_surface):
    (en, rad, mass, restmass, alpha, omega) = rotprof_surface
    conv = 1.34685326e43 # g cm^2 per km^3
    I = 1.0/6.0 * rad**4 * alpha / (omega + 1.0/3.0 * rad * alpha)
    return I * conv

def k2(tidalprof_surface):
    (en, rad, mass, restmass, beta, hfunc) = tidalprof_surface
    y = rad * beta / hfunc
    C = mass / rad
    return 8. / 5. * C**5 * (1. - 2. * C)**2 * (2. + 2. * C * (y - 1.) - y) \
      / (2. * C * (6. - 3. * y + 3. * C * (5. * y - 8.)) \
      + 4. * C**3. * (13. - 11.*y + C * (3.*y - 2.) + 2. * C**2 * (1. + y)) \
      + 3. * ( 1. - 2.*C)**2 * ( 2. - y + 2.*C * (y - 1.)) * np.log(1. - 2.*C))

def properties_of_energy(eos,en):
    ''' given a central energy density in g/cm^3, output some relavant
    properties for a tidally deformed star: mass in Solar masses, radius in
    km, dimensionless Lambda parameter, and k2.'''
    prof = tideprofile(eos, en)
    (en, rad, mass, restmass, beta, hfunc) = prof[:,-1]
    k2val = k2(prof[:,-2])
    tidelambda = 2. / 3. * k2val * (rad / mass)**5.
    return (mass/Solarmass_km, rad, tidelambda, k2val)

def lambdamarray(eos, lowenergy, highenergy, number):
    eps = np.logspace(np.log10(lowenergy),np.log10(highenergy),number)
    table = []
    for ep_c in eps:
        prof = tideprofile(eos, ep_c)
        (en, rad, mass, restmass, beta, hfunc) = prof[:,-1]
        k2val = k2(prof[:,-2])
        tidelambda = 2. / 3. * k2val * (rad / mass)**5.
        table.append((ep_c, mass/Solarmass_km, rad, tidelambda, k2val))
    return np.array(table)

