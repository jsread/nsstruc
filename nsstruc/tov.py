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

def eosderiv(eos, rm_vector, eta):
    (rad, mass) = rm_vector
    drdeta =  - rad * (rad - 2 * mass) \
              / (mass + 4 * pi * rad **3  * eos.pressfunc(eta)*Gcc) \
              / (eta+1) 
    dmdeta = 4 * pi * eos.energyfunc(eta)*Gcc * rad **2 * drdeta
    return np.array([drdeta,dmdeta])

def rotderiv(eos, rm_vector, eta):
  (rad, mass, restmass, alpha, omega) = rm_vector

  pr = eos.pressfunc(eta) * Gcc
  en = eos.energyfunc(eta) * Gcc
  rho = (pr + en)/(eta+1)

  factor = rad - 2 * mass
  coef = 4.0*pi

  drdeta =  - rad * factor / (mass + coef * rad**3 * pr) / (eta+1)
  dmdeta = coef * en * rad**2 * drdeta
  dmadeta = coef * rho * rad**2 * drdeta / ( factor / rad )**0.5
  domega = alpha * drdeta
  dalpha = drdeta * (- 4.0 * alpha / rad \
      + coef * rad * (pr + en) * (4.0 * omega + rad * alpha ) / factor )

  return array([drdeta,dmdeta, dmadeta, dalpha, domega])

def tidederiv(eos,rm_vector, eta):
    (rad, mass, restmass, beta, h) = rm_vector

    pr = eos.pressfunc(eta) * Gcc
    en = eos.energyfunc(eta) * Gcc
    func = eos.dendp(eta)
    rho = (pr + en)/(eta+1)

    factor = rad - 2 * mass
    coef = 4.0*pi

    drdeta =  - rad * factor / (mass + coef * rad**3 * pr) / (eta+1)
    dmdeta = coef * en * rad**2 * drdeta
    dmadeta = coef * rho * rad**2 * drdeta / ( factor / rad )**0.5
    dh = beta * drdeta
    dbeta = 2 * drdeta *\
        ( h * ( \
          (- 2. * pi * func * rad * (pr + en))/ factor + \
          (2. * mass**2 + \
            rad**2 * (3.0 + 2.0 * pi * rad**2 * \
              ( pr * (16.0 * pi * pr * rad**2 - 9.0 ) - 5. * en ))\
            + mass * ( - 6 * rad + 4 * pi * rad**3 * (13 *pr + 5*en) ) )\
          / rad**2 / factor**2 )
        + beta / rad / factor * ( mass - rad + 2 * pi * rad**3 * (en -
          pr)))
    return array([drdeta,dmdeta, dmadeta, dbeta, dh])


def profile(eos, en_c, number):
  # set up a core with small initial radius and mass
  initm = 4.0/3.0*pi*(initr)**3 * en_c * Gcc
  init_rm = array([initr,initm])
  eta_c = eos.etafunc(en_c)
  eta_0 = eta_surf
  # resolve core with linear energy spacing
  linetas = linspace(eta_c, eta_c/10, number+1)[:-1]
  # resolve crust with log eta spacing
  logetas = logspace(log10(linetas[-1]),log10(eta_0), number+1)
  etas = hstack((array(linetas),array(logetas[1:])))
  # define derivative for this eos
  def deriv(rm_vector,eta):
    return eosderiv(eos, rm_vector, eta)
  rms = integrate.odeint(deriv, init_rm, etas)
  emrar = vstack((eos.energyfunc(etas),transpose(rms)))
  return emrar

def rotprofile(eos, en_c, number):
  # set up a core with small initial radius and mass
  initm = 4.0/3.0*pi*(initr)**3 * en_c * Gcc
  init_rm = array([initr,initm,initm, 0.0, 1.0])
  eta_c = eos.etafunc(en_c)
  eta_0 = eta_surf
  # resolve core with linear energy spacing
  linetas = linspace(eta_c, eta_c/10, number+1)[:-1]
  # resolve crust with log eta spacing
  logetas = logspace(log10(linetas[-1]),log10(eta_0), number+1)
  etas = hstack((array(linetas),array(logetas[1:])))
  # define derivative for this eos
  def deriv(rm_vector,eta):
    return rotderiv(eos, rm_vector, eta)
  rms = integrate.odeint(deriv, init_rm, etas)
  emrar = vstack((eos.energyfunc(etas),transpose(rms)))
  return emrar

def tideprofile(eos, en_c, number):
  # set up a core with small initial radius and mass
  initm = 4.0/3.0*pi*(initr)**3 * en_c * Gcc
  init_rm = array([initr,initm,initm, 0.1* 2 * initr,0.1 *initr**2])
  eta_c = eos.etafunc(en_c)
  eta_0 = eta_surf
  # resolve core with linear energy spacing
  linetas = linspace(eta_c, eta_c/10, number+1)[:-1]
  # resolve crust with log eta spacing
  logetas = logspace(log10(linetas[-1]),log10(eta_0), number+1)
  etas = hstack((array(linetas),array(logetas[1:])))
  # define derivative for this eos
  def deriv(rm_vector,eta):
    return tidederiv(eos, rm_vector, eta)
  rms = integrate.odeint(deriv, init_rm, etas)
  emrar = vstack((eos.energyfunc(etas),transpose(rms)))
  return emrar

def moment_of_inertia( rotprof_surface):
  (en, rad, mass, restmass, alpha, omega) = rotprof_surface
  conv = 1.34685326e43 # g cm^2 per km^3
  I = 1.0/6.0 * rad**4 * alpha / (omega + 1.0/3.0 * rad * alpha)
  return I * conv


def compactness(tidalprof_surface):
  (en, rad, mass, restmass, beta, hfunc) = tidalprof_surface
  return mass / rad


def k2(tidalprof_surface):
  (en, rad, mass, restmass, beta, hfunc) = tidalprof_surface
  y = rad * beta / hfunc
  C = mass / rad
  return 8. / 5. * C**5 * (1. - 2. * C)**2 * (2. + 2. * C * (y - 1.) - y) \
      / (2. * C * (6. - 3. * y + 3. * C * (5. * y - 8.)) \
      + 4. * C**3. * (13. - 11.*y + C * (3.*y - 2.) + 2. * C**2 * (1. + y)) \
      + 3. * ( 1. - 2.*C)**2 * ( 2. - y + 2.*C * (y - 1.)) * log(1. - 2.*C))

def h2(tidalprof_surface):
  (en, rad, mass, restmass, beta, hfunc) = tidalprof_surface
  y = rad * beta / hfunc
  C = mass / rad
  Cfac = 1. - 2. * C
  return C**5 * ( 8. * Cfac * C**2 * y + 16 * (1. - C)**3 ) \
      * (8. * C**5 * (y + 1.) + 4. * C**4 * ( 3. * y - 2.) \
        + 4. * C**3 * (13. - 11. * y ) + 6. * C**2 * (5. * y - 8) \
        + 6. * C * (2. - y) \
        + 3 * Cfac**2 * (2. * C * (y-1.) - y + 2) * log(Cfac) \
        )**(-1)

def energy_of_mass(eos, mass, lowenergy, highenergy):
  def massdiff(x, eos, mass):
    return profile(eos, x, 100)[2,-1] - mass*Solarmass_km
  return optimize.zeros.brenth(massdiff, lowenergy, highenergy, \
      (eos,mass))

def energy_maxmass(eos, lowenergy, highenergy):
  massen = lambda x: 10.0 - profile(eos, x, 100)[2,-1]
  maxen = optimize.brent(massen, brack=(lowenergy,highenergy),
        tol=1e-4)
  maxpr = profile(eos,maxen, 100)[:,-1]
  return maxen, maxpr[2]/Solarmass_km, maxpr[1]

def properties_of_energy(eos,en):
    prof = tideprofile(eos, en, 100)
    (en, rad, mass, restmass, beta, hfunc) = prof[:,-1]
    k2val = k2(prof[:,-2])
    tidelambda = 2. / 3. * k2val * (rad / mass)**5.
    return (mass/Solarmass_km, rad, tidelambda, k2val)

def mrarray(eos, lowenergy, highenergy,number):
  eps = logspace(log10(lowenergy),log10(highenergy),number)
  table = []
  for ep_c in eps:
    prof = profile(eos, ep_c, 100)
    (en, rad, mass) = prof[:,-1]
    table.append((ep_c, mass/Solarmass_km,rad))
  return array(table)

def lambdamarray(eos, lowenergy, highenergy,number):
  eps = logspace(log10(lowenergy),log10(highenergy),number)
  table = []
  for ep_c in eps:
    prof = tideprofile(eos, ep_c, 100)
    (en, rad, mass, restmass, beta, hfunc) = prof[:,-1]
    k2val = k2(prof[:,-2])
    tidelambda = 2. / 3. * k2val * (rad / mass)**5.
    table.append((ep_c, mass/Solarmass_km, rad, tidelambda, k2val))
  return array(table)

