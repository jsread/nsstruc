"""
.. module:: eos
   :synopsis: Define an eos object for TOV solving
.. moduleauthor:: Jocelyn Read <jread@fullerton.edu>

"""
import warnings
import numpy as np
from scipy import integrate
from scipy import interpolate
from scipy import sort
from scipy import optimize

refpressure = 10**14.7

class Eos:
    """ A class that defines cold equations of state for TOV and augmented TOV
    equations.
    
    Current parameters: either a list of functions, or a data file
    containing rho, 

    Based on eta = exp(enthalpy) - 1.0

    dens, press, energy are in units of g/cm^3
    (dens = number density in fm^-3 * m_b for m_b in conversions.py)

    enthalpy is dimensionless

    At minimum, defines:

        density: rest mass density as a function of pseudoenthalpy

        pressure: pressure as a function of pseudoenthalpy

        energy: energy density as a function of pseudoenthalpy

        dprden: derivative of pressure wrt energy as a function of pseudoenthalpy,
        also the speed of sound squared

        enthalpy: pseudoenthalpy as a function of pressure

        range: a tuple of (enthalpy_low, eta_hi)

    From table filename, also defines:

        densitypoints: read from table

        energypoints: read from table

        pressurepoints: read from table

        enthalpypoints: numerically integrated pseudoenthalpy
    
        """
    def __init__(self, funcs=None, filename=None):
        if funcs:
            # just define the functions explicitly
            if filename: warnings.warn("given funcs; ignoring filename")
            (density, pressure, energy, dprden, enthalpy) = funcs
            self.density = density
            self.pressure = pressure
            self.energy = energy
            self.dprden= dprden
            self.enthalpy = enthalpy 
        elif filename:
            # load table
            table = np.loadtxt(filename)
            table = sort(table,0) # consistent lowest-to-highest-order

            #enforce monotonic increase in pressure
            while not all(np.diff(table[:,1]) >= 0): # not monotonic increase
                warnings.warn("Removing non-monotonic pressure points")
                # find difference in pressure with next value 
                dp = np.hstack((np.diff(table[:,1]),[0]))
                # delete the lower-density point in a non-monotonic pair
                table = [entry for (i,entry) in enumerate(table) if dp[i] <= 0]

            self.densitypoints = table[:,0]
            self.pressurepoints = table[:,1]
            self.energypoints = table[:,2]
            
            # define a thermodynamically consistent pseudo-enthalpy
            # using the given pressure and energy values
            # and the derivative defined in the thermodynamics notes
            denthalpy = 1.0 / (self.pressurepoints + self.energypoints) 
            enthalpy = integrate.cumtrapz(denthalpy,self.pressurepoints)
            # we are assuming the lowest density point in the table is the
            # surface where pseudoenthalpy vanishes
            self.enthalpypoints = np.hstack((0.0, enthalpy))
            self.range = (self.enthalpypoints[0],self.enthalpypoints[-1])


            smoothing = 0
            energyrep = interpolate.splrep(self.enthalpypoints, 
                                            np.log10(self.energypoints), 
                                            k=3, s=smoothing)
            pressrep = interpolate.splrep(self.enthalpypoints, 
                                            np.log10(self.pressurepoints), 
                                            k=3, s=smoothing)

            def pressure(enthalpy):
                '''pressure over c^2 in g/cm^3 as a function of pseudoenthalpy
                calculated from tabled energy and pressure and numerical
                integration of pseudoenthalpy'''
                return 10**(interpolate.splev(enthalpy, pressrep, der=0))
            self.pressure = np.vectorize(pressure)
            def energy(enthalpy):
                '''energy density in g/cm^3 as a function of pseudoenthalpy
                calculated from tabled energy and pressure and numerical
                integration of pseudoenthalpy'''
                return 10**(interpolate.splev(enthalpy, energyrep, der=0))
            self.energy = np.vectorize(energy)

            def dprden(enthalpy):
                '''dimensionless derivative of pressure with respect to energy 
                as a function of pseudoenthalpy; the speed of sound squared''' 
                # use the derivatives of the interpolated energy and pressure
                # interpolations are of log energy and log pressure
                # so need to also calculate energy and pressure to convert
                dlogen = interpolate.splev(enthalpy, energyrep, der=1)
                dlogp = interpolate.splev(enthalpy, pressrep, der=1)
                en = 10**interpolate.splev(enthalpy, energyrep, der=0)
                pr = 10**interpolate.splev(enthalpy, pressrep, der=0)
                # these are not natural logs, but the base factors cancel 
                return pr/en*dlogp/dlogen
            self.dprden= np.vectorize(dprden)

            # Using thermodynamically consistent density. Varies from the tabled
            # densities by ~1% due to the approximate integration of enthalpy
            # (or inconsistent tables)
            # TODO: option for using interpolated table density instead?
            def density(enthalpy):
                '''rest mass density in g/cm^3 as a function of pseudoenthalpy
                calculated from tabled energy and pressure and numerical
                integration of pseudoenthalpy'''
                return ( 10**(interpolate.splev(enthalpy, pressrep, der=0)) 
                        + 10**(interpolate.splev(enthalpy, energyrep, der=0))
                        ) / np.exp(enthalpy)
            self.density = np.vectorize(density)
            
            enthrep = interpolate.splrep(np.log10(self.energypoints),
                                         self.enthalpypoints,
                                         k=3, s=smoothing)
            def enthalpy(energy):
                '''pseudoenthalpy as a function of energy density in g/cm^3
                calculated from tabled energy and pressure and numerical
                integration of pseudoenthalpy'''
                return interpolate.splev(np.log10(energy), enthrep)
            self.enthalpy = np.vectorize(enthalpy)


def createpolytrope(p2, g):
    ''' Create function set for EOS object using a reference pressure, specified
        by  log10(p/c^2)  at rest mass density 10**14.7 and gamma, store in
        generic piecewise polytrope format''' 
    g = float(g)
    z = p2 - g * 14.7 # Determine log K 
    K = 10**z
    n = 1 / (g - 1)
    print( K, g)
    def density(enthalpy):
        '''rest mass density in g/cm^3 as a function of pseudoenthalpy'''
        eta = np.exp(enthalpy) - 1
        return (eta * (g - 1) / (K * g) )**(n)
    def pressure (enthalpy):
        '''energy density in g/cm^3 as a function of pseudoenthalpy'''
        eta = np.exp(enthalpy) - 1
        return K * ( eta * (g - 1) / ( K * g ) )**(n * g)
    def energy(enthalpy):
        '''energy density in g/cm^3 as a function of pseudoenthalpy'''
        # convert from pseudoenthalpy into parameter8 eta = specific enthalpy - 1
        eta = np.exp(enthalpy) - 1
        if (eta != 0):
            return (g + eta) / g * ( K * g * n / eta )**(-n)
        else:
            return 0
    def dprden(enthalpy):
        '''dimensionless derivative of pressure with respect to energy 
        as a function of pseudoenthalpy; the speed of sound squared''' 
        eta = np.exp(enthalpy) - 1
        if (eta != 0):
            return n * (1 + eta) / eta
        else:
            return 0
    def enthalpy(energy):
        '''pseudoenthalpy as a function of energy density in g/cm^3'''
        # No nice analytic form
        def diffen(eta):
            if (eta != 0):
                return energy - (g + eta) / g * ( K * g * n / eta )**(-n) 
            else:
                return 0
        return optimize.zeros.brentq(diffen, 10**(-10), 5)
    return (np.vectorize(density), np.vectorize(pressure),
            np.vectorize(energy),  np.vectorize(dprden), np.vectorize(enthalpy))
            

def expandedparams(params):
  """ for piecewise polytropes
      given a parameter set of [[log rhodivs][log K][gamma]]
      return a parameter set of [[ens][etas][K][gamma][a][n]]
  """
  Kis = 10**array(params[1])
  rhois = 10**array(params[0])
  gis = array(params[2])
  ais = [0.0]
  for i in range(len(rhois)-1):
    ai = ais[i] + Kis[i] / (gis[i] - 1) * rhois[i+1]**(gis[i] - 1)  \
          - Kis[i+1] / (gis[i+1] - 1) * rhois[i+1]**(gis[i+1] - 1)
    ais.append(ai)
  ais = array(ais)
  etas = ais + Kis*gis/(gis-1)*rhois**(gis-1)
  nis = 1/(gis - 1)
  enis = ( (etas - ais) / (Kis * (nis + 1) ))**nis \
       * ( 1 + (ais + nis* etas) / ( nis + 1) )
  return vstack((enis, etas, Kis, gis, ais, nis))
