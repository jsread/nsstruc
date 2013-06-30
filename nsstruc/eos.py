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
    
    Current parameters: either a tuple/list of functions, or a data file
    containing rho (rest mass density), pressure density / c^2, and energy
    density all in g/cm^3.

    density, pressure, and energy are returned in units of g/cm^3 
    (dens = number density in fm^-3 * m_b for m_b in conversions.py)

    enthalpy is dimensionless

    At minimum, defines:

        density: rest mass density as a function of pseudoenthalpy

        pressure: pressure as a function of pseudoenthalpy

        energy: energy density as a function of pseudoenthalpy

        dprden: derivative of pressure wrt energy as a function of pseudoenthalpy,
        also the speed of sound squared

        enthalpy: pseudoenthalpy as a function of pressure

        enthrange: a tuple of (enthalpy_low, enthalpy_high)

    From table filename, also defines:

        densitypoints: read from table

        energypoints: read from table

        pressurepoints: read from table

        enthalpypoints: numerically integrated pseudoenthalpy
    
        """
    def __init__(self, funcs=None, enthrange=None, filename=None):
        if funcs:
            # just define the functions explicitly
            if filename: warnings.warn("given funcs; ignoring filename")
            (density, pressure, energy, dprden, enthalpy) = funcs
            self.density = density
            self.pressure = pressure
            self.energy = energy
            self.dprden= dprden
            self.enthalpy = enthalpy 
            if enthrange:
                self.enthrange = enthrange
            else:
                self.enthrange = (1e-40,3.0)
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
            self.enthrange = (self.enthalpypoints[0],self.enthalpypoints[-1])


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


def polytropefuncs(p2, g, enthalpybounds=None):
    ''' Create function set for EOS object using gamma and a reference pressure
    p2, specified as log10(p/c^3 in g/cm^3), at rest mass density 10**14.7 also
    in g/cm^3.

    To invert enthalpy as a function of energy, bounds are required. Default
    is (1e-40,2). Can specify tuple as enthalpybounds=(min,max) if the energies
    you desire are outside this range.  ''' 
    g = float(g)
    z = p2 - g * 14.7 # Determine log K 
    K = 10**z
    n = 1 / (g - 1)

    if enthalpybounds:
        (lowenth, highenth) = enthalpybounds
    else:
        (lowenth,highenth) = (1e-40,3.0)

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
        try:
            return optimize.zeros.brentq(diffen, lowenth, highenth)
        except ValueError:
            return 0
    return (np.vectorize(density), np.vectorize(pressure),
            np.vectorize(energy),  np.vectorize(dprden), np.vectorize(enthalpy))
            
