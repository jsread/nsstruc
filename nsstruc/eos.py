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
                self.enthrange = (1e-40,2.0)
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
            # enthalpy range where sensibly defined by table
            self.enthrange = (self.enthalpypoints[1],self.enthalpypoints[-1])


            logpressure = interpolate.InterpolatedUnivariateSpline(
                                            self.enthalpypoints,
                                            np.log10(self.pressurepoints),k=1)
            
            def pressure(enthalpy):
                if enthalpy < self.enthrange[-1]:
                    if enthalpy > self.enthrange[0]:
                        warnings.warn(
                                "Below bound in interpolation of tabled values")
                        return 10**logpressure(enthalpy).astype(np.float32)
                    else:
                        return 0.0
                else:
                    warnings.warn("Above bound in interpolation of tabled values")
                    return self.pressurepoints[-1]
            self.pressure = np.vectorize(pressure, 
                        doc='''pressure over c^2 in g/cm^3 as a
                      function of pseudoenthalpy; calculated from tabled
                      energy and pressure and numerical integration of
                      pseudoenthalpy''')
            
            logenergy = interpolate.InterpolatedUnivariateSpline(
                                            self.enthalpypoints,
                                            np.log10(self.energypoints),k=1)

            def energy(enthalpy):
                if enthalpy < self.enthrange[-1]:
                    if enthalpy > self.enthrange[0]:
                        warnings.warn(
                                "Below bound in interpolation of tabled values")
                        return 10**logenergy(enthalpy).astype(np.float32)
                    else:
                        return 0.0
                else:
                    warnings.warn("Above bound in interpolation of tabled values")
                    return self.energypoints[-1]
            self.energy=np.vectorize(energy,
                            doc='''energy density in g/cm^3 as a function
                      of pseudoenthalpy; calculated from tabled energy and
                      pressure and numerical integration of
                      pseudoenthalpy''')

            dlogen = logenergy.derivatives(enthalpy[0])[1]
            dlogp = logpressure.derivatives(enthalpy[0])[1]
            en = float(10**logenergy(enthalpy[0])) 
            pr = float(10**logpressure(enthalpy[0]))
            dpdr0 = pr / en * dlogp/dlogen
 
            def dprden(enthalpy):
                # use the derivatives of the interpolated energy and pressure
                # interpolations are of log energy and log pressure
                # so need to also calculate energy and pressure to convert
                if enthalpy < self.enthrange[-1]:
                    if enthalpy > self.enthrange[0]:
                        warnings.warn(
                                "Below bound in interpolation of tabled values")
                        dlogen = logenergy.derivatives(enthalpy)[1]
                        dlogp = logpressure.derivatives(enthalpy)[1]
                        en = float(10**logenergy(enthalpy)) 
                        pr = float(10**logpressure(enthalpy))
                        # these are not natural logs, but the base factors cancel 
                        return pr/en*dlogp/dlogen
                    else:
                        return dpdr0
                else:
                    warnings.warn("Above bound in interpolation of tabled values")
                    return 1.0
            self.dprden = np.vectorize(dprden,
                        doc='''dimensionless derivative of pressure
                      with respect to energy, as a function of
                      pseudoenthalpy; the speed of sound squared''')

            # Using thermodynamically consistent density. Varies from the tabled
            # densities by ~1% due to the approximate integration of enthalpy
            # (or inconsistent tables)
            # TODO: option for using interpolated table density instead?
            def density(enthalpy):
                if enthalpy < self.enthrange[-1]:
                    if enthalpy > self.enthrange[0]:
                        warnings.warn(
                                "Below bound in interpolation of tabled values")
                        return ( 10**logenergy(enthalpy) +
                         10**logpressure(enthalpy) ) / np.exp(enthalpy)
                    else:
                        return 0.0
                else:
                    warnings.warn("Above bound in interpolation of tabled values")
                    return self.densitypoints[-1]
            self.density = np.vectorize(density,
                            doc= '''rest mass density in g/cm^3 as a
                      function of pseudoenthalpy calculated from tabled
                      energy and pressure and numerical integration of
                      pseudoenthalpy''')
 
            
            enthalpy_of_log = interpolate.InterpolatedUnivariateSpline(
                                          np.log10(self.energypoints),
                                          self.enthalpypoints, k=1)
            def enthalpy(energy):
                return enthalpy_of_log(np.log10(energy))
            self.enthalpy = enthalpy
#            self.enthalpy.__doc__ = '''pseudoenthalpy as a function of
#                      energy density in g/dcm^3 calculated from tabled
#                      energy and pressure and numerical integration of
#                      pseudoenthalpy'''


def polytropefuncs(p_ref, g, enthalpybounds=None):
    ''' Create function set for EOS object using gamma and a reference pressure
    p_ref, specified as log10(p/c^3 in g/cm^3), at rest mass density 10**14.7 also
    in g/cm^3. Typical values are ~ 13-14.

    p = K rho ^ g ; rho is rest mass density

    To invert enthalpy as a function of energy, bounds are required. Default
    is (0,2). Can specify tuple as enthalpybounds=(min,max) if the energies
    you desire are outside this range.  ''' 
    g = float(g) # enforce floating point math on g
    z = p_ref - g * 14.7 # Determine log K by reference pressure p_ref
    K = 10**z
    n = 1 / (g - 1)
    print(g,z, K, n)
    if enthalpybounds:
        (lowenth, highenth) = enthalpybounds
    else:
        (lowenth,highenth) = (0,2.0)
    def density(enthalpy):
        '''rest mass density in g/cm^3 as a function of pseudoenthalpy'''
        eta = np.exp(enthalpy) - 1.
        return (eta / n / K / g)**(n)
    def pressure (enthalpy):
        '''energy density in g/cm^3 as a function of pseudoenthalpy'''
        eta = np.exp(enthalpy) - 1.
        return K * ( eta / n /  K / g )**(n * g)
    def energy(enthalpy):
        '''energy density in g/cm^3 as a function of pseudoenthalpy'''
        # convert from pseudoenthalpy into parameter8 eta = specific enthalpy - 1
        eta = max(np.exp(enthalpy) - 1.,0.)
        if (eta != 0):
            return (1. + eta / g) * (eta / n / K / g )**n
        else:
            return 0.
    def dprden(enthalpy):
        '''dimensionless derivative of pressure with respect to energy 
        as a function of pseudoenthalpy; the speed of sound squared''' 
        eta = max(np.exp(enthalpy) - 1.,0)
        return eta / (1. + eta) / n 
    def enthalpy(energy):
        '''pseudoenthalpy as a function of energy density in g/cm^3'''
        # No nice analytic form
        def diffen(eta):
            return energy - (1. + eta / g) * (eta * (g - 1.) / (K * g) )**n
        higheta = np.exp(highenth) - 1.
        eta = optimize.zeros.brentq(diffen,0.0, higheta)
        return np.log(eta + 1.)
    return (np.vectorize(density), np.vectorize(pressure),
            np.vectorize(energy),  np.vectorize(dprden), np.vectorize(enthalpy))

    
lowsim = np.transpose(np.array([np.loadtxt("lowsim.params")]))
loweos = np.transpose(np.loadtxt("loweos.params"))

def peicewise(lowp, p2, g1, g2, g3):
    ''' Create an upper density equation of state for fiducial densities
    14.7 and 15.0 based on pressures at density 14.7 and gamma in each
    region'''
    g = np.array((g1, g2, g3))
    # Define log of K in each region based on pressure and continuity
    z1 = p2 - g1 * 14.7 
    z2 = p2 - g2 * 14.7
    z3 = z2 + (g2 - g3)*15.0
    zis = np.array((z1, z2, z3))
    # Find dividing density between lowp and core
    x0 = (p2 - g1 * 14.7 - lowp[1][-1]) / (lowp[2][-1] - g1)
    xis = np.array((x0, 14.7, 15.0))

    # Join crust and core
    params = np.hstack((lowp, np.array((xis, zis, g)) ))
    rhois = 10**np.array(params[0])

    # Crust to zero density
    rhois[0]=0.0

    K = 10**np.array(params[1])
    g = np.array(params[2])
    n = 1/(g - 1)

    # thermodynamic consistency term in energy density
    a = [0.0]
    for i in range(len(rhois)-1):
        ai = a[i] + K[i] / (g[i] - 1) * rhois[i+1]**(g[i] - 1)  \
            - K[i+1] / (g[i+1] - 1) * rhois[i+1]**(g[i+1] - 1)
        a.append(ai)
    a = np.array(a)

    etas = a + K * g * n * rhois**(g-1)
    ens = ( (etas - a) / K / g )**n \
              * ( 1 + (a + n*etas) / ( n + 1) )
    def density(enthalpy):
        eta = np.exp(enthalpy) - 1.
        i = etas.searchsorted(eta) - 1
        return ( (eta - a[i]) / (K[i] * (n[i] + 1) ))**n[i]
    density = np.vectorize(density)
    def energy(enthalpy):
        eta = np.exp(enthalpy) - 1.
        i = etas.searchsorted(eta) - 1
        return ( (eta - a[i]) / (K[i] * (n[i] + 1) ))**n[i] \
           * ( 1 + (a[i] + n[i] * eta) / ( n[i] + 1) )
    energy = np.vectorize(energy)
    def pressure(enthalpy):
        eta = np.exp(enthalpy) - 1.
        i = etas.searchsorted(eta) - 1
        return  K[i] * \
            ( ( eta - a[i] ) / ( K[i] * (n[i] + 1) ) )**( n[i] + 1 )
    pressure = np.vectorize(pressure)
    def dprden(enthalpy):
        eta = np.exp(enthalpy) - 1.
        i = etas.searchsorted(eta) - 1
        return  (eta - a[i]) / n[i] / ( 1 + eta )
    dprden = np.vectorize(dprden)
    def enthalpy(energy):
        i = ens.searchsorted(energy) - 1
        if i < len(ens) - 1:
            def diffen(et):
                eten = ((et - a[i]) / (K[i] * (n[i] + 1) ))**n[i] \
                      * ( 1 + (a[i] + n[i] * et) /(n[i]+1))
                return eten - energy
            eta = optimize.zeros.brentq(diffen, etas[i], etas[i+1])
            return np.log(eta + 1.)
        else:
            def diffen(et):
                eten = ((et - a[i]) / (K[i] * (n[i] + 1) ))**n[i] \
                      * ( 1 + (a[i] + n[i] * et) /(n[i]+1))
                return eten - energy
            eta = optimize.zeros.brentq(diffen, etas[i], 1000*etas[i])
            return np.log(eta + 1.)
    enthalpy = np.vectorize(enthalpy)
    return (density, pressure, energy, dprden, enthalpy) 


