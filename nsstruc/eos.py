"""
.. module:: eos
   :synopsis: Define an eos object for polytropes
.. moduleauthor:: Jocelyn Read <jread@fullerton.edu>

"""

class Eos:
    """ A class that defines cold equations of state for TOV and augmented TOV
    equations.

    Based on eta = exp(enthalpy) - 1.0

    dens, press, energy are in units of g/cm^3
    dens = number density in fm^-3 * m_b for m_b in conversions.py
    eta is dimensionless

    Supplies:
        dens: rest mass density as a function of eta
        press: pressure as a function of eta
        energy: energy density as a function of eta
        dendp: de/dp derivative as a function of eta
        eta: eta (enthalpy as a function of energy
        range: a tuple of (eta_low, eta_hi)
        """
    # FIXME this is not the right way to choose how to initiate the object
    def __init__(self, *args, **kwargs):
     # if len(args) == 1:
     #     self.initFromTable(args[0])
      if len(args) == 5:
          self.initFromFuncs(*args)
     # if len(args) == 2:
     #     self.initFromParam(*args)

    def initFromFuncs(self, densfunc,pressfunc,energyfunc,etafunc, dendp):
        self.densfunc = densfunc
        self.pressfunc = pressfunc
        self.energyfunc = energyfunc
        self.etafunc = etafunc
        self.dendp = dendp

def rns_output(eos, filename):
    '''output Eos object to file format suitable for the rns code'''
    (loweta,higheta) = eos.range
    etas = logspace(log10(loweta),log10(higheta), 400)
    endens = eos.energy(etas)
    press = eos.pressure(etas) * Erg_per_Gram
    numdens = eos.density(etas) / M_b_g
    enth = enthalpy* C**2 *100**2
    out_table = pylab.transpose((endens,press, enth, numdens))
    out_table = out_table[-200:]
    f = open(filename, 'w')
    print >>f, len(out_table)
    print >>f, str(out_table).replace('[',' ').replace(']', ' ')
    f.close()

def createpolytrope(p2, g):
    ''' Create a polytrope using a reference pressure at density 14.7 and
    gamma, store in generic piecewise polytrope format'''
    z = p2 - g * 14.7 # Determine K  in log scale
    return transpose(array([[-2.,z,g]]))


def expandedparams(params):
  """ given a parameter set of [[log rhodivs][log K][gamma]]
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
