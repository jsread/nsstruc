.. _thermodynamics:

Thermodynamics
==============

The structure of the zero-temperature equation of state in terms of
densities comes from the First Law of Thermodynamics:

.. math:: 
  \mathop{d \epsilon} = T \mathop{d s} + \frac{\epsilon + p}{n} \mathop{d n}
  :label: full-first

where :math:`\epsilon` is the energy density, :math:`T` is the temperature,
:math:`s` is the entropy, :math:`p` is pressure density, :math:`n` is number
density. We have set the speed of light :math:`c=1` such that
:math:`\epsilon` and :math:`p` are in the same units.

We set :math:`\mathop{d s} = 0`, and we will also use the average mass per
baryon :math:`m_b` to express in terms of the rest mass density
:math:`\rho` which also has the units of :math:`\epsilon`.

.. math:: 
  \mathop{d \epsilon} = \frac{\epsilon + p}{\rho} \mathop{d \rho}
  :label: first

This implies a single-parameter equation of state, where we express
thermodynamic quantities in terms of a single parameter, for example
:math:`\epsilon(\rho)` and :math:`p(\rho)`; in fact, only one of these is
required, as we can use the first law to determine one from the other. For
example, given :math:`p(\rho)`, we can use the reformulation [#f1]_ of :eq:`first`
as 

.. math::
  d\left(\frac{\epsilon}{\rho}\right) = - p  d\frac{1}{\rho},
  :label: epsilonder

along with the surface boundary condition, to determine
:math:`\epsilon(\rho)`.

We will later make use of the specific enthalpy density :math:`h`,

.. math::
  h = \frac{\epsilon + p}{\rho}.
  :label: enthalpy

One can show [#f2]_ for the single-parameter equation of state that:

.. math::
  \frac{ dh}{dp} = \frac{h}{\epsilon + p}.
  :label: enthalphyder

This enthalpy will vanish at the stellar surface, and given this known boundary
value it is a convenient variable for solving the stellar structure equations.
If we specify :math:`\epsilon(h)` and :math:`p(h)`, then we can determine
:math:`\rho` using:

.. math::
  \rho = (\epsilon + p) e^{-h}.
  :label: densidydef

Polytropes
==============

Polytropes are defined using the equation of state

.. math::
  p = K \rho^\Gamma
  :label: polytrope

and the integral of :eq:`epsilondir` along with the condition that
:math:`\epsilon \rightarrow \rho` as :math:`\rho \rightarrow 0` gives

.. math::
  \epsilon = \rho + \frac{K}{\Gamma - 1}\rho^{\Gamma} 
  :label: polyeps

The specific enthalpy is

.. math::
  h = 1 + \frac{K}{\Gamma - 1} \Gamma \rho^{\Gamma - 1}
  :label: polyenth

which allows the equations for :math:`{\rho, p, \epsilon}` to be inverted:

.. math::
  \rho &= \left(\frac{ (h-1) (\Gamma - 1)}{K \Gamma}\right)^{1 / (\Gamma - 1)}\\
  p &= K \left(\frac{ (h-1) (\Gamma - 1)}{K \Gamma}
                \right)^{\Gamma / (\Gamma - 1)}\\
  \epsilon &= \left(1 + \frac{h-1}{\Gamma}\right) 
        \left(\frac{ (h-1) (\Gamma - 1)}{K \Gamma}\right)^{1 / (\Gamma - 1)}\\
  :label: entheos

  
.. rubric:: Footnotes

.. [#f1] :ref:`apptherm1`
.. [#f2] :ref:`apptherm2`
