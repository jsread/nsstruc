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

.. rubric:: Footnotes

.. [#f1] :ref:`apptherm1`
.. [#f2] :ref:`apptherm2`
