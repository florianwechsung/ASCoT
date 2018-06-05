"""
ASCOT is an 'Automated Stability COndition Tester'
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

# Check requirements
# FIXME


# Expose the parameters
from ascot.parameters import ascot_parameters

#from dolfin import parameters
#parameters.add(ascot_parameters)

# Expose the test_stability
from ascot.stabilitytester import test_stability

# Expose the individual constant computations
from ascot.conditions import compute_brezzi_infsup, \
     compute_brezzi_coercivity, compute_babuska_infsup

# Expose the stability classes
from ascot.infsup import InfSupCollection, StabilityResult

# Expose the test_stability
from ascot.space_factory import create_spaces, create_inner_products


