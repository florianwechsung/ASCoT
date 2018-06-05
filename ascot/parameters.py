"""
ascot specific parameters
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import Parameters

ascot_parameters = Parameters("ascot")
ascot_parameters.add("eps", 1.e-05)
ascot_parameters.add("magic_rate", 0.1)
ascot_parameters.add("only_stable", False)
ascot_parameters.add("check_continuity", False)
ascot_parameters.add("inf", 1.e10)
ascot_parameters.add("number_of_eigenvalues", 1)

# Eigensolver parameters for the different conditions:
bip = Parameters("brezzi_infsup")
bip.add("solver", "krylov-schur")
bip.add("spectral_transform", "shift-and-invert")
bip.add("spectrum", "target magnitude")
bip.add("spectral_shift", -0.1)

bcp = Parameters("brezzi_coercivity")
bcp.add("solver", "lapack")

bab = Parameters("babuska")
bab.add("solver", "krylov-schur")
bab.add("spectral_transform", "shift-and-invert")
bab.add("spectral_shift", 1.e-06)
bab.add("spectrum", "target magnitude")

ascot_parameters.add(bip)
ascot_parameters.add(bcp)
ascot_parameters.add(bab)

