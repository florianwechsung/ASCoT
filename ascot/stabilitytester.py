"""
Module containing methods for testing Babuska and Brezzi stability
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import info, warning
from ascot.conditions import *
from ascot.parameters import ascot_parameters
from ascot.infsup import InfSupCollection, StabilityResult

def test_continuity(form, inner_products, spaces):
    cs = [compute_babuska_infsup(form, inner_products, V_h, inverse=True)
          for V_h in spaces]
    cs = StabilityResult(InfSupCollection(cs, "c_h"))
    return cs

def test_stability(forms, inner_products, spaces, bcs=None):
    """
    Test Babuska-Brezzi stability for a variational problem specified
    by (a set of) forms, inner products and a discretization spaces
    parameterized over a mesh family.
    """

    if ascot_parameters["check_continuity"]:
        info("Checking continuity of form")
        continuous = test_continuity(forms, inner_products, spaces)
        if not continuous.is_stable():
            warning("The form does not seems to be bounded!. Check your inner products!")

    if is_saddle_point(forms):
        return _test_brezzi_stability(forms, inner_products, spaces, bcs)

    return _test_babuska_stability(forms, inner_products, spaces, bcs)

def _test_babuska_stability(c, m, spaces, bcs=None):
    """
    For a given form c and an inner product m and discretization
    spaces parameterized over a mesh family, check compute the Babuska
    infsup constants
    """

    info("Testing Babuska stability")

    if bcs is None:
        bcs = [None for W_h in spaces]

    # Compute Babuska constants for each space/mesh
    gamma_hs = [compute_babuska_infsup(c, m, W_h, bc=bc)
                for (W_h, bc) in zip(spaces, bcs)]

    # Collect all constants in a stability result
    gamma_hs = StabilityResult(InfSupCollection(gamma_hs, "gamma_h"))

    # Return stability result
    return gamma_hs

def _test_brezzi_stability((a, b), (m_V, m_Q), spaces, bcs=None):
    """
    For given forms a and b and inner products m_V and m_Q and
    discretization spaces parameterized over a mesh family, check
    compute the Brezzi infsup and coercivity constants
    """

    info("Testing Brezzi stability")

    if bcs is None:
        bcs = [None for W_h in spaces]

    # Compute Brezzi inf-sups
    beta_hs = [compute_brezzi_infsup(b, (m_V, m_Q), W_h, bc)
               for (W_h, bc) in zip(spaces, bcs)]
    beta_hs = InfSupCollection(beta_hs, "beta_h")

    # If inf-sup is non-singular, compute coercivity constant
    if not beta_hs.is_singular():
        alpha_hs = [compute_brezzi_coercivity((a, b), m_V, W_h, bc)
                    for (W_h, bc) in zip(spaces, bcs)]
        alpha_hs = InfSupCollection(alpha_hs, "alpha_h")
    else:
        info("Not computing Brezzi coercivity constants because of singularity")
        alpha_hs = InfSupCollection([], "alpha_h")

    return StabilityResult((beta_hs, alpha_hs))

def is_saddle_point(forms):
    """
    True if two forms are specified.
    """

    if isinstance(forms, (list, tuple)):
        assert(len(forms) == 2)
        return True

    return False

