"""
Module containing methods for computing the Babuska inf-sup and Brezzi
inf-sup and coercivity conditions
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

from math import sqrt
from dolfin import TestFunctions, TrialFunctions
from ascot.eigenvalues import EigenProblem
from ascot.infsup import InfSupConstant
from ascot.parameters import ascot_parameters

__all__ = ["compute_brezzi_infsup", "compute_brezzi_coercivity",
           "compute_babuska_infsup"]

def compute_brezzi_infsup(b, xxx_todo_changeme, W, bc=None):
    """
    For a given form b: V x Q \rightarrow \R and inner products m and
    n defining V and Q respectively and a function space W = V_h x
    Q_h, compute the Brezzi inf-sup constant.
    """
    (m, n) = xxx_todo_changeme
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)
    lhs = m(v, u) + b(v, p) + b(u, q)
    rhs = - n(q, p)

    # Get parameters
    params = ascot_parameters["brezzi_infsup"]
    num = ascot_parameters["number_of_eigenvalues"]

    # Compute eigenvalues
    eigenvalues = EigenProblem(lhs, rhs, params, bc).solve(num)
    return InfSupConstant(W.mesh().hmax(), eigenvalues, operator=sqrt)

def compute_brezzi_coercivity(xxx_todo_changeme1, m, W, bc=None):
    """
    For given forms a: V x V \rightarrow \R and b: V x Q \rightarrow
    \R and an inner product m defining V and a function space W =
    V_h x Q_h, compute the Brezzi coercivity constant.
    """
    (a, b) = xxx_todo_changeme1
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)
    lhs = a(v, u) + b(v, p) + b(u, q)
    rhs = m(v, u)

    # Compute eigenvalues
    params = ascot_parameters["brezzi_coercivity"]
    num = ascot_parameters["number_of_eigenvalues"]
    eigenvalues = EigenProblem(lhs, rhs, params, bc).solve(num)
    return InfSupConstant(W.mesh().hmax(), eigenvalues)

def compute_babuska_infsup(c, m, W, bc=None, inverse=False):
    """
    For a given form c : V x V \rightarrow \R and an inner product m
    defining V and a element space function space W, compute the
    Babuska inf-sup constant.

    The boolean 'inverse' is for ... It defaults to False.
    """

    u = TrialFunctions(W)
    v = TestFunctions(W)

    if len(u) == 1:
        u, v = u[0], v[0]

    lhs = c(v, u)
    rhs = m(v, u)

    if inverse:
        lhs, rhs = rhs, lhs

    # Compute eigenvalues
    params = ascot_parameters["babuska"]
    num = ascot_parameters["number_of_eigenvalues"]
    eigenvalues = EigenProblem(lhs, rhs, params, bc).solve(num)
    return InfSupConstant(W.mesh().hmax(), eigenvalues)
