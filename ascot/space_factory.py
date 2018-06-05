"""
Module containing methods for creating a many different element families
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import FunctionSpace, VectorFunctionSpace
from dolfin import dot, inner, grad, curl, div, dx
import itertools

def create_spaces(meshes, value_dimensions, spaces, degrees):
    """
    Create a list of discretization families, each parameterized over
    the mesh based on the given specifications
    """

    families = [_generate_function_spaces(space, dim)[1]
                for space, dim in zip(spaces, value_dimensions)]

    a = itertools.product(*families)
    combos = list(itertools.product(a, degrees))

    W_h = []
    for (space, deg) in combos:
        try:
            W_h += [[space[0](mesh, deg[0]) * space[1](mesh, deg[1])
                     for mesh in meshes]]
        except:
            pass
    return W_h

def create_inner_products(spaces):
    return [_generate_function_spaces(space,0)[0] for space in spaces]

def _generate_function_spaces(regularity, dim):
    """
    Return a list of lambda functions, that, given a mesh and a
    degree, defines a conforming finite element spaces for the given
    space name
    """
    if dim == 1:
        vFunctionSpace = FunctionSpace
    else:
        vFunctionSpace = VectorFunctionSpace

    spaces = {"h1": (lambda u, v: (dot(u,v) + inner(grad(u),grad(v)))*dx,
                     [lambda mesh, deg: vFunctionSpace(mesh, "CG", deg)]),

              "hdiv": (lambda u, v: (dot(u, v) + div(u)*div(v))*dx,
                       [lambda mesh, deg: vFunctionSpace(mesh, "CG", deg),
                        lambda mesh, deg: FunctionSpace(mesh, "RT", deg),
                        lambda mesh, deg: FunctionSpace(mesh, "BDM", deg)]),

              "hcurl": (lambda u, v: (dot(u, v) + inner(curl(u), curl(v)))*dx,
                        [lambda mesh, deg: vFunctionSpace(mesh, "CG", deg),
                         lambda mesh, deg: FunctionSpace(mesh, "N1curl", deg)]),

              "l2": (lambda p, q: dot(p, q)*dx,
                     [lambda mesh, deg: vFunctionSpace(mesh, "CG", deg),
                      lambda mesh, deg: vFunctionSpace(mesh, "DG", deg)])
              }
    return spaces[regularity.lower()]
