"""
Code examples used in "Automated testing of saddle point stability
conditions", Chapter 36 of the FEniCS book.
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) Marie E. Rognes, 2010 --"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import *
from ascot import *

def case1():

    a = lambda u, v: dot(u, v)*dx
    b = lambda v, q: div(v)*q*dx

    Hdiv = lambda u, v: (dot(u, v) + div(u)*div(v))*dx
    L2 = lambda p, q: dot(p, q)*dx

    meshsizes = [2, 4, 6, 8, 10]
    meshes = [UnitSquareMesh(n, n) for n in meshsizes]
    W_hs = [VectorFunctionSpace(mesh, "CG", 1)*FunctionSpace(mesh, "CG", 1)
            for mesh in meshes]

    result = test_stability((a, b), (Hdiv, L2), W_hs)

    print(W_hs[0].ufl_element())
    print(result)
    for condition in result.conditions:
        print(condition)

def case2():

    a = lambda u, v: dot(u, v)*dx
    b = lambda v, q: div(v)*q*dx

    Hdiv = lambda u, v: (dot(u, v) + div(u)*div(v))*dx
    L2 = lambda p, q: dot(p, q)*dx

    meshsizes = [2, 4, 6, 8, 10]
    meshes = [UnitSquareMesh(n, n) for n in meshsizes]
    W_hs = [VectorFunctionSpace(mesh, "CG", 1)*FunctionSpace(mesh, "DG", 0)
            for mesh in meshes]

    result = test_stability((a, b), (Hdiv, L2), W_hs)

    print(W_hs[0].ufl_element())
    print(result)
    for condition in result.conditions:
        print(condition)

def case3():

    a = lambda u, v: dot(u, v)*dx
    b = lambda v, q: div(v)*q*dx

    Hdiv = lambda u, v: (dot(u, v) + div(u)*div(v))*dx
    L2 = lambda p, q: dot(p, q)*dx

    meshsizes = [2, 4, 6, 8, 10]
    meshes = [UnitSquareMesh(n, n) for n in meshsizes]
    W_hs = [FunctionSpace(mesh, "RT", 1)*FunctionSpace(mesh, "DG", 0)
            for mesh in meshes]

    result = test_stability((a, b), (Hdiv, L2), W_hs)

    print(W_hs[0].ufl_element())
    print(result)
    for condition in result.conditions:
        print(condition)

def case4():

    a = lambda u, v: dot(u, v)*dx
    b = lambda v, q: div(v)*q*dx

    Hdiv = lambda u, v: (dot(u, v) + div(u)*div(v))*dx
    L2 = lambda p, q: dot(p, q)*dx

    meshsizes = [2, 4, 6, 8, 10]
    meshes = [UnitSquareMesh(n, n, "crossed") for n in meshsizes]
    W_hs = [VectorFunctionSpace(mesh, "CG", 1)*FunctionSpace(mesh, "DG", 0)
            for mesh in meshes]

    result = test_stability((a, b), (Hdiv, L2), W_hs)

    print(W_hs[0].ufl_element())
    print(result)
    for condition in result.conditions:
        print(condition)

def case5():
    meshsizes = [2, 4, 6, 8, 10]
    meshes = [UnitSquareMesh(n, n) for n in meshsizes]
    specifications = {"value_dimensions": (2, 1),
                      "degrees":((i, j) for i in range(1, 5) for j in range(i)),
                      "spaces": ("H1", "L2")}
    spaces = create_spaces(meshes, **specifications)

    b = lambda v, q: div(v)*q*dx
    H1 = lambda u, v: (dot(u, v) + inner(grad(u), grad(v)))*dx
    L2 = lambda p, q: dot(p, q)*dx

    ascot_parameters["only_stable"] = True

    stable_elements = []
    for W_hs in spaces:
        beta_hs = [compute_brezzi_infsup(b, (H1, L2), W_h) for W_h in W_hs]
        result = StabilityResult(InfSupCollection(beta_hs, "beta_h"))
        if result.is_stable():
            stable_elements += [W_hs[0].ufl_element()]

    print()
    print("Number of stable elements: %d" % len(stable_elements))
    print("These are: ")
    for element in stable_elements:
        print(element.shortstr())

if __name__ == "__main__":

    case1()
    #case2() # glibc detected in SLEPc 3.1-p6 and later versions
    #case3() # glibc detected in SLEPc 3.1-p6 and later versions
    case4()
    case5()

