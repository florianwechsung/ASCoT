__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import *
from ascot import *

# Reduce DOLFIN verbosity
set_log_level(WARNING)

def _a(v, u):

    # Define DG parameters
    alpha = 4.0
    gamma = 8.0
    mesh = v.function_space().mesh()
    n = FacetNormal(mesh)
    h = CellSize(mesh)

    # Define bilinear form
    a = dot(grad(v), grad(u))*dx \
        - dot(avg(grad(v)), jump(u, n))*dS \
        - dot(jump(v, n), avg(grad(u)))*dS \
        + alpha/avg(h)*dot(jump(v, n), jump(u, n))*dS \
        - dot(grad(v), u*n)*ds \
        - dot(v*n, grad(u))*ds \
        + gamma/h*v*u*ds
    return a

def H1_dg(v, u):

    mesh = v.function_space().mesh()
    n = FacetNormal(mesh)
    h = CellSize(mesh)

    # Define DG inner product
    m = (inner(grad(u), grad(v)))*dx \
        + 1.0/(avg(h))*dot(jump(v, n), jump(u, n))*dS \
        + 1.0/h*dot(v*n, u*n)*ds
    return m

# Construct a family of function spaces
meshsizes = [2, 4, 8, 16, 32]
meshes = [UnitSquareMesh(n, n) for n in meshsizes]
W_hs = [FunctionSpace(mesh, "DG", 1) for mesh in meshes]

# Enable continuity checking of form:
ascot_parameters["check_continuity"] = True

# Test stability
result = test_stability(_a, H1_dg, W_hs)
print(W_hs[0].ufl_element().shortstr())
print(result)
for condition in result.conditions:
    print()
    print(condition)
