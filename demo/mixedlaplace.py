__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import *
from ascot import *

# Reduce DOLFIN verbosity
set_log_level(WARNING)

# Define a and b forms:
a = lambda u, v: dot(u, v)*dx
b = lambda v, q: div(v)*q*dx

# Define inner products:
inners = create_inner_products(['Hdiv', 'L2'])

# Construct a family of mixed function spaces
meshsizes = [2, 4, 6, 8, 10]
meshes = [UnitSquareMesh(n, n) for n in meshsizes]
W_hs = [VectorFunctionSpace(mesh, "CG", 1) * FunctionSpace(mesh, "CG", 1)
        for mesh in meshes]

# Test stability
result = test_stability((a, b), inners, W_hs)
print W_hs[0].ufl_element().shortstr()
print result
for condition in result.conditions:
    print
    print condition
