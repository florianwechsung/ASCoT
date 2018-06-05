__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import *
from ascot import *

# Reduce DOLFIN verbosity
set_log_level(WARNING)

# Shallow-water
def c((u, p), (v, q)):
    return (-u[1]*v[0] + u[0]*v[1] + div(v)*p + div(u)*q)*dx

# Poisson
#def c((u, p), (v, q)):
#    return (dot(u, v) + div(v)*p + div(u)*q)*dx

def m((u, p), (v, q)):
    return (dot(u, v) + div(u)*div(v)
            + p*q)*dx

# Construct a family of mixed function spaces
meshsizes = [2, 4, 6, 8, 10]
meshes = [UnitSquareMesh(n, n) for n in meshsizes]
W_hs = [FunctionSpace(mesh, "RT", 1) * FunctionSpace(mesh, "DG", 0)
        for mesh in meshes]

# Test stability
result = test_stability(c, m, W_hs)
print W_hs[0].ufl_element().shortstr()
print result
for condition in result.conditions:
    print
    print condition

