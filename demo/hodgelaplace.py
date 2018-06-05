__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import *
from ascot import *

# Reduce DOLFIN verbosity
set_log_level(WARNING)

# Define meshes
meshsizes = [2, 4, 8, 16, 32, 64]
meshes = [UnitSquareMesh(n, n) for n in meshsizes]

# Define function spaces (on series of meshes)
Vs = [FunctionSpace(mesh, "CG", 2) * FunctionSpace(mesh, "N1curl", 2)
      for mesh in meshes]

# Define bilinear form
def a((sigma, u), (tau, v)):
    return (- sigma*tau + dot(grad(tau), u)
            + dot(grad(sigma), v) + rot(u)*rot(v))*dx

# Define inner product on mixed space
def m((sigma, u), (tau, v)):
    return (sigma*tau + dot(grad(sigma), grad(tau))
            + dot(u, v) + rot(u)*rot(v))*dx

# Define essential boundary conditions
bcs = [DirichletBC(V, Constant((0.0, 0.0, 0.0)), "on_boundary")
       for V in Vs]

# Test stability (computing Babuska constants)
result = test_stability(a, m, spaces=Vs, bcs=bcs)

# Print result
print result

# Print info about conditions tested and values
for condition in result.conditions:
    print condition
