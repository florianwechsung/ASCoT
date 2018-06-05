__author__ = "Joachim B. Haga (jobh@simula.no), Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

# Illustrate three-field formulation arising from a Biot consolidation
# problem

from dolfin import *
from ascot import *

# Reduce DOLFIN verbosity
set_log_level(WARNING)

def biot_three_field(epsilon):

    epsilon = Constant(epsilon)
    lamda = 1.0

    def a((v,q,s), (u,p,r)):
        constitutive = dot(v, u) + lamda*div(v)*div(u) + inner(grad(v), grad(u))
        return (constitutive + div(v)*p
                + q*div(u) + q*div(r)
                + epsilon*div(s)*p + dot(s, r))*dx

    def inner_product((v,q,s), (u,p,r)):
        return (inner(v, u) + inner(grad(v), grad(u))
                + q*p
                + inner(s, r) + epsilon*div(s)*div(r))*dx

    return (a, inner_product)

# Reduce spectral shift parameter from default
# ascot_parameters["babuska"]["spectral_shift"] = 0.00001

# Construct a collection of meshes
meshsizes = [4, 8, 12, 16, 20]
meshes = [UnitSquareMesh(n, n) for n in meshsizes]

# Construct a family of function spaces
generate_space = lambda mesh: [VectorFunctionSpace(mesh, "CG", 2),
                               FunctionSpace(mesh, "DG", 0),
                               FunctionSpace(mesh, "RT", 1)]
spaces = [MixedFunctionSpace(generate_space(mesh)) for mesh in meshes]

# Check a variety of parameters
epsilons = [1.0, 1.e-2, 1.e-4, 1.e-6, 0.0]
#epsilons = [1.e-4]

for epsilon in epsilons:
    result = test_stability(*biot_three_field(epsilon), spaces=spaces)
    print "epsilon = %g: %s " % (epsilon, str(result))
    for condition in result.conditions:
        print condition

