__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import *
from ascot import *

# Reduce DOLFIN verbosity
set_log_level(WARNING)

# Define b form:
b = lambda v, q: div(v)*q*dx

# Construct a family of mixed function spaces
meshsizes = list(range(2, 10, 2))
meshes = [UnitSquareMesh(n, n) for n in meshsizes]

# Generate spaces
specifications = {"value_dimensions": (2, 1),
                  "degrees": [(i,j) for i in range(1,5) for j in range(i)],
                  "spaces": ("H1", "L2")}
spaces = create_spaces(meshes, **specifications)
W = create_inner_products(specifications["spaces"])

# Only check for inf-sup stability:
ascot_parameters["only_stable"] = True

stable_elements = []
for W_hs in spaces:
    beta_hs = [compute_brezzi_infsup(b, W, W_h) for W_h in W_hs]
    result = StabilityResult(InfSupCollection(beta_hs, "beta_h"))
    if result.is_stable():
        stable_elements += [W_hs[0].ufl_element()]

print()
print("Number of stable elements: %d" % len(stable_elements))
print("These are: ")
for element in stable_elements:
    print(element.shortstr())
