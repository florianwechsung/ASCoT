__author__ = "Jack S. Hale (j.hale09@imperial.ac.uk)"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import *
from ascot import *

set_log_level(WARNING)

def all_boundary(x, on_boundary):
    tol = 1E-14
    return on_boundary

def reissner_mindlin(t):
    t = Constant(t)
    nu = Constant(0.3)

    # Recipe for alpha
    h = CellSize(mesh)
    max_condition = ge(h,t)
    alpha = conditional(max_condition, h**(-1.0)/(5.0*(1.0 - nu)), t**(-1.99))

    # Stabilised Reissner-Mindlin problem
    # Arnold and Brezzi, Boundary Value Problems for Partial Differential Equations and Applications,
    # 1993, 289-292
    def a((theta, z_3, gamma), (eta, y_3, psi)):
        # Small strain operator
        e = lambda theta: 0.5*(grad(theta) + grad(theta).T)
        # Bending operator
        L = lambda e: ((1 - nu)*e + nu*tr(e)*Identity(2))

        # Bending and shear bilinear forms calculated from displacements
        bending = inner(L(e(theta)), e(eta))
        shear = inner(grad(z_3) - theta, grad(y_3) - eta)

        return (bending + alpha*shear + inner(grad(y_3) - eta, gamma) + inner(grad(z_3) - theta, psi)
                - t**2/(1.0 - alpha*t**2)*inner(gamma, psi))*dx

    # WARNING: Norm not bounded for t=0 (due to boundary layers)
    def inner_product((theta, z_3, gamma), (eta, y_3, psi)):
        return (inner(theta, eta) + inner(grad(theta), grad(eta))
                + inner(z_3, y_3) + inner(grad(z_3), grad(y_3))
                + t**2*inner(gamma, psi))*dx

    return (a, inner_product)

mesh_sizes = [4, 8, 12, 16, 20]
meshes = [UnitSquareMesh(n, n) for n in mesh_sizes]

# Chinosi and Lovadina, Computational Mechanics, 16 (1995) 36-44
# Note: Labels for Fig. 2 and 3 are confused in this paper
TRIA0220 = lambda mesh: [VectorFunctionSpace(mesh, "CG", 2, dim=2),
                         FunctionSpace(mesh, "CG", 2),
                         VectorFunctionSpace(mesh, "DG", 0, dim=2)]
# Chinosi and Lovadina, Computational Mechanics, 16 (1995) 36-44, also in Arnold and Brezzi
TRIA1B20 = lambda mesh: [VectorFunctionSpace(mesh, "CG", 1, dim=2) +
                         VectorFunctionSpace(mesh, "Bubble", 3, dim=2),
                         FunctionSpace(mesh, "CG", 2),
                         VectorFunctionSpace(mesh, "DG", 0, dim=2)]
# Arnold and Brezzi, Boundary Value Problems for Partial Differential Equations and Applications,
# 1993, 289-292
MINI1 = lambda mesh: [VectorFunctionSpace(mesh, "CG", 1, dim=2),
                      FunctionSpace(mesh, "CG", 1) + FunctionSpace(mesh, "Bubble", 3),
                      VectorFunctionSpace(mesh, "CG", 1, dim=2)]
# Arnold and Brezzi, Boundary Value Problems for Partial Differential Equations and Applications,
# (1993), 289-292
MINI2 = lambda mesh: [VectorFunctionSpace(mesh, "CG", 2, dim=2),
                      FunctionSpace(mesh, "CG", 2) + FunctionSpace(mesh, "Bubble", 3),
                      VectorFunctionSpace(mesh, "CG", 1, dim=2)]

families = [TRIA0220, TRIA1B20, MINI1, MINI2]
family_names = ["TRIA0220", "TRIA1B20", "MINI1", "MINI2"]

for family, family_name in zip(families, family_names):
    spaces = [MixedFunctionSpace(family(mesh)) for mesh in meshes]
    # Eliminate the three rigid body modes using CCCC boundary conditions
    bcs = [DirichletBC(MixedFunctionSpace(space.split()[0:2]), [0.0,0.0,0.0], all_boundary)
           for space in spaces]

    result = test_stability(*reissner_mindlin(0.001), spaces=spaces, bcs=bcs)
    print "\n%s\n%s" % (spaces[0].ufl_element().shortstr(), str(result))
    for condition in result.conditions:
        print condition

