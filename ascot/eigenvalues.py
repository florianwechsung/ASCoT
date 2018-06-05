"""
EigenProblem wraps the SLEPcEigensolver class for ease of working with
singular generalized eigenvalue problems
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import info, assemble
from dolfin import PETScMatrix, SLEPcEigenSolver

from ascot.parameters import ascot_parameters

class EigenProblem:
    """
    EigenProblem(a, b=None, parameters=None) holds the forms for the
    (generalized eigenvalue problem)

    Ax = \lambda B x

    with B = I if B is None.
    """
    def __init__(self, a, b=None, parameters=None, bc=None):
        """
        Initialize
        """
        self.a = a
        self.b = b
        self.bc = bc
        self.parameters = parameters
        return

    def solve(self, n=None):
        """
        solve(n=None):

        Solve the eigenvalue problem and return the converged
        eigenvalues

        n indicates the number of eigenvalues to be calculated
        """

        # Standard eigenvalue problem
        if self.b is None:
            A = PETScMatrix()
            assemble(self.a, tensor=A)

            if self.bc:
                self.bc.apply(A)

            solver = SLEPcEigenSolver(A)
            solver.solve()

        # Generalized eigenvalue problem
        else:
            A = PETScMatrix()
            assemble(self.a, tensor=A)
            B = PETScMatrix()
            assemble(self.b, tensor=B)

            # Set options for eigensolver
            solver = SLEPcEigenSolver(A, B)
            for key, key_info in self.parameters.iterdata():
                solver.parameters[key] = self.parameters[key]

            if n is None:
                n = A.size(0)

            if self.bc:
                self.bc.apply(A)
                self.bc.apply(B)

            solver.solve(n)

        # Pick real part of eigenvalues computed
        m = solver.get_number_converged()

        complex_eps = 0.001
        eigenvalues = [solver.get_eigenvalue(i)[0] for i in range(m)
                       if abs(solver.get_eigenvalue(i)[1]) < complex_eps]

        if len(eigenvalues) == 0:
            raise RuntimeError("No real-valued eigenvalues found")

        # Compute all eigenvalues if eigenvalue is zero and not only
        # testing for stability
        if (n < A.size(0)
            and abs(eigenvalues[0]) < ascot_parameters["eps"]
            and not ascot_parameters["only_stable"]):

            info("Only zero eigenvalues detected. Computing all.")
            solver.solve(A.size(0))
            m = solver.get_number_converged()
            eigenvalues = [solver.get_eigenvalue(i)[0] for i in range(m)
                           if abs(solver.get_eigenvalue(i)[1]) < 0.1]

        return eigenvalues




