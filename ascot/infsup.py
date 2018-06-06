"""
Classes for handing stability characterization and collections of
inf-sup constants
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU LGPL version 3 or any later version"


from math import log
from numpy import array, where

from ascot.parameters import ascot_parameters
from dolfin import warning, info, DEBUG, INFO

def formatted_list(values, formatting):
    """
    For a given list of numbers ('values'), return a pretty string
    representation where each number is formatted according to the
    specified 'formatting'
    """
    return "[%s]" %  ", ".join([formatting % v for v in values])

def filter_values(values, operator):
    """
    Remove inf's and nan's, and apply operator to each value in
    values. Return sorted array of nonzero values, and the number
    of zeroes.
    """
    values = array(values)

    # Remove infs and nans
    upper_bound = ascot_parameters["inf"]
    values = values[where(abs(values) < upper_bound)]
    if any(abs(values) > 1.e8):
        warning("Some eigenvalues are larger than 1.e8. Consider checking your inner products!")

    # Replace small (possibly negative numbers) by zero!
    is_nonzero = abs(values) >= ascot_parameters["eps"]
    zeros = values[abs(values) < ascot_parameters["eps"]]
    values = values[where(is_nonzero)]

    # Apply operator (sqrt) if indicated
    values = array(list(map(operator, values)))

    # Sort values
    values.sort()
    return values, sum(is_nonzero==False)

class InfSupConstant(float):
    """
    An InfSupConstant represents ...

    """
    def __new__(cls, h, values, operator=None):
        "Initialize."

        values, num_singularities = filter_values(values, operator)

        if num_singularities > 0:
            inst = float.__new__(cls, 0.0)
        else:
            inst = float.__new__(cls, min(abs(values)))

        inst.h = h
        inst.values = values
        inst.num_singularities = num_singularities
        return inst

    def reduced_constant(self):
        "Return smallest non-zero eigenvalue or 0 if no such."

        if len(self.values):
            return min(abs(self.values))
        return 0.0

    def singularities(self):
        "Return number of singular values (zeros)."
        return self.num_singularities

    def is_singular(self):
        "Return true if there are singularities (zero values)."
        return self.num_singularities > 0

    def __repr__(self):
        s = "InfSupConstant(%s, %s)" % (str(self.h), str(self.values))
        return s

class InfSupCollection:
    """
    An InfSupCollection represents ...

    """
    def __init__(self, constants, label=""):
        "Initialize."

        self.constants = constants
        self.label = label

    def is_stable(self):
        "True if all constants are non-zero and not decaying."

        if ascot_parameters["only_stable"] and self.is_singular():
                return False
        return (not self.is_singular() and not self.decays())

    def is_singular(self):
        "True if any of the constants is singular."

        return any(c_h.is_singular() for c_h in self.constants)

    def decay_rates(self):
        """
        Rates of decay for the reduced constant
        """

        # FIXME: Make this safer.
        try:
            log_vs = [log(c_h.reduced_constant(), 2) for c_h in self.constants]
            log_hs = [log(c_h.h, 2) for c_h in self.constants]
            rates = [(log_vs[i] - log_vs[i+1])/(log_hs[i] - log_hs[i+1])
                     for i in range(len(log_vs) - 1)]
        except:
            rates = []
        return array(rates)

    def decays(self):
        """
        The collection of constants decay if any of the rates of decay
        is higher than 1.0, or if any rate is higher than 0.3
        ("magic_rate") and increasing.
        """

        # FIXME: Recheck this logic.
        rates = self.decay_rates()

        if len(rates) == 0:
            return False

        # If all rates are less than a certain magic rate, the numbers
        # do NOT decay too much.
        if all(rates < ascot_parameters["magic_rate"]):
            return False

        if rates[-1] >= 1:
            return True

        diffs = rates[:-1] - rates[1:]
        rates_going_down = (all([d > 0 for d in diffs])
                            and all([r < 1 for r in rates]))

        if rates_going_down:
            return False

        return True


    def __str__(self):
        """ Pretty print"""

        if len(self.constants) == 0:
            return "Empty InfSupCollection: %s" % self.label

        s = "InfSupCollection: %s" % self.label
        if self.is_singular():

            (singularities, reduced) = list(zip(*((c.singularities(),
                                              c.reduced_constant())
                                             for c in self.constants)))
            s += "\nsingularities = %s" % formatted_list(singularities, "%d")
            s += "\nreduced = \t %s" % formatted_list(reduced, "%0.5g")

        else:
            s += "\nvalues = \t %s" % formatted_list(self.constants, "%0.5g")
        s += "\nrates = \t %s" %  formatted_list(self.decay_rates(), "%0.3g")

        return s

# Colors (Copied from UFL)
RED   = "\033[1;37;31m%s\033[0m"
BLUE  = "\033[1;37;34m%s\033[0m"
GREEN = "\033[1;37;32m%s\033[0m"

class StabilityResult:
    """
    A StabilityResult ...
    """

    def __init__(self, conditions):
        """
        Each condition is an InfSupCollection
        """
        if not isinstance(conditions, (list, tuple)):
            conditions = (conditions,)
        self.conditions = conditions
        return

    def confidence(self):
        return "Very confident."

    def is_stable(self):
        """ True if all constants are non-singular and non-decaying"""

        return all([condition.is_stable() for condition in self.conditions])

    def __str__(self):
        """ Pretty print"""

        if self.is_stable():
            return GREEN % "Discretization family is: Stable."

        s = "Discretization family is: Unstable. "
        if ascot_parameters["only_stable"]:
            return RED % s

        (singular, decaying) = list(zip(*((c.is_singular(), c.decays())
                                     for c in self.conditions)))
        if any(singular):
            s += "Singular. "
        if any(decaying):
            s += "Decaying. "
        else:
            s += "Reduced stable. "

        return RED % s

