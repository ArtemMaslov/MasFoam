from collections.abc import Iterable

import numpy

def DerivRightO1(f, x, h):
    """
    (f(x + h) - f(x)) / h
    """
    return (f(x + h) - f(x)) / h

def DerivLeftO1(f, x, h):
    """
    (f(x) - f(x - h)) / h
    """
    return (f(x) - f(x - h)) / h

def DerivCentralO2(f, x, h):
    """
    (f(x + h) - f(x - h)) / (2h)
    """
    return (f(x + h) - f(x - h)) / (2 * h)

def DerivCentralO4(f, x, h):
    """
    2/(3h) * (f(x + h) - f(x - h)) - 1/(12h) * (f(x + 2h) - f(x - 2h))
    """
    return 4 * (f(x +     h) - f(x -     h)) / (3 * 2 * h) - \
               (f(x + 2 * h) - f(x - 2 * h)) / (3 * 4 * h)

def DerivCentralO6(f, x, h):
    """

    """
    return 3 * (f(x +     h) - f(x -     h)) / (2 * 2 * h) - \
           3 * (f(x + 2 * h) - f(x - 2 * h)) / (5 * 4 * h) + \
               (f(x + 3 * h) - f(x - 3 * h)) / (10 * 6 * h)

def CalcJacobiMatrix(F, x, deriv_funct, x_step):
    derF = numpy.zeros((len(x), len(x)))
    for columnIndex in range(0, len(x)):
        partialX = numpy.copy(x)
        def helper(x):
            nonlocal partialX
            partialX[columnIndex] = x
            return F(partialX)
        derColumn = deriv_funct(helper, x[columnIndex], x_step)
        for rowIndex in range(0, len(x)):
            derF[rowIndex][columnIndex] = derColumn[rowIndex]
    return derF