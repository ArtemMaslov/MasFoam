from collections.abc import Iterable

import numpy

def DerivRightO1(f, x, h):
    """
    (f(x + h) - f(x)) / h
    """
    if (h == 0):
        raise Exception()
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
    return 3 * (f(x +     h) - f(x -     h)) / (2 * 2 * h) - \
           3 * (f(x + 2 * h) - f(x - 2 * h)) / (5 * 4 * h) + \
               (f(x + 3 * h) - f(x - 3 * h)) / (10 * 6 * h)

def DerivRightO2(f, x, h):
    return (-3/2 * f(x) + 2 * f(x + h) - 1/2 * f(x + 2*h)) / h

def DerivLeftO2(f, x, h):
    return (+3/2 * f(x) - 2 * f(x - h) + 1/2 * f(x - 2*h)) / h

def Deriv2CentralO2(f, x, h):
    return (f(x + h) - 2 * f(x) + f(x - h)) / (2*h**2)

def Deriv2RightO2(f, x, h):
    return (2 * f(x) - 5 * f(x + h) + 4 * f(x + 2*h) - f(x + 3*h)) / h**2

def Deriv2LeftO2(f, x, h):
    return (2 * f(x) - 5 * f(x - h) + 4 * f(x - 2*h) - f(x - 3*h)) / h**2

def DerivCentralO2A(f, i, h):
    """
    (f[i + 1] - f[i - 1]) / (2h)
    """
    return (f[i + 1] - f[i - 1]) / (2 * h)

def DerivRightO2A(f, i, h):
    return (-3/2 * f[i] + 2 * f[i + 1] - 1/2 * f[i + 2]) / h

def DerivLeftO2A(f, i, h):
    return (+3/2 * f[i] - 2 * f[i - 1] + 1/2 * f[i - 2]) / h

def Deriv2CentralO2A(f, i, h):
    return (f[i + 1] - 2 * f[i] + f[i - 1]) / (h**2)

def Deriv2RightO2A(f, i, h):
    return (2 * f[i] - 5 * f[i + 1] + 4 * f[i + 2] - f[i + 3]) / h**2

def Deriv2LeftO2A(f, i, h):
    return (2 * f[i] - 5 * f[i - 1] + 4 * f[i - 2] - f[i - 3]) / h**2

def CalcJacobiMatrix(F, x, deriv_funct, x_step, dtype = numpy.float64):
    if (not isinstance(x, Iterable)):
        return deriv_funct(F, x, x_step)

    derF = numpy.zeros((len(x), len(x)), dtype=dtype)
    for columnIndex in range(0, len(x)):
        partialX = numpy.copy(x)
        def helper(x):
            nonlocal partialX
            partialX[columnIndex] = x
            #print(f"Jacobi: f{partialX} == {F(partialX)}")
            return F(partialX)
        derColumn = deriv_funct(helper, x[columnIndex], x_step)
        for rowIndex in range(0, len(x)):
            derF[rowIndex][columnIndex] = derColumn[rowIndex]
    return derF