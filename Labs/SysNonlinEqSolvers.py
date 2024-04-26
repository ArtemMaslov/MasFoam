import numpy
from collections.abc import Iterable

import LinAlg

def ApplyFunctionVector(funct, xn):
    res = numpy.empty(len(funct))
    for st in range(0, len(res)):
        res[st] = funct[st](xn)
    return res

def ApplyFunctionMatrix(funct, xn):
    res = numpy.empty((len(funct), len(funct)))
    for st1 in range(0, len(res)):
        for st2 in range(0, len(res)):
            res[st1][st2] = funct[st1][st2](xn)
    return res

def SolveSystemFPI(F, x0, eps, norm, residuals = None, maxIters = 100):
    xn = x0
    iterNum = 1
    while (True):
        xn1 = ApplyFunctionVector(F, xn)
        r = norm(xn1 - xn)

        iterNum += 1
        xn = xn1

        if (residuals is not None):
            residuals.append(r)

        if (abs(r) < eps):
            return xn1

        if (iterNum > maxIters):
            return xn1

def SolveSystemNewthon(F, derF, x0, eps, norm, residuals = None, maxIters = 100):
    xn = x0
    iterNum = 1
    while (True):
        A = ApplyFunctionMatrix(derF, xn)
        f = ApplyFunctionVector(F, xn)
        # A*res = f
        
        #res = SolveGauss(A, f)
        #res = numpy.dot(MatrixInverse3(A), -ApplyFunctionVector(F, xn))

        #xn1 = res + xn
        xn1 = xn - numpy.dot(LinAlg.MatrixInverse3(A), f)
        res = xn1 - xn
        r = norm(res)

        if (False):
            print("xn\n", xn)
            print("A\n", A)
            print("f\n", f)
            print("Gauss:\n", numpy.dot(A, res) - f)
            print("res\n", res)
            print("xn1\n", xn1)

        iterNum += 1
        xn = xn1

        if (residuals is not None):
            residuals.append(r)

        if (abs(r) < eps):
            return xn1

        if (iterNum > maxIters):
            return xn1
        
def SolveSystemFPI2(F, x0, eps, norm, residuals = None, maxIters = 30):
    xn = x0
    iterNum = 1
    while (True):
        xn1 = F(xn)
        r = norm(xn1 - xn)

        iterNum += 1
        xn = xn1

        if (residuals is not None):
            residuals.append(r)

        if (abs(r) < eps):
            #print(f"FPI sys solver. ItersCount = {iterNum}")
            return xn1

        if (iterNum > maxIters):
            #print(f"FPI sys solver. Exceed ItersCount = {iterNum}")
            return xn1

def SolveSystemNewton2(F, derF, x0, eps, norm, residuals = None, maxIters = 30):
    xn = x0
    iterNum = 1

    printDebug = False

    while (True):
        A = derF(xn)

        if (printDebug):
            print(f"\nNewton sys solver. derF = {A}\n")

        f = F(xn)

        if (printDebug):
            print(f"\nNewton sys solver. f = {f}\n")

        if (isinstance(A, Iterable)):
            try:
                invA = numpy.linalg.inv(A)
            except numpy.linalg.LinAlgError:
                return SolveSystemFPI2(F, xn, eps, norm, residuals, maxIters)
        else:
            invA = 1/A
        
        if (printDebug):
            print(f"Newton sys solver. inv A = {invA}")
            print(f"Newton sys solver. dx = {-numpy.dot(invA, f)}")
            print(f"Newton sys solver. new x = {xn - numpy.dot(invA, f)}\n")
            print(f"Newton sys solver. res = {xn - numpy.dot(invA, f) - xn}\n")

        xn1 = xn - numpy.dot(invA, f)
        res = xn1 - xn
        r = norm(res)
            
        iterNum += 1
        xn = xn1

        if (residuals is not None):
            residuals.append(r)

        if (abs(r) < eps):
            #print(f"Newthon sys solver. ItersCount = {iterNum}")
            return xn1

        if (iterNum > maxIters):
            print(f"Newthon sys solver. Exceed ItersCount = {iterNum}")
            return xn1