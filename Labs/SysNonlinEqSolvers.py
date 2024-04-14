import numpy

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
        
