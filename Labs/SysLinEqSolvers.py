import numpy
import LinAlg

def GetLDU(A):
    D = numpy.zeros((len(A), len(A)))
    L = numpy.zeros((len(A), len(A)))
    U = numpy.zeros((len(A), len(A)))
    
    D[0][0] = A[0][0]
    for st1 in range(1, len(A)):
        for st2 in range(0, st1):
            L[st1][st2] = A[st1][st2]
            U[st2][st1] = A[st2][st1]
        D[st1][st1] = A[st1][st1]
    return (L, D, U)

def SolveIterate(M1, f2, s0, eps, norm, residuals, maxIters):
    sn = numpy.empty(len(M1))
    numpy.copyto(sn, s0)
    iterNum = 1
    while (True):
        f1  = numpy.dot(M1, sn)
        sn1 = f1 + f2
        r   = norm(sn1 - sn)
        sn = sn1

        if (residuals is not None):
            residuals.append(r)

        if (r < eps):
            return sn
        
        iterNum += 1
        if (iterNum > maxIters):
            return sn

def SolveRelaxation(A, f, s0, eps, tau, norm, residuals : list = None, maxIters = 100):
    (L, D, U) = GetLDU(A)
    
    M1 = numpy.dot(
        -LinAlg.MatrixInverse3(D + tau * L),
        (tau - 1) * D + tau * U)
    f2 = numpy.dot(tau * LinAlg.MatrixInverse3(D + tau * L), f)

    return SolveIterate(M1, f2, s0, eps, norm, residuals, maxIters)

def SolveSeidel(A, f, s0, eps, norm, residuals : list = None, maxIters = 100):
    (L, D, U) = GetLDU(A)
    
    M2 = LinAlg.MatrixInverse1(L + D)
    M1 = -numpy.dot(M2, U)
    f2 = numpy.dot(M2, f)

    return SolveIterate(M1, f2, s0, eps, norm, residuals, maxIters)

def SolveJacob(A, f, s0, eps, norm, residuals = None, maxIters = 100):
    (L, D, U) = GetLDU(A)
    
    M2 = LinAlg.MatrixInverse1(D)
    M1 = -numpy.dot(M2, L + U)
    f2 = numpy.dot(M2, f)

    return SolveIterate(M1, f2, s0, eps, norm, residuals, maxIters)