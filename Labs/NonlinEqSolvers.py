import Derivs

def SolveFPI(F, x0, eps, residuals = None, maxIters = 100):
    """
    Метод простой итерации (Fixed-point iterations).
    Решает уравнение F(x) = 0.
    """
    xn = x0
    iterNum = 0
    while (True):
        xn1 = F(xn)
        r = xn1 - xn

        iterNum += 1
        xn = xn1

        if (residuals is not None):
            residuals.append(r)

        if (abs(r) < eps):
            return xn1

        if (iterNum > maxIters):
            return xn1

def SolveNewthon(F, derF, x0, eps, residuals = None, maxIters = 100):
    newF = lambda xn: xn - F(xn)/derF(xn)
    return SolveFPI(newF, x0, eps, residuals, maxIters)

def SolveNewthon2(F, x0, h, eps, residuals = None, maxIters = 10, warnMaxIters = True):
    xn = x0
    iterNum = 0
    while (True):
        deriv = Derivs.DerivRightO1(F, xn, h)

        if (deriv == 0):
            return SolveFPI(F, xn1, eps, residuals, maxIters - iterNum)
        
        xn1 = xn - F(xn) / deriv
        r = xn1 - xn

        iterNum += 1
        xn = xn1

        if (residuals is not None):
            residuals.append(r)

        if (abs(r) < eps):
            return xn1

        if (iterNum > maxIters):
            if (warnMaxIters):
                print("SolveNewton. Max iters reached.")
            return xn1