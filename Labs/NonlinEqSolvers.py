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
    """
    Метод простой итерации (Fixed-point iterations).
    Решает уравнение x = F(x).
    """
    newF = lambda xn: xn - F(xn)/derF(xn)
    return SolveFPI(newF, x0, eps, residuals, maxIters)

def SolveNewthon2(F, x0, h, eps, residuals = None, maxIters = 10):
    """
    Метод простой итерации (Fixed-point iterations).
    Решает уравнение x = F(x).
    """
    xn = x0
    iterNum = 0
    print("SolveNewthon:")
    while (True):
        print(f"h = {h}, xn= {xn}, F(xn) = {F(xn)}, der = {Derivs.DerivRightO1(F, xn, h)}")
        print(f"F(xn+1) = {F(xn + h)}, dF = {F(xn + h) - F(xn)}, frac = {(F(xn + h) - F(xn)) / h}")
        xn1 = xn - F(xn) / Derivs.DerivRightO1(F, xn, h)
        r = xn1 - xn

        iterNum += 1
        xn = xn1

        print(f"\t#{iterNum}: xn = {xn1}, r = {abs(r)}")

        if (residuals is not None):
            residuals.append(r)

        if (abs(r) < eps):
            return xn1

        if (iterNum > maxIters):
            return xn1