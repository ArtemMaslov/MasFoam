import NonlinEqSolvers
import math

###################################################################################################################
# Методы Эйлера.
###################################################################################################################

def DifEqEulerSolver1a(f, t, u0):
    """
    Явный метод Эйлера 1 порядка.
    """
    # du/dt = f(t, u)
    # (u_{n} - u_{n - 1}) / tau = f(t_{n - 1}, u_{n - 1}), where t_{n} = tau * n
    # u_{n} = u_{n - 1} + tau * f(t_{n - 1}, u_{n - 1})
    u = []
    u.append(u0) # t[0]
    var = {"u_{n-1}" : 0, "u_{n}" : 0}
    for n in range(1, len(t)):
        var["u_{n-1}"] = u[n - 1]
        tau = t[n] - t[n - 1]
        var["u_{n}"] = var["u_{n-1}"] + tau * f(t[n - 1], var["u_{n-1}"])
        u.append(var["u_{n}"])
    return u

def DifEqEulerSolver1b(f, t, u0):
    """
    Явный метод Эйлера 1 порядка.
    """
    # du/dt = f(t, u)
    # (u_{n} - u_{n - 1}) / tau = f(t_{n - 1}, u_{n - 1}), where t_{n} = tau * n
    # u_{n} = u_{n - 1} + tau * f(t_{n}, u_{n - 1})
    u = []
    u.append(u0) # t[0]
    var = {"u_{n-1}" : 0, "u_{n}" : 0}
    for n in range(1, len(t)):
        var["u_{n-1}"] = u[n - 1]
        tau = t[n] - t[n - 1]
        var["u_{n}"] = var["u_{n-1}"] + tau * f(t[n], var["u_{n-1}"])
        u.append(var["u_{n}"])
    return u

def DifEqEulerSolver2a(f, t, u0, eps=1e-6):
    """
    Неявный метод Эйлера 1 порядка.
    """
    print("DifEqEulerSolver2a")
    # du/dt = f(t, u)
    # (u_{n} - u_{n - 1}) / tau = f(t_{n - 1}, u_{n}), where t_{n} = tau * n
    # u_{n} = u_{n - 1} + tau * f(t_{n}, u_{n})
    u = []
    u.append(u0) # t[0]
    var = {"u_{n-1}" : 0, "u_{n}" : 0}
    for n in range(1, len(t)):
        var["u_{n-1}"] = u[n - 1]
        tau = t[n] - t[n - 1]
        funct = lambda un: var["u_{n-1}"] + tau * f(t[n], un)
        # уравнение: u_{n} = funct(u_{n}) = u_{n-1} + tau * f(t_{n}, u_{n}). Метод Ньютона.
        var["u_{n}"] = NonlinEqSolvers.SolveNewthon2(funct, var["u_{n-1}"], math.sqrt(eps) * var["u_{n-1}"], eps)
        u.append(var["u_{n}"])
    return u

def DifEqEulerSolver2b(f, t, u0, eps=1e-6):
    """
    Неявный метод Эйлера 1 порядка.
    """
    # du/dt = f(t, u)
    # (u_{n} - u_{n - 1}) / tau = f(t_{n - 1}, u_{n}), where t_{n} = tau * n
    # u_{n} = u_{n - 1} + tau * f(t_{n - 1}, u_{n})
    u = []
    u.append(u0) # t[0]
    var = {"u_{n-1}" : 0, "u_{n}" : 0}
    for n in range(1, len(t)):
        var["u_{n-1}"] = u[n - 1]
        tau = t[n] - t[n - 1]
        funct = lambda un: var["u_{n-1}"] + tau * f(t[n - 1], un)
        # уравнение: u_{n} = funct(u_{n}) = u_{n-1} + tau * f(t_{n - 1}, u_{n}). Метод Ньютона.
        var["u_{n}"] = NonlinEqSolvers.SolveNewthon2(funct, var["u_{n-1}"], math.sqrt(eps) * var["u_{n-1}"], eps)
        u.append(var["u_{n}"])
    return u

###################################################################################################################
###################################################################################################################