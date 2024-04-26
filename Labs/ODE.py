from collections.abc import Iterable
import numpy
import math

import Norms
import SysNonlinEqSolvers
import Derivs

def RungeKuttaExplicit(t0, y0, tau, stepsCount, funct, ButchersTable, ButchersY, ButchersTau):
    """
    Решает ОДУ или систему ОДУ методом Рунге-Кутты с помощью таблицы Бутчера.

    Аргументы:
        t0, y0     - начальное приближение (скаляры для ОДУ, векторы для системы).
        tau        - шаг по времени.
        stepsCount - число шагов по времени.
        funct      - функция, задающая уравнение, или вектор-функция для системы.
        ButchersTable, ButchersY, ButchersTau - таблица Бутчера:
        B |
        u |
        t |
        c |
        h |
        e |       ButchersTable
        r |
        s |
        T |
        a |
        u |
        ------------------------------
          |         ButchersY
        Пример:
            Метод Эйлера 1 порядка:
                0 | 0
                -----
                  | 1

                ButchersTable = [[0]]
                ButchersTau   =  [0]
                ButchersY     =  [1]

            Метод Рунге-Кутты 3 порядка:
                0   |  0   0   0
                1/2 | 1/2  0   0
                1   |  0   1   0
                -----------------
                    | 1/6 2/3 1/6

                ButchersTable = [[ 0,   0,   0], 
                                 [1/2,  0,   0],
                                 [ 0    1,   0]]
                ButchersTau   =  [ 0,  1/2,  1]
                ButchersY     =  [1/6, 2/3, 1/6]

    Возвращает t, y - решение (либо массивы скаляров для ОДУ, либо массивы векторов для системы ОДУ).
    """

    assert(len(ButchersTable) == len(ButchersY))
    assert(len(ButchersTable[0]) == len(ButchersY))
    assert(len(ButchersTau) == len(ButchersY))

    tn = t0
    yn = y0
    t = [t0]
    y = [y0]

    printDebug = False

    stagesCount = len(ButchersY)
    for n in range(0, stepsCount):
        
        if (printDebug):
            print(f"\nRungeKutta. t = {tn}\n\n")

        k = []
        for st in range(0, stagesCount):
            ksum = 0
            for st1 in range(0, st):
                ksum += ButchersTable[st][st1] * k[st1]
            k = numpy.array(list(k) + [funct(tn + ButchersTau[st] * tau, yn + tau * ksum)])
            if (printDebug):
                print("k:\n", k)
                print("ksum:\n", ksum)
                print("funct:\n", funct(tn + ButchersTau[st] * tau, yn + tau * ksum))
                

        if (printDebug):
            print("k:\n", k)
        ysum = 0
        for st in range(0, len(ButchersY)):
            ysum += ButchersY[st] * k[st]
            
        yn1 = yn + tau * ysum
        tn1 = tn + tau

        t.append(tn1)
        y.append(yn1)

        tn = tn1
        yn = yn1

    return numpy.array(t), numpy.array(y)

def Euler1(t0, y0, tau, stepsCount, funct):
    ButchersTable = numpy.array([[0]])
    ButchersTau   = numpy.array([0])
    ButchersY     = numpy.array([1])
    return RungeKuttaExplicit(t0, y0, tau, stepsCount, funct, ButchersTable, ButchersY, ButchersTau)

def Euler2(t0, y0, tau, stepsCount, funct):
    ButchersTable = numpy.array([[0, 0],
                                 [1/2, 0]])
    ButchersTau   = numpy.array([0, 1/2])
    ButchersY     = numpy.array([0, 1])
    return RungeKuttaExplicit(t0, y0, tau, stepsCount, funct, ButchersTable, ButchersY, ButchersTau)

def Heun3(t0, y0, tau, stepsCount, funct):
    ButchersTable = numpy.array([[0, 0, 0],
                                 [1/3, 0, 0],
                                 [0, 2/3, 0]])
    ButchersTau   = numpy.array([0, 1/3, 2/3])
    ButchersY     = numpy.array([1/4, 0, 3/4])
    return RungeKuttaExplicit(t0, y0, tau, stepsCount, funct, ButchersTable, ButchersY, ButchersTau)

def RungeKutta3(t0, y0, tau, stepsCount, funct):
    ButchersTable = numpy.array([[0, 0, 0],
                                 [1/2, 0, 0],
                                 [0, 1, 0]])
    ButchersTau   = numpy.array([0, 1/2, 1])
    ButchersY     = numpy.array([1/6, 2/3, 1/6])
    return RungeKuttaExplicit(t0, y0, tau, stepsCount, funct, ButchersTable, ButchersY, ButchersTau)

def RungeKutta4(t0, y0, tau, stepsCount, funct):
    ButchersTable = numpy.array([[0, 0, 0, 0],
                                 [1/2, 0, 0, 0],
                                 [0, 1/2, 0, 0],
                                 [0, 0, 1, 0]])
    ButchersTau   = numpy.array([0, 1/2, 1/2, 1])
    ButchersY     = numpy.array([1/6, 2/6, 2/6, 1/6])
    return RungeKuttaExplicit(t0, y0, tau, stepsCount, funct, ButchersTable, ButchersY, ButchersTau)

def Butchers6(t0, y0, tau, stepsCount, funct):
    ButchersTable = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                                 [1/2, 0, 0, 0, 0, 0, 0],
                                 [2/9, 4/9, 0, 0, 0, 0, 0],
                                 [7/36, 2/9, -1/12, 0, 0, 0, 0],
                                 [-35/144, -55/36, 35/48, 15/8, 0, 0, 0],
                                 [-1/360, -11/36, -1/8, 1/2, 1/10, 0, 0],
                                 [-41/260, 22/13, 43/156, -118/39, 32/195, 80/39, 0]])
    ButchersTau   = numpy.array([0, 1/2, 2/3, 1/3, 5/6, 1/6, 1])
    ButchersY     = numpy.array([13/200, 0, 11/40, 11/40, 4/25, 4/25, 13/200])
    return RungeKuttaExplicit(t0, y0, tau, stepsCount, funct, ButchersTable, ButchersY, ButchersTau)

def Butchers6(t0, y0, tau, stepsCount, funct):
    ButchersTable = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                                 [1/2, 0, 0, 0, 0, 0, 0],
                                 [2/9, 4/9, 0, 0, 0, 0, 0],
                                 [7/36, 2/9, -1/12, 0, 0, 0, 0],
                                 [-35/144, -55/36, 35/48, 15/8, 0, 0, 0],
                                 [-1/360, -11/36, -1/8, 1/2, 1/10, 0, 0],
                                 [-41/260, 22/13, 43/156, -118/39, 32/195, 80/39, 0]])
    ButchersTau   = numpy.array([0, 1/2, 2/3, 1/3, 5/6, 1/6, 1])
    ButchersY     = numpy.array([13/200, 0, 11/40, 11/40, 4/25, 4/25, 13/200])
    return RungeKuttaExplicit(t0, y0, tau, stepsCount, funct, ButchersTable, ButchersY, ButchersTau)

def RungeKuttaImlicit(t0, y0, tau, stepsCount, funct, ButchersTable, ButchersY, ButchersTau, eps=1e-6):
    """
    Решает ОДУ или систему ОДУ неявным методом Рунге-Кутты с помощью таблицы Бутчера.

    Аргументы:
        t0, y0     - начальное приближение (скаляры для ОДУ, векторы для системы).
        tau        - шаг по времени.
        stepsCount - число шагов по времени.
        funct      - функция, задающая уравнение, или вектор-функция для системы.
        ButchersTable, ButchersY, ButchersTau - таблица Бутчера:
        B |
        u |
        t |
        c |
        h |
        e |       ButchersTable
        r |
        s |
        T |
        a |
        u |
        ------------------------------
          |         ButchersY
        Пример:
            Метод Эйлера 1 порядка:
                0 | 0
                -----
                  | 1

                ButchersTable = [[0]]
                ButchersTau   =  [0]
                ButchersY     =  [1]

            Метод Рунге-Кутты 3 порядка:
                0   |  0   0   0
                1/2 | 1/2  0   0
                1   |  0   1   0
                -----------------
                    | 1/6 2/3 1/6

                ButchersTable = [[ 0,   0,   0], 
                                 [1/2,  0,   0],
                                 [ 0    1,   0]]
                ButchersTau   =  [ 0,  1/2,  1]
                ButchersY     =  [1/6, 2/3, 1/6]

    Возвращает t, y - решение (либо массивы скаляров для ОДУ, либо массивы векторов для системы ОДУ).
    """

    assert(len(ButchersTable) == len(ButchersY))
    assert(len(ButchersTable[0]) == len(ButchersY))
    assert(len(ButchersTau) == len(ButchersY))
    assert(not isinstance(y0, Iterable))

    tn = t0
    yn = y0
    t = [t0]
    y = [y0]

    printDebug = False

    stagesCount = len(ButchersY)
    for n in range(0, stepsCount):
        if (printDebug):
            print(f"\nRungeKutta. t = {tn}\n\n")
            
        def F_helper(k):
            nonlocal tn
            nonlocal yn
            delta_k = []
            for st in range(0, stagesCount):
                ksum = 0
                for st1 in range(0, stagesCount):
                    if (printDebug):
                        print(f"F_helper: k[st1] == {k[st1]}")
                    ksum += ButchersTable[st][st1] * k[st1]
                delta_k.append(k[st] - funct(tn + ButchersTau[st] * tau, yn + tau * ksum))
                if (printDebug):
                    print(f"F_helper: ksum == {ksum}")
                    print(f"F_helper: funct({tn + ButchersTau[st] * tau}, {yn + tau * ksum})")
                    print(f"F_helper: delta_k = {delta_k}\n")
            if (printDebug):
                print(f"F_helper: delta_k = {delta_k}\n")
            return numpy.array(delta_k)
        
        def derF_helper(k):
            return Derivs.CalcJacobiMatrix(F_helper, k, Derivs.DerivCentralO2, tau/10)
        
        k = SysNonlinEqSolvers.SolveSystemNewton2(F_helper, derF_helper, numpy.ones(stagesCount), eps, Norms.NormV1, maxIters=20)

        ysum = 0
        for st in range(0, len(ButchersY)):
            ysum += ButchersY[st] * k[st]
        
        yn1 = yn + tau * ysum
        tn1 = tn + tau

        t.append(tn1)
        y.append(yn1)

        tn = tn1
        yn = yn1
    return numpy.array(t), numpy.array(y)

def RungeKuttaImlicitSystem(t0, y0, tau, stepsCount, funct, ButchersTable, ButchersY, ButchersTau, eps=1e-6):
    assert(len(ButchersTable) == len(ButchersY))
    assert(len(ButchersTable[0]) == len(ButchersY))
    assert(len(ButchersTau) == len(ButchersY))
    assert(isinstance(y0, Iterable)) # переменная должна быть массивом переменных.

    tn = t0
    yn = y0
    t = [t0]
    y = [y0]
    var_dim = len(y0)

    printDebug = False

    def GetVect(k):
        veck = []
        for st in range(0, stagesCount):
            ar = []
            for st1 in range(0, var_dim):
                ar.append(k[st * var_dim + st1])
            veck.append(ar)
        return numpy.array(veck)

    stagesCount = len(ButchersY)
    for n in range(0, stepsCount):
        if (printDebug):
            print(f"\nRungeKutta. t = {tn}\n\n")

        def F_helper(k):
            delta_k = []
            k = GetVect(k)
            for st in range(0, stagesCount):
                ksum = numpy.zeros(var_dim)
                for st1 in range(0, stagesCount):
                    if (printDebug):
                        print(f"F_helper: k[st1] == {k[st1]}")
                    ksum += ButchersTable[st][st1] * k[st1]

                elem = k[st] - funct(tn + ButchersTau[st] * tau, yn + tau * ksum)
                delta_k += list(elem)

                if (printDebug):
                    print(f"F_helper: ksum == {ksum}")
                    print(f"F_helper: funct({tn + ButchersTau[st] * tau}, {yn + tau * ksum}) = {funct(tn + ButchersTau[st] * tau, yn + tau * ksum)}")
                    print(f"F_helper: delta_k = {delta_k}\n")
            if (printDebug):
                print(f"F_helper: delta_k = {delta_k}\n")
            return numpy.array(delta_k)
        
        def derF_helper(k):
            return Derivs.CalcJacobiMatrix(F_helper, k, Derivs.DerivCentralO2, tau/10)
        
        kvec = SysNonlinEqSolvers.SolveSystemNewton2(F_helper, derF_helper, numpy.zeros(stagesCount * var_dim), eps, Norms.NormV1, maxIters=20)
        k = GetVect(kvec)

        ysum = 0
        for st in range(0, stagesCount):
            ysum += ButchersY[st] * k[st]
        yn1 = yn + tau * ysum
        tn1 = tn + tau

        t.append(tn1)
        y.append(yn1)

        tn = tn1
        yn = yn1
    return numpy.array(t), numpy.array(y)

def GetIEulerCoefs():
    ButchersTable = numpy.array([[1]])
    ButchersTau   = numpy.array([1])
    ButchersY     = numpy.array([1])
    return (ButchersTable, ButchersY, ButchersTau)

def GetIRungeKutta2Coefs():
    ButchersTable = numpy.array([[1/4, 3/4], 
                                 [-5/12, 3/4]])
    ButchersTau   = numpy.array([1, 1/3])
    ButchersY     = numpy.array([1/4, 3/4])
    return (ButchersTable, ButchersY, ButchersTau)

def GetIRungeKutta3Coefs():
    ButchersTable = numpy.array([[1/4, 3/4], 
                                 [-5/12, 3/4]])
    ButchersTau   = numpy.array([1, 1/3])
    ButchersY     = numpy.array([1/4, 3/4])
    return (ButchersTable, ButchersY, ButchersTau)

def GetIRungeKutta4Coefs():
    ButchersTable = numpy.array([[1/6, 2/6, 2/6, 1/6],
                                 [-2/6, 2/6, 2/6, 1/6],
                                 [1/6, -1/6, 2/6, 1/6],
                                 [1/6, 2/6, -4/6, 1/6]])
    ButchersTau   = numpy.array([1, 1/2, 1/2, 0])
    ButchersY     = numpy.array([1/6, 2/6, 2/6, 1/6])
    return (ButchersTable, ButchersY, ButchersTau)

def GetIGauss2Coefs():
    ButchersTable = numpy.array([[1/2]])
    ButchersTau   = numpy.array([1/2])
    ButchersY     = numpy.array([1])
    return (ButchersTable, ButchersY, ButchersTau)

def GetIHammerHollingsworth4Coefs():
    ButchersTable = numpy.array([[1/4, 1/4 - math.sqrt(3)/6],
                                 [1/4 + math.sqrt(3)/6, 1/4]])
    ButchersTau   = numpy.array([1/2 - math.sqrt(3)/6, 1/2 + math.sqrt(3)/6])
    ButchersY     = numpy.array([1/2, 1/2])
    return (ButchersTable, ButchersY, ButchersTau)

def GetIGauss6Coefs():
    ButchersTable = numpy.array([[5/36, 2/9 - math.sqrt(15)/15, 5/36 - math.sqrt(15)/30],
                                [5/36 + math.sqrt(15)/24, 2/9, 5/36-math.sqrt(15)/24],
                                 [5/36 + math.sqrt(15)/30, 2/9 + math.sqrt(15)/15, 5/36]])
    ButchersTau   = numpy.array([1/2 - math.sqrt(15)/10, 1/2, 1/2 + math.sqrt(15)/10])
    ButchersY     = numpy.array([5/18, 4/9, 5/18])
    return (ButchersTable, ButchersY, ButchersTau)

def GetIRado3Coefs():
    ButchersTable = numpy.array([[5/12, -1/12],
                                 [3/4, 1/4]])
    ButchersTau   = numpy.array([1/3, 1])
    ButchersY     = numpy.array([3/4, 1/4])
    return (ButchersTable, ButchersY, ButchersTau)

def ILobatto2Coefs():
    ButchersTable = numpy.array([[0, 0],
                                 [1/2, 12]])
    ButchersTau   = numpy.array([0, 1])
    ButchersY     = numpy.array([1/2, 1/2])
    return (ButchersTable, ButchersY, ButchersTau)

def ILobatto4Coefs():
    ButchersTable = numpy.array([[0, 0, 0],
                                 [5/24, 1/3, -1/24],
                                 [1/6, 2/3, 1/6]])
    ButchersTau   = numpy.array([0, 1/2, 1])
    ButchersY     = numpy.array([1/6, 2/3, 1/6])
    return (ButchersTable, ButchersY, ButchersTau)

def GetIGauss6Coefs():
    ButchersTable = numpy.array([[1/2]])
    ButchersTau   = numpy.array([1/2])
    ButchersY     = numpy.array([1])
    return (ButchersTable, ButchersY, ButchersTau)

def IEuler1(t0, y0, tau, stepsCount, funct, eps = 1e-6):
    if (isinstance(y0, Iterable)):
        return RungeKuttaImlicitSystem(t0, y0, tau, stepsCount, funct, *GetIEulerCoefs(), eps)
    else:
        return RungeKuttaImlicit(t0, y0, tau, stepsCount, funct, *GetIEulerCoefs(), eps)

def IRungeKutta2(t0, y0, tau, stepsCount, funct, eps = 1e-6):
    if (isinstance(y0, Iterable)):
        return RungeKuttaImlicitSystem(t0, y0, tau, stepsCount, funct, *GetIRungeKutta2Coefs(), eps)
    else:
        return RungeKuttaImlicit(t0, y0, tau, stepsCount, funct,*GetIRungeKutta2Coefs(), eps)

def IRungeKutta3(t0, y0, tau, stepsCount, funct, eps = 1e-6):
    if (isinstance(y0, Iterable)):
        return RungeKuttaImlicitSystem(t0, y0, tau, stepsCount, funct, *GetIRungeKutta3Coefs(), eps)
    else:
        return RungeKuttaImlicit(t0, y0, tau, stepsCount, funct, *GetIRungeKutta3Coefs(), eps)

def IRungeKutta4(t0, y0, tau, stepsCount, funct, eps = 1e-6):
    if (isinstance(y0, Iterable)):
        return RungeKuttaImlicitSystem(t0, y0, tau, stepsCount, funct, *GetIRungeKutta4Coefs(), eps)
    else:
        return RungeKuttaImlicit(t0, y0, tau, stepsCount, funct, *GetIRungeKutta4Coefs(), eps)
    
def Adams(t0  : [float], 
          y0  : [float], 
          tau : float, 
          stepsCount : int,
          funct, 
          AdamsCoefficients : numpy.array(float)):
    """
    Решает ОДУ или систему ОДУ методом Адамса.

    Аргументы:
        t0, y0            - начальное приближение (массив точек для ОДУ).
        tau               - шаг по времени.
        stepsCount        - число шагов по времени.
        funct             - функция, задающая уравнение, или вектор-функция для системы.
        AdamsCoefficients - коэффициенты метода Адамса (k штук).
                            u_{n+1} = u_{n} + tau * sum_{j=0}^{k-1} AdamsCoefficients_j * funct(t_n - tau * j)
        Пример:
            Метод Эйлера 1 порядка:
                AdamsCoefficients = [1]

            Метод Адамса 3 порядка:
                AdamsCoefficients = [23/12, -16/12, 5/12]

    Возвращает t, y - решение (либо массивы скаляров для ОДУ, либо массивы векторов для системы ОДУ).
    """
    assert(len(t0) == len(y0))
    assert(len(t0) == len(AdamsCoefficients))

    t = [ti for ti in t0]
    y = [yi for yi in y0]
    f = [funct(t0[st], y0[st]) for st in range(len(t0) - 1, -1, -1)]

    tn = t0[len(t0) - 1]
    yn = y0[len(y0) - 1]

    for n in range(0, stepsCount):
        sum = 0
        for st in range(0, len(AdamsCoefficients)):
            sum += AdamsCoefficients[st] * f[st]
        yn1 = yn + tau * sum
        tn1 = tn + tau

        t.append(tn1)
        y.append(yn1)

        f.pop()
        f.insert(0, funct(tn1, yn1))
        tn = tn1
        yn = yn1

    return numpy.array(t), numpy.array(y)

def Adams2(t0, y0, tau, stepsCount, funct):
    AdamsCoefficients = numpy.array([3/2, -1/2])
    t, y = Euler2(t0, y0, tau, len(AdamsCoefficients) - 1, funct)
    return Adams(t, y, tau, stepsCount - len(AdamsCoefficients) + 1, funct, AdamsCoefficients)

def Adams3(t0, y0, tau, stepsCount, funct):
    AdamsCoefficients = numpy.array([23/12, -16/12, 5/12])
    t, y = RungeKutta3(t0, y0, tau, len(AdamsCoefficients) - 1, funct)
    return Adams(t, y, tau, stepsCount - len(AdamsCoefficients) + 1, funct, AdamsCoefficients)

def Adams4(t0, y0, tau, stepsCount, funct):
    AdamsCoefficients = numpy.array([55/24, -59/24, 37/24, -9/24])
    t, y = RungeKutta4(t0, y0, tau, len(AdamsCoefficients) - 1, funct)
    return Adams(t, y, tau, stepsCount - len(AdamsCoefficients) + 1, funct, AdamsCoefficients)
def IAdams(t0  : [float], 
           y0  : [float], 
           tau : float, 
           stepsCount : int,
           funct, 
           AdamsCoefficients : numpy.array(float),
           eps = 1e-6):
    """
    Решает ОДУ или систему ОДУ методом Адамса.

    Аргументы:
        t0, y0            - начальное приближение (массив точек для ОДУ).
        tau               - шаг по времени.
        stepsCount        - число шагов по времени.
        funct             - функция, задающая уравнение, или вектор-функция для системы.
        AdamsCoefficients - коэффициенты неявного метода Адамса (k штук).
        Пример:
            Неявный метод Адамса 3 порядка:
                AdamsCoefficients = [5/12, 8/12, -1/12]
                u_{n+1} = u_n + tau * (5/12 f_{n+1} + 8/12 f_n - 1/12 f_{n-1})

    Возвращает t, y - решение (либо массивы скаляров для ОДУ, либо массивы векторов для системы ОДУ).
    """
    assert(len(t0) == len(y0))
    assert(len(t0) == len(AdamsCoefficients) - 1)

    t = [ti for ti in t0]
    y = [yi for yi in y0]
    f = [funct(t0[st], y0[st]) for st in range(len(t0) - 1, -1, -1)]

    tn = t0[len(t0) - 1]
    yn = y0[len(y0) - 1]

    debug = False

    if (debug):
        print(f"f0 = {f}\n")

    for n in range(0, stepsCount):
        tn1 = tn + tau

        def F(yn1):
            sum = AdamsCoefficients[0] * funct(tn1, yn1)
            for st in range(1, len(AdamsCoefficients)):
                if (debug):
                    print(f"st = {st}, f = {f[st - 1]} * coef = {AdamsCoefficients[st]}\n")
                sum += AdamsCoefficients[st] * f[st - 1]
            delta = yn1 - (yn + tau * sum)
            if (debug):
                print(f"yn1 = {yn1}\n")
                print(f"delta = {delta}\n")
            return delta
        
        def derF(yn1):
            J = Derivs.CalcJacobiMatrix(F, yn1, Derivs.DerivCentralO2, tau/10)
            if (debug):
                print(f"J = {J}\n")
            return J
        
        yn1 = SysNonlinEqSolvers.SolveSystemNewton2(F, derF, yn, eps, Norms.NormV1)

        t.append(tn1)
        y.append(yn1)

        f.pop()
        f.insert(0, funct(tn1, yn1))
        tn = tn1
        yn = yn1

    return numpy.array(t), numpy.array(y)

def GetIAdams2Coefs():
    return numpy.array([1/2, 1/2])

def GetIAdams3Coefs():
    return numpy.array([5/12, 8/12, -1/12])

def GetIAdams4Coefs():
    return numpy.array([9/24, 19/24, -5/24, 1/24])

def IAdams2(t0, y0, tau, stepsCount, funct):
    AdamsCoefficients = GetIAdams2Coefs()
    return IAdams([t0], [y0], tau, stepsCount - len(AdamsCoefficients) + 1, funct, AdamsCoefficients)

def IAdams3(t0, y0, tau, stepsCount, funct):
    AdamsCoefficients = GetIAdams3Coefs()
    t, y = IRungeKutta3(t0, y0, tau, len(AdamsCoefficients) - 2, funct)
    return IAdams(t, y, tau, stepsCount - len(AdamsCoefficients) + 2, funct, AdamsCoefficients)

def IAdams4(t0, y0, tau, stepsCount, funct):
    AdamsCoefficients = GetIAdams4Coefs()
    t, y = IRungeKutta4(t0, y0, tau, len(AdamsCoefficients) - 2, funct)
    return IAdams(t, y, tau, stepsCount - len(AdamsCoefficients) + 2, funct, AdamsCoefficients)

def BDF(t0  : [float],
        y0  : [float],
        tau : float,
        stepsCount : int,
        funct,
        BDFCoefficients : numpy.array(float)):
    """
    Backward differentiation formula.
    
    Решает ОДУ или систему ОДУ по явным формулам дифференцирования назад.

    Аргументы:
        t0, y0          - начальное приближение (массив точек для ОДУ).
        tau             - шаг по времени.
        stepsCount      - число шагов по времени.
        funct           - функция, задающая уравнение, или вектор-функция для системы.
        BDFCoefficients - коэффициенты формул дифференцирования назад (k штук).
                            sum_{j=0}^{k-1} BDFCoefficients_j * y_{n - tau * j} = tau * f_n
        Пример:
            ФДН 2 порядка:
                BDFCoefficients = [1/2, -1/2]

            ФДН 3 порядка (неустойчива):
                BDFCoefficients = [1/3, 1/2, -1, 1/6]

    Возвращает t, y - решение (либо массивы скаляров для ОДУ, либо массивы векторов для системы ОДУ).
    """
    assert(len(t0) == len(y0))
    assert(len(t0) == len(BDFCoefficients) - 1)

    _yn = [y0[st] for st in range(len(y0) - 1, -1, -1)]

    t = [ti for ti in t0]
    y = [yi for yi in y0]
    
    tn = t0[len(t0) - 1]
    yn = y0[len(y0) - 1]

    for n in range(0, stepsCount):
        sum = 0
        for j in range(1, len(BDFCoefficients)):
            sum += BDFCoefficients[j] * _yn[j - 1]
        yn1 = 1/BDFCoefficients[0] * (tau * funct(tn, yn) - sum)
        tn1 = tn + tau

        t.append(tn1)
        y.append(yn1)

        _yn.pop()
        _yn.insert(0, yn1)
        tn = tn1
        yn = yn1

    return numpy.array(t), numpy.array(y)
    
def BDF2(t0, y0, tau, stepsCount, funct):
    BDFCoefficients = numpy.array([1/2, 0, -1/2])
    t, y = Euler2(t0, y0, tau, len(BDFCoefficients) - 2, funct)
    return BDF(t, y, tau, stepsCount - len(BDFCoefficients) + 2, funct, BDFCoefficients)

def BDF3(t0, y0, tau, stepsCount, funct):
    BDFCoefficients = numpy.array([1/3, 1/2, -1, 1/6])
    t, y = RungeKutta3(t0, y0, tau, len(BDFCoefficients) - 2, funct)
    return BDF(t, y, tau, stepsCount - len(BDFCoefficients) + 2, funct, BDFCoefficients)

# More methods https://en.wikipedia.org/w/index.php?title=Backward_differentiation_formula&oldid=1166087118

def IBDF(t0  : [float],
        y0  : [float],
        tau : float,
        stepsCount : int,
        funct,
        BDFCoefficients : numpy.array(float),
        eps = 1e-6):
    """
    Backward differentiation formula.
    
    Решает ОДУ или систему ОДУ по явным формулам дифференцирования назад.

    Аргументы:
        t0, y0          - начальное приближение (массив точек для ОДУ).
        tau             - шаг по времени.
        stepsCount      - число шагов по времени.
        funct           - функция, задающая уравнение, или вектор-функция для системы.
        BDFCoefficients - коэффициенты формул дифференцирования назад (k штук).
        Пример:
            Неявная ФДН 3 порядка:
                BDFCoefficients = [11/6, -3, 3/2, -1/3]

    Возвращает t, y - решение (либо массивы скаляров для ОДУ, либо массивы векторов для системы ОДУ).
    """
    assert(len(t0) == len(y0))
    assert(len(t0) == len(BDFCoefficients) - 1)

    _yn = [y0[st] for st in range(len(y0) - 1, -1, -1)]

    t = [ti for ti in t0]
    y = [yi for yi in y0]
    
    tn = t0[len(t0) - 1]
    yn = y0[len(y0) - 1]

    for n in range(0, stepsCount):
        tn1 = tn + tau

        def F(yn1):
            sum = 0
            for j in range(1, len(BDFCoefficients)):
                sum += BDFCoefficients[j] * _yn[j - 1]
            new_yn1 = 1/BDFCoefficients[0] * (tau * funct(tn1, yn1) - sum)
            return yn1 - new_yn1
        
        def derF(yn1):
            return Derivs.CalcJacobiMatrix(F, yn1, Derivs.DerivCentralO2, tau/10)
        
        yn1 = SysNonlinEqSolvers.SolveSystemNewton2(F, derF, yn, eps, Norms.NormV1)

        t.append(tn1)
        y.append(yn1)

        _yn.pop()
        _yn.insert(0, yn1)
        tn = tn1
        yn = yn1

    return numpy.array(t), numpy.array(y)

def GetIBDF2Coefs():
    return numpy.array([3/2, -2, 1/2])

def GetIBDF3Coefs():
    return numpy.array([11/6, -3, 3/2, -1/3])

def GetIBDF4Coefs():
    return numpy.array([25/12, -4, 3, -4/3, 1/4])

def GetIBDF5Coefs():
    return numpy.array([137/60, -5, 5, -10/3, 5/4, -1/5])

def GetIBDF6Coefs():
    return numpy.array([147/60, -6, 15/2, -20/3, 15/4, -6/5, 1/6])

def IBDF2(t0, y0, tau, stepsCount, funct):
    BDFCoefficients = GetIBDF2Coefs()
    t, y = IRungeKutta2(t0, y0, tau, len(BDFCoefficients) - 2, funct)
    return IBDF(t, y, tau, stepsCount - len(BDFCoefficients) + 2, funct, BDFCoefficients)

def IBDF3(t0, y0, tau, stepsCount, funct):
    BDFCoefficients = GetIBDF3Coefs()
    t, y = IRungeKutta3(t0, y0, tau, len(BDFCoefficients) - 2, funct)
    return IBDF(t, y, tau, stepsCount - len(BDFCoefficients) + 2, funct, BDFCoefficients)

def IBDF4(t0, y0, tau, stepsCount, funct):
    BDFCoefficients = GetIBDF4Coefs()
    t, y = IRungeKutta4(t0, y0, tau, len(BDFCoefficients) - 2, funct)
    return IBDF(t, y, tau, stepsCount - len(BDFCoefficients) + 2, funct, BDFCoefficients)

def IBDF5(t0, y0, tau, stepsCount, funct):
    BDFCoefficients = GetIBDF5Coefs()
    t, y = IRungeKutta4(t0, y0, tau, len(BDFCoefficients) - 2, funct)
    return IBDF(t, y, tau, stepsCount - len(BDFCoefficients) + 2, funct, BDFCoefficients)

def IBDF6(t0, y0, tau, stepsCount, funct):
    BDFCoefficients = GetIBDF6Coefs()
    t, y = IRungeKutta4(t0, y0, tau, len(BDFCoefficients) - 2, funct)
    return IBDF(t, y, tau, stepsCount - len(BDFCoefficients) + 2, funct, BDFCoefficients)

def CROS2_s1(t0  : [float],
             y0  : [float],
             tau : float,
             stepsCount : int,
             funct,
             eps = 1e-6):
    """
    y_{n+1} = y_n + tau * Re{k}
    (E - tau * (1+i)/2 * df(y_n)/dy) * k = f(y_n)
    """
    assert(not isinstance(y0, Iterable))

    t = [t0]
    y = [y0]

    tn = t0
    yn = y0

    for n in range(0, stepsCount):
        def funct_helper(yn):
            return funct(tn, yn)
        
        J = Derivs.DerivCentralO2(funct_helper, yn, tau/10, numpy.complex128)

        def F(k):
            return (1 - tau * (1+1j)/2 * J) * k - funct(tn, yn)
            
        def derF(k):
            return Derivs.CalcJacobiMatrix(F, k, Derivs.DerivCentralO2, tau/10, numpy.complex128)

        k = SysNonlinEqSolvers.SolveSystemNewton2(F, derF, complex(0), eps, Norms.NormV1)

        yn1 = yn + tau * k.real
        tn1 = tn + tau

        t.append(tn1)
        y.append(yn1)

        tn = tn1
        yn = yn1

    return numpy.array(t), numpy.array(y)

def SysCROS2_s1(t0  : [float],
                y0  : [float],
                tau : float,
                stepsCount : int,
                funct,
                eps = 1e-6):
    """
    y_{n+1} = y_n + tau * Re{k}
    (E - tau * (1+i)/2 * df(y_n)/dy) * k = f(y_n)
    """
    assert(isinstance(y0, Iterable))

    t = [t0]
    y = [y0]

    tn = t0
    yn = y0

    dim = len(y0)

    for n in range(0, stepsCount):
        def funct_helper(yn):
            return funct(tn, yn)

        J = Derivs.CalcJacobiMatrix(funct_helper, yn, Derivs.DerivCentralO2, tau/10, numpy.complex128)

        def F(k):
            E = numpy.eye(dim)
            return (E - tau * (1+1j)/2 * J) @ k - funct(tn, yn)
        
        def derF(k):
            return Derivs.CalcJacobiMatrix(F, k, Derivs.DerivCentralO2, tau/10, numpy.complex128)

        k = SysNonlinEqSolvers.SolveSystemNewton2(F, derF, numpy.zeros(dim, dtype=complex), eps, Norms.NormV1)

        yn1 = yn + tau * k.real
        tn1 = tn + tau

        t.append(tn1)
        y.append(yn1)

        tn = tn1
        yn = yn1

    return numpy.array(t), numpy.array(y)

def BinomCoef(n, k):
    C_n_k = 1
    for st in range(1, k + 1):
        C_n_k *= (n - st + 1) / (st)
    return C_n_k

def GetE1(k1):
    e1 = numpy.zeros((1, k1))
    e1[0][1] = 1
    return e1

def GetP(k1):
    P = numpy.zeros((k1, k1))
    for j in range(0, k1):
        for i in range(0, j + 1):
            P[i][j] = BinomCoef(j, i)
    return P

def SolveNordsieckEquation(t0  : [float],
                           y0  : [float],
                           tau : float,
                           stepsCount : int,
                           funct,
                           l,
                           eps = 1e-6):
    zn = numpy.zeros(len(l))
    zn[0] = y0
    zn[1] = tau * funct(t0, y0)

    def derY_T(t, y, k, tau):
        if (k == 1):
            return funct(t, y)
        
        derF_t = (derY_T(t + tau, y, k - 1, tau) - derY_T(t - tau, y, k - 1, tau)) / (2*tau)
        derF_y = (derY_T(t, y + tau, k - 1, tau) - derY_T(t, y - tau, k - 1, tau)) / (2*tau)
        derY_t = funct(t, y)
        return derF_t + derF_y * derY_t
    
    mn = tau
    for st in range(2, len(l)):
        mn *= tau / st
        zn[st] = mn * derY_T(t0, y0, st, tau/10)

    tn = t0
    yn = y0

    t = [t0]
    y = [y0]

    P = GetP(len(l))
    e1 = GetE1(len(l))

    for n in range(0, stepsCount):
        tn1 = tn + tau

        zn_corr = P @ zn - l[:, numpy.newaxis] @ e1 @ P @ zn

        def F(yn1):
            zn1 = l * tau * funct(tn1, yn1) + zn_corr
            return zn1[0] - yn1
        
        def derF(yn1):
            return Derivs.CalcJacobiMatrix(F, yn1, Derivs.DerivCentralO2, tau/10)
        
        yn1 = SysNonlinEqSolvers.SolveSystemNewton2(F, derF, yn, eps, Norms.NormV1)

        t.append(tn1)
        y.append(yn1)

        zn = l * tau * funct(tn1, yn1) + zn_corr
        tn = tn1
        yn = yn1
    
    return numpy.array(t), numpy.array(y)

def TensorDot(A, B):
    res = numpy.zeros((len(A) * len(B), len(A[0]) * len(B[0])))
    for aRowIndex in range(0, len(A)):
        for aColumnIndex in range(0, len(A[aRowIndex])):
            subMatr = A[aRowIndex][aColumnIndex] * B
            for st1 in range(0, len(subMatr)):
                for st2 in range(0, len(subMatr[st1])):
                    res[aRowIndex * len(B) + st1][aColumnIndex * len(B) + st2] = subMatr[st1][st2]
    return res

def SolveNordsieckSystem(t0  : [float],
                         y0  : [float],
                         tau : float,
                         stepsCount : int,
                         funct,
                         l,
                         eps = 1e-7):
    zn = numpy.zeros((len(l), len(y0)))
    zn[0] = y0
    zn[1] = tau * funct(t0, y0)

    def derY_T(t, y, k, tau):
        if (k == 1):
            return funct(t, y)
        
        derF_t = (derY_T(t + tau, y, k - 1, tau) - derY_T(t - tau, y, k - 1, tau)) / (2*tau)
        derF_y = 0
        derY_t = funct(t, y)

        for st in range(0, len(y)):
            yp = 1.0*y
            ym = 1.0*y
            yp[st] = yp[st] + tau
            ym[st] = ym[st] - tau
            derF_y += (derY_T(t, yp, k - 1, tau) - derY_T(t, ym, k - 1, tau)) / (2*tau) * derY_t[st]

        return derF_t + derF_y
    
    mn = tau
    for st in range(2, len(l)):
        mn *= tau / st
        zn[st] = mn * derY_T(t0, y0, st, tau)

    tn = t0
    yn = y0

    t = [t0]
    y = [y0]
    
    P = GetP(len(l))
    E = numpy.eye(len(y0))
    e1 = GetE1(len(l))

    size = len(l) * len(y0)
    zn = zn.reshape((size, 1))

    for n in range(0, stepsCount):
        tn1 = tn + tau

        zn_corr1 = TensorDot(P, E) @ zn
        zn_corr2 = -TensorDot(l[:, numpy.newaxis], E) @ TensorDot(e1 @ P, E) @ zn

        def F(yn1):
            fn1 = funct(tn1, yn1).reshape((len(y0), 1))
            zn1 = tau * TensorDot(l[:, numpy.newaxis], E) @ fn1 + zn_corr1 + zn_corr2
            delta = numpy.zeros(len(yn1))
            for st in range(0, len(delta)):
                delta[st] = zn1[st][0] - yn1[st]
            return delta
        
        def derF(yn1):
            return Derivs.CalcJacobiMatrix(F, yn1, Derivs.DerivCentralO2, tau/10)
        
        yn1 = SysNonlinEqSolvers.SolveSystemNewton2(F, derF, yn, eps, Norms.NormV1, maxIters = 100)

        t.append(tn1)
        y.append(yn1)

        fn1 = funct(tn1, yn1).reshape((len(y0), 1))
        zn = tau * TensorDot(l[:, numpy.newaxis], E) @ fn1 + zn_corr1 + zn_corr2
        tn = tn1
        yn = yn1
    
    return numpy.array(t), numpy.array(y)
def GetNordsieckAdams1L():
    return numpy.array([1/2, 1])

def GetNordsieckAdams2L():
    return numpy.array([5/12, 1, 1/2])

def GetNordsieckAdams3L():
    return numpy.array([3/8, 1, 3/4, 1/6])

def GetNordsieckAdams4L():
    return numpy.array([251/720, 1, 11/12, 1/3, 1/24])

def GetNordsieckAdams5L():
    return numpy.array([95/288, 1, 25/24, 35/72, 5/48, 1/120])

def GetNordsieckAdams6L():
    return numpy.array([19087/60480, 1, 137/120, 5/8, 17/96, 1/40, 1/720])

def NordsieckAdams1(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams1L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckAdams2(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams2L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckAdams3(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams3L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckAdams4(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams4L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckAdams5(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams5L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckAdams6(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams6L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemAdams1(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams1L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemAdams2(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams2L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemAdams3(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams3L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemAdams4(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams4L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemAdams5(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams5L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemAdams6(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckAdams6L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def GetNordsieckBDF1L():
    return numpy.array([1, 1])

def GetNordsieckBDF2L():
    return numpy.array([2/3, 1, 1/3])

def GetNordsieckBDF3L():
    return numpy.array([6/11, 1, 6/11, 1/11])

def GetNordsieckBDF4L():
    return numpy.array([12/25, 1, 7/10, 1/5, 1/50])

def GetNordsieckBDF5L():
    return numpy.array([60/137, 1, 225/274, 85/274, 15/274, 1/274])

def GetNordsieckBDF6L():
    return numpy.array([20/49, 1, 58/63, 5/12, 25/252, 1/84, 1/1764])

def NordsieckBDF1(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF1L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckBDF2(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF2L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckBDF3(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF3L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckBDF4(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF4L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckBDF5(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF5L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckBDF6(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF6L()
    return SolveNordsieckEquation(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemBDF1(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF1L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemBDF2(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF2L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemBDF3(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF3L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemBDF4(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF4L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemBDF5(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF5L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)

def NordsieckSystemBDF6(t0  : [float],
                    y0  : [float],
                    tau : float,
                    stepsCount : int,
                    funct):
    l = GetNordsieckBDF6L()
    return SolveNordsieckSystem(t0, y0, tau, stepsCount, funct, l)