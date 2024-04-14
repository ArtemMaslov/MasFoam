###################################################################################################################
# Модуль численного интегрирования.
#
# Версия: 1.
# Дата изменения: 01.03.2024, MasFoam.
###################################################################################################################

from Interpolation import LinearInterp

def IntegrateO1Right(x, f):
    result = 0
    for st in range(1, len(f)):
        result += (x[st] - x[st - 1]) * f[st]
    return result

def IntegrateO1Right2(x1, x2, dx, funct):
    result = 0
    x = x1
    while (x < x2):
        result += funct(x + dx) * dx
        x += dx
    return result

def IntegrateO1Left(x, f):
    result = 0
    for st in range(1, len(f)):
        result += f[st - 1] * (x[st] - x[st - 1])
    return result

def IntegrateO1Left2(x1, x2, dx, funct):
    result = 0
    x = x1
    while (x < x2):
        result += funct(x) * dx
        x += dx
    return result

def IntegrateO1Central(x, f):
    result = 0
    for st in range(1, len(f)):
        result += (x[st] - x[st - 1]) * 1/2 * (f[st - 1] + f[st])
    return result

def IntegrateO1Left2(x1, x2, dx, funct):
    result = 0
    x = x1
    while (x < x2):
        result += funct(x + dx/2) * dx
        x += dx
    return result

def IntegrateTrapezoid(x, f):
    step = x[1] - x[0]
    integral = 1/2 * (f[0] + f[len(f) - 1])
    for st in range(1, len(f) - 1):
        integral += f[st]
    return step * integral
    
def IntegrateTrapezoid2(x1, x2, dx, funct):
    result = 1/2 * (funct(x1) + funct(x2)) * dx
    x = x1 + dx
    while (x < x2 - dx):
        result += funct(x) * dx
        x += dx
    return result

def IntegrateTrapInterp(interp : LinearInterp, x1, x2, step, startIndex, stopIndex):
    integral = 1/2 * (interp.Interpolate2(x1, startIndex, stopIndex) + interp.Interpolate2(x2, startIndex, stopIndex))
    x = x1 + step
    while (x < x2):
        integral += interp.Interpolate2(x, startIndex, stopIndex)
        x += step
        #print(f"x = {x} < x2 = {x2}")
    return integral * step

def IntegrateSimpson(x, f):
    step = x[1] - x[0]
    integral = f[0] + f[len(f) - 1]
    for st in range(1, len(f) - 1):
        if (st % 2 == 1):
            integral += 4 * f[st]
        else:
            integral += 2 * f[st]
    return 1/3 * step * integral

def IntegrateSimpson2(x1, x2, dx, funct):
    result = (funct(x1) + funct(x2)) * dx
    x = x1 + dx
    st = 1
    while (x < x2 - dx):
        if (st % 2 == 1):
            result += 4/3 * funct(x) * dx
        else:
            result += 2/3 * funct(x) * dx
        st += 1
        x += dx
    return result

# p - порядок функции integrationMethod
def IntegrateRichardson(x, f, integrationMethod, p):
    assert (len(x) % 2 == 1)

    x2 = []
    f2 = []
    for st in range(0, len(x), 2):
        x2.append(x[st])
        f2.append(f[st])

    Ih  = integrationMethod(x, f)
    I2h = integrationMethod(x2, f2)
    
    return Ih + (Ih - I2h) / (2**p - 1)

###################################################################################################################
###################################################################################################################