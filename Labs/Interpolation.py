import numpy
import Gauss

class Spline:
    def __init__(self, x, f):
        def FunctEq(A, b, i, rowIndex, x_i, f_i):
            columnIndex = 4 * i
            # a_i
            A[rowIndex[0]][columnIndex]     = x_i**3
            # b_i
            A[rowIndex[0]][columnIndex + 1] = x_i**2
            # c_i
            A[rowIndex[0]][columnIndex + 2] = x_i
            # d_i
            A[rowIndex[0]][columnIndex + 3] = 1

            b[rowIndex[0]] = f_i
            rowIndex[0] += 1
        
        def FunctDerEq(A, b, i, rowIndex, x_i):
            columnIndex = 4 * i
            # a_i
            A[rowIndex[0]][columnIndex]     = 3 * x_i**2
            # b_i
            A[rowIndex[0]][columnIndex + 1] = 2 * x_i
            # c_i
            A[rowIndex[0]][columnIndex + 2] = 1

            # a_i+1
            A[rowIndex[0]][columnIndex + 4] = -3 * x_i**2
            # b_i+1
            A[rowIndex[0]][columnIndex + 5] = -2 * x_i
            # c_i+1
            A[rowIndex[0]][columnIndex + 6] = -1

            b[rowIndex[0]] = 0
            rowIndex[0] += 1

        def FunctDer2Eq(A, b, i, rowIndex, x_i):
            columnIndex = 4 * i
            # a_i
            A[rowIndex[0]][columnIndex]     = 6 * x_i
            # b_i
            A[rowIndex[0]][columnIndex + 1] = 2

            # a_i+1
            A[rowIndex[0]][columnIndex + 4] = -6 * x_i
            # b_i+1
            A[rowIndex[0]][columnIndex + 5] = -2

            b[rowIndex[0]] = 0
            rowIndex[0] += 1

        def FunctDer2Eq0(A, b, i, rowIndex, x_i):
            columnIndex = 4 * i
            # a_i+1
            A[rowIndex[0]][columnIndex]     = 6 * x_i
            # b_i+1
            A[rowIndex[0]][columnIndex + 1] = 2

            b[rowIndex[0]] = 0
            rowIndex[0] += 1

        # P_3(x) = ax^3 + bx^2 + cx + d.
        # a1, b1, c1, d1, a2, b2, c2, d2, ..., an, bn, cn, dn.
        n = len(x) - 1
        A = numpy.zeros((4*n, 4*n))
        b = numpy.empty(4 * n)

        rowIndex = [0]
        # Равенство значений кубических интерполянтов в точках
        FunctEq(A, b, 0, rowIndex, x[0], f[0])
        for i in range(1, len(x) - 1):
            FunctEq(A, b, i-1, rowIndex, x[i], f[i])
            FunctEq(A, b, i, rowIndex, x[i], f[i])
        FunctEq(A, b, n - 1, rowIndex, x[len(x) - 1], f[len(x) - 1])

        # Равенство первых производных
        for i in range(1, len(x) - 1):
            FunctDerEq(A, b, i-1, rowIndex, x[i])
        
        # Равенство вторых производных
        for i in range(1, len(x) - 1):
            FunctDer2Eq(A, b, i-1, rowIndex, x[i])
        
        FunctDer2Eq0(A, b, 0, rowIndex, x[0])
        FunctDer2Eq0(A, b, n - 1, rowIndex, x[len(x) - 1])

        self.coefs = Gauss.SolveGauss(A, b)
        self.x = x

    def Interpolate(self, x):
        def Calc(coefs, i, x):
            startI = 4 * i
            return coefs[startI] * x**3 + coefs[startI + 1] * x**2 + coefs[startI + 2] * x + coefs[startI + 3]
        
        if (x <= self.x[0]):
            return Calc(self.coefs, 0, x)
        elif (x >= self.x[len(self.x) - 2]):
            return Calc(self.coefs, len(self.x) - 2, x)
        else:
            for st in range(0, len(self.x) - 1):
                if (self.x[st] <= x and x <= self.x[st + 1]):
                    return Calc(self.coefs, st, x)
                
class MLS:
    def __init__(self, x, f):
        A = numpy.empty((len(x), 4))
        _f = f.copy()
        for rowIndex in range(0, len(x)):
            A[rowIndex][0] = x[rowIndex]**3
            A[rowIndex][1] = x[rowIndex]**2
            A[rowIndex][2] = x[rowIndex]
            A[rowIndex][3] = 1
        _A = numpy.dot(numpy.transpose(A), A)
        _f = numpy.dot(numpy.transpose(A), _f)
        self.coefs = Gauss.SolveGauss(_A, _f)

    def Interpolate(self, x):
        return self.coefs[0] * x**3 + self.coefs[1] * x**2 + self.coefs[2] * x + self.coefs[3]
    
class LinearInterp:
    # ax + b
    class LinInt:
        def __init__(self, x1, x2, f1, f2):
            if (x1 < x2):
                self.x1 = x1
                self.x2 = x2
            else:
                self.x1 = x2
                self.x2 = x1

            self.a  = (f2 - f1) / (x2 - x1)
            self.b  = f2 - self.a * x2
            
        def ContainsX(self, x):
            return self.x1 <= x and x <= self.x2
            
        def Interpolate(self, x):
            return self.a * x + self.b

    def Compute(x, f):
        assert(len(x) == len(f))
        assert(len(x) >= 2)
        linInts = []
        
        x1 = x[0]
        for st in range(1, len(x)):
            linInt = LinearInterp.LinInt(x[st - 1], x[st], f[st - 1], f[st])
            linInts.append(linInt)
        
        return linInts

    def __init__(self, x, f):
        self.x = x
        self.f = f
        self.coefs = LinearInterp.Compute(x, f)
       
    def GetX(self, index):
        return self.x[index]
        
    def GetF(self, index):
        return self.f[index]
    
    def Interpolate1(self, x):
        return self.Interpolate2(x, 0, len(self.coefs))
    
    def Interpolate2(self, x, startIndex, stopIndex):
        first = self.coefs[0]
        last  = self.coefs[len(self.coefs) - 1]
        
        if (x <= first.x1):
            return 0 # The pressure equal to 'p_inf' before shock wave appears.
        elif (x >= last.x2):
            return last.Interpolate(last.x2)
        else:
            for st in range(startIndex, stopIndex):
                if (self.coefs[st].ContainsX(x)):
                    return self.coefs[st].Interpolate(x)
            
            print(f"LinearInterp logic error. Missed interval:\n"
                  f"    startIndex = {startIndex}\n"
                  f"    stopIndex  = {stopIndex}\n"
                  f"    x          = {x}\n"
                  f"    startX1    = {self.coefs[startIndex].x1}\n"
                  f"    startX2    = {self.coefs[startIndex].x2}\n"
                  f"    stopX1     = {self.coefs[stopIndex].x1}\n"
                  f"    stopX2     = {self.coefs[stopIndex].x2}\n")

            raise Exception("LinearInterp logic error.")
            
    # Solves equation "ax + b = Interpolate(x)".
    def SolveEquation(self, a, b):
        roots = []
        for coef in self.coefs:
            if (coef.a == a):
                continue
            x_inter = - (b - coef.b) / (a - coef.a)
            if (coef.x1 <= x_inter and x_inter <= coef.x2):
                roots.append(x_inter)
        return roots
