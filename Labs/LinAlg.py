import numpy
import Gauss

# Находит ошибку решения системы уравнений: residual = Af - s.
def CheckSystemSolution(A, f, s):
    residual = numpy.empty(len(A))
    for rowIndex in range(0, len(A)):
        residual[rowIndex] = numpy.dot(A[rowIndex], s) - f[rowIndex]
    return residual

# Вычисление детерминанта по формуле разложения по строке.
# Медленно работает из-за постоянного выделения памяти для миноров.
def MatrixDeterminantVeryVerySlow(A):
    if (len(A) == 1):
        return A[0][0]
    if (len(A) == 2):
        return A[0][0] * A[1][1] - A[1][0] * A[0][1]
    
    det = 0
    sign = 1
    for columnIndex in range(0, len(A)):
        cols = numpy.array(list(range(0, columnIndex)) + list(range(columnIndex + 1, len(A))))
        rows = numpy.arange(1, len(A))
        minor = A[cols[:,numpy.newaxis], rows]
        det += sign * A[0][columnIndex] * MatrixDeterminantVeryVerySlow(minor)
        sign = -sign
    return det

# Методом Гаусса приводим к верхне-диагональному виду, затем находим определить как произведение диагональных элементов верхне-треугольной матрицы.
def MatrixDeterminant1(A, findMaxElem):
    det = 1
    # От первой до предпоследней строки. Последняя строка уже приведена к нужному виду.
    for columnIndex in range(0, len(A) - 1):
        rowIndex = columnIndex
        
        # Ищем максимальный элемент.
        # Строки матрицы можно переставлять, значение определителя поменяет знак.
        if (findMaxElem):
            maxRowIndex = rowIndex
            maxElem = A[rowIndex][columnIndex]
            for _rowIndex in range(rowIndex + 1, len(A)):
                if (A[_rowIndex][columnIndex] > maxElem):
                    maxElem = A[_rowIndex][columnIndex]
                    maxRowIndex = _rowIndex
            if (rowIndex != maxRowIndex):
                # Меняем строки местами.
                A[[rowIndex, maxRowIndex]] = A[[maxRowIndex, rowIndex]]
                det = -det

        # Прямой ход метода Гаусса.
        a = A[rowIndex] / A[rowIndex][columnIndex]
        # Обнуляем элементы под диагональным.
        for _rowIndex in range(rowIndex + 1, len(A)):
            A[_rowIndex] -= a * A[_rowIndex][columnIndex]
            # Присваиваем точное значение, из-за возможной ошибки при округлении.
            A[_rowIndex][columnIndex] = 0
        # Детерминант равен произведению диагональных элементов.
        det *= A[rowIndex][columnIndex]
    det *= A[len(A) - 1, len(A) - 1]
    return det

# Вычисление детерминанта по формуле полного разложения (через перестановки) слишком вычислительно затратный.
# Для формулы полного разложения 100! ~9,332621544×10^157 - слишком большое число. 
# В тоже время метод Гаусса требует ~(2/3 * 100^3 + 100^2) итераций - что приемлемо.

# Нахождение обратной матрицы методом Гаусса.
# Не изменяет исходную матрицу A.
def MatrixInverse1(A):
    invA = numpy.identity(len(A))
    srcA = numpy.empty((len(A), len(A)))
    numpy.copyto(srcA, A)

    Gauss.MakeUpDiagonal(srcA, invA, False)
    Gauss.MakeDiagonal(srcA, invA)
    return invA

def MatrixInverse3(A):
    return numpy.linalg.inv(A)