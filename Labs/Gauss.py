def MakeDiagonal(A, f):
    for rowIndex in range(len(A)-1, -1, -1):
        columnIndex = rowIndex

        # Уже сделали это в MakeUpDiagonal.
        #f[rowIndex] /= A[rowIndex][columnIndex]
        #A[rowIndex][columnIndex] = 1

        for _rowIndex in range(0, rowIndex):
            f[_rowIndex] -= f[rowIndex] * A[_rowIndex][columnIndex]
            A[_rowIndex][columnIndex] = 0

def MakeUpDiagonal(A, f, findMaxElem):
    for columnIndex in range(0, len(A)):
        rowIndex = columnIndex
        
        # Ищем максимальный элемент.
        if (findMaxElem):
            maxRowIndex = rowIndex
            maxElem = A[rowIndex][columnIndex]
            for _rowIndex in range(rowIndex + 1, len(A)):
                # Если главный элемент равен 0, его нужно заменить.
                if (A[_rowIndex][columnIndex] != 0 and
                    (A[_rowIndex][columnIndex] > maxElem or
                     maxElem == 0)):
                    maxElem = A[_rowIndex][columnIndex]
                    maxRowIndex = _rowIndex
            if (rowIndex != maxRowIndex):
                # Меняем строки местами.
                A[[rowIndex, maxRowIndex]] = A[[maxRowIndex, rowIndex]]
                f[[rowIndex, maxRowIndex]] = f[[maxRowIndex, rowIndex]]

        # Прямой ход метода Гаусса.

        # Делаем диагональный элемент равным 1.
        a = A[rowIndex][columnIndex]
        A[rowIndex] /= a
        f[rowIndex] /= a
        # Присваиваем точное значение, из-за возможной ошибки при округлении.
        A[rowIndex][columnIndex] = 1
        # Обнуляем элементы под диагональным.
        for _rowIndex in range(rowIndex + 1, len(A)):
            f[_rowIndex] -= f[rowIndex] * A[_rowIndex][columnIndex]
            A[_rowIndex] -= A[rowIndex] * A[_rowIndex][columnIndex]
            # Присваиваем точное значение, из-за возможной ошибки при округлении.
            A[_rowIndex][columnIndex] = 0
 
# Решает систему линейных уравнений Ax = f методом Гаусса. 
# Внимание! Изменяет исходные значения A и f.
def SolveGauss(A, f, findMaxElem = True):
    MakeUpDiagonal(A, f, findMaxElem)
    MakeDiagonal(A, f)
    return f