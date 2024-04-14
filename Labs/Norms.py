import numpy
import math

def NormV1(x):
    return numpy.max(numpy.abs(x))

def NormV2(x):
    return numpy.sum(numpy.abs(x))

def NormV3(x):
    return math.sqrt(numpy.dot(x, x))

def NormM1(A):
    # Максимальная сумма по строке модулей элементов.
    result = 0
    for rowIndex in range(0, len(A)):
        sum = 0
        for columnIndex in range(0, len(A)):
            sum += abs(A[rowIndex][columnIndex])
        if (sum > result):
            result = sum
    return result

def NormM2(A):
    # Максимальная сумма по столбцам модулей элементов.
    result = 0
    for columnIndex in range(0, len(A)):
        sum = 0
        for rowIndex in range(0, len(A)):
            sum += abs(A[rowIndex][columnIndex])
        if (sum > result):
            result = sum
    return result