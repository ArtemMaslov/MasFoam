#pragma once

#include <MasFoam/Data/Vector.h>
#include <MasFoam/Data/Matrix.h>

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

namespace MasFoam
{
    template <typename T>
    void ApplyForwardGaussMethod(Matrix<T>& A, Vector<T>& f, bool findMaxElem = true)
    {
        const size_t rowsCount    = A.GetRowsCount();
        const size_t columnsCount = A.GetColumnsCount();

        for (size_t rowIndex = 0; rowIndex < rowsCount; rowIndex++)
        {
            T maxElem = 0;
            if (findMaxElem)
            {
                // Ищем максимальный элемент ниже текущего диагонального.
                size_t maxElemRowIndex = 0;
                maxElem = A.GetElement(rowIndex);
                for (size_t findMaxIndex = rowIndex + 1; findMaxIndex < rowsCount; findMaxIndex++)
                {
                    T currentElem = A.GetElement(findMaxIndex, rowIndex);
                    if (currentElem > maxElem)
                    {
                        maxElem = currentElem;
                        maxElemRowIndex = findMaxIndex;
                    }
                }
                // Меняем строки местами.
                A.SwapRows(rowIndex, maxElemRowIndex);
                std::swap(f[rowIndex], f[maxElemRowIndex]);
            }
            else
                maxElem = A.GetElement(rowIndex, columnsCount);
            
            // Делим текущую строку на значения диагонального элемента. Диагональный элемент становится равным 1.
            for (size_t columnIndex = rowIndex + 1; columnIndex < columnsCount; columnIndex++)
                A.GetElement(rowIndex, columnIndex) /= maxElem;
            f[rowIndex] /= maxElem;
            // Делаем диагональный элемент равным 1.
            A.GetElement(rowIndex, rowIndex) = 1;
            
            // Обнуляем элементы под текущим диагональным.
            for (size_t nullRowIndex = rowIndex + 1; nullRowIndex < rowsCount; nullRowIndex++)
            {
                T rowNullKoef = A.GetElement(nullRowIndex, rowIndex);
                for (size_t columnIndex = rowIndex + 1; columnIndex < columnsCount; columnIndex++)
                    A.GetElement(nullRowIndex, columnIndex) -= A.GetElement(rowIndex, columnIndex) * rowNullKoef;
                A.GetElement(nullRowIndex, rowIndex) = 0;
                f[nullRowIndex] -= f[rowIndex] * rowNullKoef;
            }
        }
    }

    /**
     * @brief  Применяет обратный ход метода Гаусса к системе линейных уравнений.
     * Система: Ax = f, где 
     *     A - квадратная матрица системы.
     *     f - вектор свободных коэффициентов.
     *     x - вектор искомых переменных.
     * 
     * Внимание! Матрица системы A должна быть квадратной.
     * 
     * Внимание! Матрица системы A должна быть в верхне-треугольном виде с произвольными ненулевыми диагональными 
     * элементами.
     * 
     * В ходе работы метода для оптимизации матрица системы не изменяется. Вектор свободных коэффициентов в конце
     * работы алгоритма будет равен вектору решений системы линейных уравнений.
     * 
     * @tparam T Тип используемых чисел. Обычно либо double, либо float. Целочисленные типы имеют смысл, если нет 
     *           деления.
     * 
     * @param A Матрица системы.
     * @param f Вектор свободных коэффициентов.
    */
    template <typename T>
    void ApplyBackwardGaussMethod(const Matrix<T>& A, Vector<T>& f)
    {
        const size_t rowsCount = A.GetRowsCount();

        // Обрабатываем все диагональные элементы.
        for (ptrdiff_t rowIndex = rowsCount - 1; rowIndex > 0; rowIndex--)
        {
            // На каждом шаге элемент на диагонали [rowIndex, rowIndex] возможно не равен 1,
            // Все элементы в строке правее диагонального равны 0.

            // Делаем диагональный элемент равным 1. f[rowIndex] будет равен корню №rowIndex уравнения.
            f[rowIndex] /= A.GetElement(rowIndex, rowIndex);

            // Обнуляем все элементы в столбце выше диагонального.
            for (ptrdiff_t nullRowIndex = rowIndex - 1; nullRowIndex > 0; nullRowIndex--)
                f[nullRowIndex] -= f[rowIndex] * A.GetElement(rowIndex, nullRowIndex);
        }
    }

    /**
     * @brief  Решает систему линейных уравнений методом Гаусса.
     * Система: Ax = f, где
     *     A - матрица системы.
     *     f - вектор свободных коэффициентов.
     *     x - вектор искомых переменных.
     * 
     * @tparam T Тип используемых чисел. Обычно либо double, либо float. Целочисленные типы имеют смысл, если нет деления.
     * 
     * @param A Матрица системы.
     * @param f Вектор свободных коэффициентов.
     * @param solution Вектор искомых переменных.
    */
    template <typename T>
    void SolveGauss(Matrix<T>& A, Vector<T>& f, Vector<T>& solution, bool findMaxElem = true)
    {
        ApplyForwardGaussMethod<T>(A, f, findMaxElem);
        ApplyBackwardGaussMethod<T>(A, f);
        solution = f;
    }
}

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///