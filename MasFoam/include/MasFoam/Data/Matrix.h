#pragma once

#include <cstddef>
#include <vector>

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

namespace MasFoam
{
    template <typename T>
    class Vector;

    template <typename T>
    class Matrix
    {
    public:
        Matrix(size_t columnCount, size_t rowCount);

        Matrix(const Matrix<T>& matrix);
        Matrix(Matrix<T>&& matrix);
        Matrix<T>& operator = (const Matrix<T>& matrix);
        Matrix<T>& operator = (Matrix<T>&& matrix);

        ~Matrix();

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

        T& GetElement(size_t index);
        const T& GetElement(size_t index) const;

        T& GetElement(size_t rowIndex, size_t columnIndex);
        const T& GetElement(size_t rowIndex, size_t columnIndex) const;

        Vector<T> GetRow(size_t rowIndex);
        const Vector<T> GetRow(size_t rowIndex) const;

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    private:
        size_t GetRowStartIndex(size_t rowIndex) const;

    public:

        size_t GetSize() const;

        size_t GetColumnsCount() const;

        size_t GetRowsCount() const;

        void SwapRows(size_t row1, size_t row2);

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

        Matrix<T>& operator *= (const T number);
        Matrix<T>& operator += (const T number);
        Matrix<T>& operator -= (const T number);
        Matrix<T>& operator += (const Matrix<T>& matrix);
        Matrix<T>& operator -= (const Matrix<T>& matrix);

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    private:
        // Количество столбцов.
        size_t ColumnsCount;
        // Количество строк.
        size_t RowsCount;

        // Элементы матрицы.
        T* Data;

        /**
         * Задаёт порядок строк матрицы A.
         * 
         * Для нормальной матрицы RowPhysicalIndexes = [0, 1, 2, ..., RowsCount - 1].
         * Для оптимизации, чтобы при перестановках строк выполнялась не реальная перестановка, которая 
         * занимает много времени для больших матриц, меняются местами соответствующие элементы RowPhysicalIndexes.
         * То есть RowPhysicalIndexes указывает на номер физической строки матрицы A, когда для пользователя все
         * индексы строк являются логическими.
        */
        std::vector<size_t> RowPhysicalIndexes;

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    template <typename U>
    friend Vector<U> operator * (const Matrix<U>& matrix, const Vector<U>& vector);

    template <typename U>
    friend Vector<U> operator * (const Vector<U>& vector, const Matrix<U>& matrix);

    template <typename U>
    friend Matrix<U> operator * (const Matrix<U>& m1, const Matrix<U>& m2);
    };

    template <typename T>
    Vector<T> operator * (const Matrix<T>& matrix, const Vector<T>& vector);

    template <typename T>
    Vector<T> operator * (const Vector<T>& vector, const Matrix<T>& matrix);

    template <typename T>
    Matrix<T> operator * (const Matrix<T>& m1, const Matrix<T>& m2);

    using MatrixD = Matrix<double>;
}

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

#include "_matrix_imp.h"

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///