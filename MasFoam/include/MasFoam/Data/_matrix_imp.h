#pragma once

#include <cassert>
#include <cstring>
#include <algorithm>

#include <MasFoam/Data/Matrix.h>
#include <MasFoam/Data/Vector.h>

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

namespace MasFoam
{
    template <typename T>
    Matrix<T>::Matrix(size_t columnCount, size_t rowCount) :
        ColumnsCount(columnCount),
        RowsCount(rowCount),
        Data(nullptr),
        RowPhysicalIndexes(RowsCount)
    {
        Data = new T[GetSize()];
        memset(Data, 0, GetSize());
        for (size_t st = 0; st < RowsCount; st++)
            RowPhysicalIndexes[st] = st;
    }

    template <typename T>
    Matrix<T>::Matrix(const Matrix<T>& matrix) : 
        ColumnsCount(matrix.ColumnsCount),
        RowsCount(matrix.RowsCount),
        Data(nullptr),
        RowPhysicalIndexes(matrix.RowPhysicalIndexes)
    {
        Data = new T[GetSize()];
        for (size_t st = 0; st < GetSize(); st++)
            Data[st] = matrix.Data[st];
    }

    template <typename T>
    Matrix<T>::Matrix(Matrix<T>&& matrix) :
        ColumnsCount(matrix.ColumnsCount),
        RowsCount(matrix.RowsCount),
        Data(matrix.Data),
        RowPhysicalIndexes(std::move(matrix.RowPhysicalIndexes))
    {
        matrix.ColumnsCount = 0;
        matrix.RowsCount    = 0;
        matrix.Data         = nullptr;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator = (const Matrix<T>& matrix)
    {
        if (this == &matrix)
            return *this;
        
        this->~Matrix();

        ColumnsCount = matrix.ColumnsCount;
        RowsCount    = matrix.RowsCount;
        Data         = new T[GetSize()];
        for (size_t st = 0; st < GetSize(); st++)
            Data[st] = matrix.Data[st];
        RowPhysicalIndexes = matrix.RowPhysicalIndexes;

        return *this;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator = (Matrix<T>&& matrix)
    {
        if (this == &matrix)
            return *this;
        
        this->~Matrix();

        ColumnsCount = matrix.ColumnsCount;
        RowsCount    = matrix.RowsCount;
        Data         = matrix.Data;
        RowPhysicalIndexes = std::move(matrix.RowPhysicalIndexes);

        matrix.ColumnsCount = 0;
        matrix.RowsCount    = 0;
        matrix.Data         = nullptr;

        return *this;
    }

    template <typename T>
    Matrix<T>::~Matrix()
    {
        if (Data)
            delete [] Data;
        ColumnsCount = 0;
        RowsCount    = 0;
        Data         = nullptr;
    }

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    template <typename T>
    T& Matrix<T>::GetElement(size_t index)
    {
        return Data[index];
    }

    template <typename T>
    const T& Matrix<T>::GetElement(size_t index) const
    {
        // Не константная версия метода не меняет внутреннее состояние объекта.
        return const_cast<Matrix<T>*>(this)->GetElement(index);
    }

    template <typename T>
    T& Matrix<T>::GetElement(size_t rowIndex, size_t columnIndex)
    {
        return Data[GetRowStartIndex(rowIndex) + columnIndex];
    }

    template <typename T>
    const T& Matrix<T>::GetElement(size_t rowIndex, size_t columnIndex) const
    {
        // Не константная версия метода не меняет внутреннее состояние объекта.
        return const_cast<Matrix<T>*>(this)->GetElement(rowIndex, columnIndex);
    }

    template <typename T>
    Vector<T> Matrix<T>::GetRow(size_t rowIndex)
    {
        return Vector{ColumnsCount, Data + GetRowStartIndex(rowIndex)};
    }

    template <typename T>
    const Vector<T> Matrix<T>::GetRow(size_t rowIndex) const
    {
        // Не константная версия метода не меняет внутреннее состояние объекта.
        return const_cast<Matrix<T>*>(this)->GetRow(rowIndex);
    }

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    template <typename T>
    size_t Matrix<T>::GetRowStartIndex(size_t rowIndex) const
    {
        return RowPhysicalIndexes[rowIndex] * ColumnsCount;
    }

    template <typename T>
    size_t Matrix<T>::GetSize() const
    {
        return ColumnsCount * RowsCount;
    }

    template <typename T>
    size_t Matrix<T>::GetColumnsCount() const
    {
        return ColumnsCount;
    }

    template <typename T>
    size_t Matrix<T>::GetRowsCount() const
    {
        return RowsCount;
    }

    template <typename T>
    void Matrix<T>::SwapRows(size_t row1, size_t row2)
    {
        // size_t value = RowPhysicalIndexes[row1];
        // RowPhysicalIndexes[row1] = RowPhysicalIndexes[row2];
        // RowPhysicalIndexes[row2] = value;
        std::swap(RowPhysicalIndexes[row1], RowPhysicalIndexes[row2]);
    }

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    template <typename T>
    Matrix<T>& Matrix<T>::operator *= (const T number)
    {
        const size_t size = GetSize();
        for (size_t st = 0; st < size; st++)
            Data[st] *= number;
        return *this;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator += (const T number)
    {
        const size_t size = GetSize();
        for (size_t st = 0; st < size; st++)
            Data[st] += number;
        return *this;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator -= (const T number)
    {
        const size_t size = GetSize();
        for (size_t st = 0; st < size; st++)
            Data[st] -= number;
        return *this;
    }
    
    template <typename T>
    Matrix<T>& Matrix<T>::operator += (const Matrix<T>& matrix)
    {
        for (size_t rowIndex = 0; rowIndex < RowsCount; rowIndex++)
        {
            T* resMatrixRow = Data + GetRowStartIndex(rowIndex);
            T* addMatrixRow = matrix.Data + matrix.GetRowStartIndex(rowIndex);
            for (size_t columnIndex = 0; columnIndex < ColumnsCount; columnIndex++)
                resMatrixRow[columnIndex] += addMatrixRow[columnIndex];
        }
        return *this;
    }

    template <typename T>
    Matrix<T>& Matrix<T>::operator -= (const Matrix<T>& matrix)
    {
        for (size_t rowIndex = 0; rowIndex < RowsCount; rowIndex++)
        {
            T* resMatrixRow = Data + GetRowStartIndex(rowIndex);
            T* addMatrixRow = matrix.Data + matrix.GetRowStartIndex(rowIndex);
            for (size_t columnIndex = 0; columnIndex < ColumnsCount; columnIndex++)
                resMatrixRow[columnIndex] -= addMatrixRow[columnIndex];
        }
        return *this;
    }

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    template <typename T>
    Vector<T> operator * (const Matrix<T>& matrix, const Vector<T>& vector)
    {
        // матрица * вектор = вектор.
        // Матрица и вектор должны быть согласованы.
        assert(matrix.ColumnsCount == vector.Size);

        const size_t vectorSize = vector.Size;
        const size_t resultVectorSize = matrix.RowsCount;
        Vector<T> result{resultVectorSize};

        for (size_t resultIndex = 0; resultIndex < resultVectorSize; resultIndex++)
        {
            result[resultIndex] += vector * matrix.GetRow(resultIndex);
        }

        return result;
    }

    template <typename T>
    Vector<T> operator * (const Vector<T>& vector, const Matrix<T>& matrix)
    {
        // вектор * матрица = вектор.
        // Матрица и вектор должны быть согласованы.
        assert(vector.Size == matrix.RowsCount);

        const size_t vectorSize = vector.Size;
        const size_t resultVectorSize = matrix.ColumnsCount;
        Vector<T> result{resultVectorSize};

        for (size_t resultIndex = 0; resultIndex < resultVectorSize; resultIndex++)
        {
            for (size_t scalarMulIndex = 0; scalarMulIndex < vectorSize; scalarMulIndex++)
                result[resultIndex] += vector[scalarMulIndex] * matrix.GetElement(scalarMulIndex, resultIndex);
        }

        return std::move(result);
    }

    template <typename T>
    Matrix<T> operator * (const Matrix<T>& m1, const Matrix<T>& m2)
    {
        // матрица * матрица = матрица.
        // Матрицы должны быть согласованы.
        assert(m1.ColumnsCount == m2.RowsCount);

        const size_t resRowsCount    = m1.RowsCount;
        const size_t resColumnsCount = m2.ColumnsCount;
        const size_t scalarSize      = m1.ColumnsCount;
        Matrix<T> result{resRowsCount, resColumnsCount};

        for (size_t rowIndex = 0; rowIndex < resRowsCount; rowIndex++)
        {
            const Vector<T> row = m1.GetRow(rowIndex);
            for (size_t columnIndex = 0; columnIndex < resColumnsCount; columnIndex++)
            {
                T elem = 0;
                for (size_t scalarIndex = 0; scalarIndex < scalarSize; scalarIndex++)
                    elem += row[scalarIndex] + m2.GetElement(scalarIndex, columnIndex);
                result.GetElement(rowIndex, columnIndex) = elem;
            }
        }

        return std::move(result);
    }

}

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///