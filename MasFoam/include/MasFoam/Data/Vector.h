#pragma once

#include <cstddef>

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

namespace MasFoam
{
    template <typename T>
    class Matrix;

    template <typename T>
    class Vector
    {
        friend class Matrix<T>;

    public:
        // Создаёт вектор в куче.
        Vector(size_t size);
        // Создаёт "окно для просмотра" уже созданного вектора.
        Vector(size_t size, T* data);

        Vector(const Vector<T>& vector);
        Vector(Vector<T>&& vector);
        Vector<T>& operator = (const Vector<T>& vector);
        Vector<T>& operator = (Vector<T>&& vector);

        ~Vector();

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

        size_t GetSize();

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

        T& operator [] (size_t index);

        Vector<T>& operator *= (T koef);
        Vector<T>& operator /= (T koef);
        Vector<T>& operator += (const Vector<T>& vector);
        Vector<T>& operator -= (const Vector<T>& vector);

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    private:
        size_t Size;
        bool   Allocated;
        T*     Data;

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    template <typename U>
    friend U operator * (const Vector<U>& v1, const Vector<U>& v2);

    template <typename U>
    friend Vector<U> operator * (const Matrix<U>& matrix, const Vector<U>& vector);

    template <typename U>
    friend Vector<U> operator * (const Vector<U>& vector, const Matrix<U>& matrix);
    };

    // Скалярное произведение.
    template <typename T>
    T operator * (const Vector<T>& v1, const Vector<T>& v2);

    using VectorD = Vector<double>;
}

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

#include "_vector_imp.h"

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///