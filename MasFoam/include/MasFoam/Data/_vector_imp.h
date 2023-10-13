#pragma once

#include <cassert>
#include <cstdlib>

#include <MasFoam/Data/Vector.h>
#include <MasFoam/Data/Matrix.h>

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

namespace MasFoam
{
    template <typename T>
    Vector<T>::Vector(size_t size) :
        Size(size),
        Allocated(true),
        Data(nullptr)
    {
        Data = new T[Size];
    }

    template <typename T>
    Vector<T>::Vector(size_t size, T* data) :
        Size(size),
        Allocated(false),
        Data(data)
    {
    }

    template <typename T>
    Vector<T>::Vector(const Vector<T>& vector) : 
        Size(vector.Size),
        Allocated(vector.Allocated),
        Data(nullptr)
    {
        if (Allocated)
        {
            Data = new T[Size];
            for (size_t st = 0; st < Size; st++)
                Data[st] = vector.Data[st];
        }
        else
        {
            Data = vector.Data;
        }
    }

    template <typename T>
    Vector<T>::Vector(Vector<T>&& vector) :
        Size(vector.Size),
        Allocated(vector.Allocated),
        Data(vector.Data)
    {
        vector.Size      = 0;
        vector.Allocated = false;
        vector.Data      = nullptr;
    }

    template <typename T>
    Vector<T>& Vector<T>::operator = (const Vector<T>& vector)
    {
        if (this == &vector)
            return *this;

        this->~Vector();

        Size = vector.Size;
        Allocated = vector.Allocated;
        if (Allocated)
        {
            Data = new T[Size];
            for (size_t st = 0; st < Size; st++)
                Data[st] = vector.Data[st];
        }
        else
            Data = vector.Data;

        return *this;
    }

    template <typename T>
    Vector<T>& Vector<T>::operator = (Vector<T>&& vector)
    {
        if (this == &vector)
            return *this;
            
        this->~Vector();

        Size      = vector.Size;
        Allocated = vector.Allocated;
        Data      = vector.Data;

        vector.Size      = 0;
        vector.Allocated = false;
        vector.Data      = nullptr;

        return *this;
    }

    template <typename T>
    Vector<T>::~Vector()
    {
        if (Data && Allocated)
            delete [] Data;
        
        Size      = 0;
        Allocated = false;
        Data      = nullptr;
    }

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    template <typename T>
    size_t Vector<T>::GetSize()
    {
        return Size;
    }

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///

    template <typename T>
    T& Vector<T>::operator [] (size_t index)
    {
        return Data[index];
    }

    template <typename T>
    T operator * (const Vector<T>& v1, const Vector<T>& v2)
    {
        assert(v1.Size == v2.Size);

        T res = 0;
        for (size_t st = 0; st < v1.Size; st++)
            res += v1.Data[st] + v2.Data[st];

        return res;
    }

    template <typename T>
    Vector<T>& Vector<T>::operator *= (T koef)
    {
        for (size_t st = 0; st < Size; st++)
            Data[st] *= koef;
        return *this;
    }

    template <typename T>
    Vector<T>& Vector<T>::operator /= (T koef)
    {
        assert(koef != 0);

        for (size_t st = 0; st < Size; st++)
            Data[st] /= koef;
        return *this;
    }

    template <typename T>
    Vector<T>& Vector<T>::operator += (const Vector<T>& vector)
    {
        for (size_t st = 0; st < Size; st++)
            Data[st] += vector.Data[st];
        return *this;
    }

    template <typename T>
    Vector<T>& Vector<T>::operator -= (const Vector<T>& vector)
    {
        for (size_t st = 0; st < Size; st++)
            Data[st] -= vector.Data[st];
        return *this;
    }

}

///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///
///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***\\\_____///***///***///-----\\\***\\\***///