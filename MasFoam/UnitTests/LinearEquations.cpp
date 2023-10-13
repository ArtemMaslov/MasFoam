#include <stdio.h>

#include <MasFoam/LinearEquations/GaussSolver.h>

using namespace MasFoam;

int main()
{
    MatrixD A{100, 100};

    for (size_t st = 0; st < 100; st++)
        A.GetElement(0, st) = st;

    for (size_t st = 1; st < 99; st++)
    {
        A.GetElement(st, st - 1) = 1;
        A.GetElement(st, st)     = 10;
        A.GetElement(st, st + 1) = 1;
    }
    A.GetElement(99, 98) = 1;
    A.GetElement(99, 99) = 1;

    VectorD f(100);

    for (size_t st = 0; st < 100; st++)
        f[st] = 100 - st;

    MatrixD _A = A;
    VectorD _f = f;

    VectorD solution{100};
    SolveGauss(A, f, solution);

    VectorD residual = _A * solution;
    residual -= _f;

    return 0;
}