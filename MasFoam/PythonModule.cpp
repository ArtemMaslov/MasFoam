#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>

#include <MasFoam/Data/Vector.h>
#include <MasFoam/Data/Matrix.h>
#include <MasFoam/Derivative/Derivative.h>
#include <MasFoam/LinearEquations/GaussSolver.h>

using namespace MasFoam;

namespace py = pybind11;

///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///
///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///

PYBIND11_MODULE(MasFoam, m)
{
    m.doc() = "MasFoam is a computational mathematic C++ library.";
    
    // Производные.
    m.def("DerivateRightO1",
          &MasFoam::DerivateRightO1<double, const std::function<double(double)>&>, 
          "Правая односторонняя 1 производная 1 порядка, O(step).");
    m.def("DerivateLeftO1",
          &MasFoam::DerivateLeftO1<double, const std::function<double(double)>&>,
          "Левая односторонняя 1 производная 1 порядка, O(step).");
    m.def("DerivateCentralO2",
          &MasFoam::DerivateCentralO2<double, const std::function<double(double)>&>,
          "1 производная 2 порядка, O(step^2).");
    m.def("DerivateCentralO4",
          &MasFoam::DerivateCentralO4<double, const std::function<double(double)>&>,
          "1 производная 4 порядка, O(step^4).");
    m.def("DerivateCentralO6",
          &MasFoam::DerivateCentralO6<double, const std::function<double(double)>&>,
          "1 производная 6 порядка, O(step^6).");

    // Вектор.
    py::class_<VectorD>(m, "VectorD")
        .def(py::init<size_t>())
        .def("GetSize", &VectorD::GetSize)
        .def("__getitem__", &VectorD::operator[])
        .def("__setitem__", [](VectorD &self, unsigned index, double value) { self[index] = value; })
        
        .def(py::self *= double())
        .def(py::self /= double())
        .def(py::self += py::self)
        .def(py::self -= py::self);
        
    // Матрица.
    py::class_<MatrixD>(m, "MatrixD")
        .def(py::init<size_t, size_t>())

        .def("GetSize", &MatrixD::GetSize)
        .def("GetElement", static_cast<double& (MatrixD::*)(size_t)>(&MatrixD::GetElement))
        .def("SetElement", [](MatrixD &self, size_t index, double value) 
            { self.GetElement(index) = value; })
        .def("GetElement", static_cast<double& (MatrixD::*)(size_t, size_t)>(&MatrixD::GetElement))
        .def("SetElement", [](MatrixD &self, size_t rowIndex, size_t columnIndex, double value) 
            { self.GetElement(rowIndex, columnIndex) = value; })
        .def("GetRow", static_cast<VectorD (MatrixD::*)(size_t)>(&MatrixD::GetRow))

        .def("GetSize", &MatrixD::GetSize)
        .def("GetColumnsCount", &MatrixD::GetColumnsCount)
        .def("GetRowsCount", &MatrixD::GetRowsCount)
        .def("SwapRows", &MatrixD::SwapRows)

        .def(py::self *= double())
        .def(py::self += double())
        .def(py::self -= double())
        .def(py::self += py::self)
        .def(py::self -= py::self);

    // Метод Гаусса решения СЛАУ.
    m.def("ApplyForwardGaussMethod",
          &MasFoam::ApplyForwardGaussMethod<double>,
          "Применяет прямой ход метода Гаусса к системе линейных уравнений.",
          py::arg("A"),
          py::arg("f"),
          py::kw_only(),
          py::arg("findMaxElem") = true);
    m.def("ApplyBackwardGaussMethod",
          &MasFoam::ApplyBackwardGaussMethod<double>);
    m.def("SolveGauss",
          &MasFoam::SolveGauss<double>,
          py::arg("A"),
          py::arg("f"),
          py::arg("solution"),
          py::kw_only(),
          py::arg("findMaxElem") = true);
}

///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///
///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///