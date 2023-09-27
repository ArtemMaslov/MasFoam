#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

#include <MasFoam/Derivative/Derivative.h>

using namespace MasFoam;
using namespace MasFoam::Derivative;

///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///
///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///

PYBIND11_MODULE(MasFoam, m)
{
    m.doc() = "MasFoam is a computational mathematic C++ library.";
    
    // Производные.
    m.def("DerivateRightO1",
          &MasFoam::Derivative::DerivateRightO1<double, const std::function<double(double)>&>, 
          "Правая односторонняя 1 производная 1 порядка, O(step).");
    m.def("DerivateLeftO1",
          &MasFoam::Derivative::DerivateLeftO1<double, const std::function<double(double)>&>,
          "Левая односторонняя 1 производная 1 порядка, O(step).");
    m.def("DerivateCentralO2",
          &MasFoam::Derivative::DerivateCentralO2<double, const std::function<double(double)>&>,
          "1 производная 2 порядка, O(step^2).");
    m.def("DerivateCentralO4",
          &MasFoam::Derivative::DerivateCentralO4<double, const std::function<double(double)>&>,
          "1 производная 4 порядка, O(step^4).");
    m.def("DerivateCentralO6",
          &MasFoam::Derivative::DerivateCentralO6<double, const std::function<double(double)>&>,
          "1 производная 6 порядка, O(step^6).");
}

///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///
///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///