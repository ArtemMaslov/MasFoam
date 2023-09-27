#pragma once

#include <functional>

///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///
///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///

namespace MasFoam::Derivative
{
    template <typename T, typename functType>
    T DerivateRightO1(functType funct, T arg, T step)
    {
        return (funct(arg + step) - funct(arg)) / step;
    }

    template <typename T, typename functType>
    T DerivateLeftO1(functType funct, T arg, T step)
    {
        return (funct(arg) - funct(arg - step)) / step;
    }

    template <typename T, typename functType>
    T DerivateCentralO2(functType funct, T arg, T step)
    {
        return (funct(arg + step) - funct(arg - step)) / (2 * step);
    }

    template <typename T, typename functType>
    T DerivateCentralO4(functType funct, T arg, T step)
    {
        return 4 * (funct(arg +     step) - funct(arg -     step)) / (3 * 2 * step) -
                   (funct(arg + 2 * step) - funct(arg - 2 * step)) / (3 * 4 * step);
    }

    template <typename T, typename functType>
    T DerivateCentralO6(functType funct, T arg, T step)
    {
        return 3 * (funct(arg +     step) - funct(arg -     step)) / (2 * 2 * step) -
               3 * (funct(arg + 2 * step) - funct(arg - 2 * step)) / (5 * 4 * step) +
                   (funct(arg + 3 * step) - funct(arg - 3 * step)) / (10 * 6 * step);
    }
}

///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///
///***///***///---\\\***\\\***\\\___///***___***\\\___///***///***///---\\\***\\\***///