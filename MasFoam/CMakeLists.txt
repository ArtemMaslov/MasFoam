cmake_minimum_required(VERSION 3.22)
project(MasFoam
        VERSION 0.0.1)

include(GNUInstallDirs)
include(CTest)
find_package(pybind11 REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC -Wno-unused-variable")

set(IncludesDirPath ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(Includes 
    Data/Matrix.h
    Data/_matrix_imp.h
    Data/Vector.h
    Data/_vector_imp.h
    Derivative/Derivative.h
    LinearEquations/GaussSolver.h
)

set(PythonModuleSources
    ${Sources}
    PythonModule.cpp)

# Модуль Python.
pybind11_add_module(MasFoam ${PythonModuleSources})
target_include_directories(MasFoam PRIVATE ${IncludesDirPath})
set_target_properties(MasFoam PROPERTIES VERSION ${PROJECT_VERSION})

install(TARGETS MasFoam
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Тесты.
enable_testing()
add_subdirectory(UnitTests)

# Динамическая библиотека C++.
#[===[
add_library(masfoam ${Sources})
target_include_directories(masfoam PRIVATE ${IncludesDirPath})
set_target_properties(masfoam PROPERTIES VERSION ${PROJECT_VERSION})

install(TARGETS masfoam
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
]===]