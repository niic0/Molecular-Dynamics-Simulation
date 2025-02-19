cmake_minimum_required(VERSION 3.10)
project(SimulationProject)

# 🔹 Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 🔹 Compilation flags for different build types
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -fno-omit-frame-pointer -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -funroll-loops -march=native -fno-omit-frame-pointer -g3 -fopenmp")

# 🔹 Include the header directory
include_directories(include)

# 🔹 Define source files
set(SOURCES 
    src/simulation.cpp
    src/utils.cpp
    main.cpp
)

# 🔹 Create the executable
add_executable(simulation ${SOURCES})

