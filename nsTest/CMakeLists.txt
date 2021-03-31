cmake_minimum_required(VERSION 3.5)
project(NekSenseiTest)

set(sources test.cxx ../DataAdaptor.cxx ../Bridge.cxx)
set(libs mpi sensei)

find_package(MPI REQUIRED)
find_package(SENSEI REQUIRED)

if($ENV{UNDERSCORE})
    add_definitions(-DUNDERSCORE)
endif()
if($ENV{UPCASE})
    add_definitions(-DUPCASE)
endif()

set(CMAKE_BUILD_TYPE Debug)

add_executable(nekSenseiTest ${sources})
target_link_libraries(nekSenseiTest ${libs})
