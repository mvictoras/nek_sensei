cmake_minimum_required(VERSION 3.5)
project(NekSenseiDriver)

set(sources DataAdaptor.cxx Bridge.cxx)
set(libs mpi sensei)

find_package(MPI REQUIRED)
find_package(SENSEI REQUIRED)

if($ENV{UNDERSCORE})
    add_definitions(-DUNDERSCORE)
endif()
if($ENV{UPCASE})
    add_definitions(-DUPCASE)
endif()

add_library(neksensei STATIC ${sources})
target_link_libraries(neksensei ${libs})

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
set(CMAKE_BUILD_TYPE Debug)
install(TARGETS neksensei DESTINATION lib)
install(FILES Bridge.h DESTINATION include)
