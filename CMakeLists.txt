cmake_minimum_required(VERSION 3.5)
project(NekSenseiDriver)

list(APPEND CMAKE_MODULE_PATH "${SENSEI_DIR}")
list(APPEND CMAKE_MODULE_PATH "${SENSEI_DIR}/lib/cmake")

set(sources DataAdaptor.cxx Bridge.cxx)

find_package(MPI REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development)
find_package(SENSEI REQUIRED)

set(libs mpi ${Python3_LIBRARIES} sensei)
message(${Python3_LIBRARIES})

if($ENV{UNDERSCORE})
    add_definitions(-DUNDERSCORE)
endif()
if($ENV{UPCASE})
    add_definitions(-DUPCASE)
endif()

add_definitions()
add_library(neksensei STATIC ${sources})
target_link_libraries(neksensei ${libs})


set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
set(CMAKE_BUILD_TYPE Release)
install(TARGETS neksensei DESTINATION lib)
install(FILES Bridge.h DESTINATION include)
