cmake_minimum_required(VERSION 3.28)
project(centro_symmetric)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_executable(run_mc run_mc_centro_symmetric.cpp ./lattice/lattice.cpp ./lattice/lattice.hpp
        ./mc_subroutine/mc_read_load_compute.cpp ./mc_subroutine/mc_read_load_compute.hpp)


set_source_files_properties(
        mc_subroutine/mc_read_load_compute.cpp
        PROPERTIES

        COMPILE_OPTIONS "-O3"
)


set_source_files_properties(
        lattice/lattice.cpp
        PROPERTIES

        COMPILE_OPTIONS "-O3"
)

set(Python3_INCLUDE_DIRS "/usr/include/python3.12")
set(Python3_LIBRARIES "/usr/lib/x86_64-linux-gnu/libpython3.12.so")
include_directories(${Python3_INCLUDE_DIRS})
set(Boost_NO_SYSTEM_PATHS YES)
if (Boost_NO_SYSTEM_PATHS)
    set(Boost_ROOT "~/cppLocal/boostBuild")
    set(Boost_INCLUDE_DIR "~/cppLocal/boostBuild/include")
    set(Boost_LIBRARY_DIR "~/cppLocal/boostBuild/lib")
endif (Boost_NO_SYSTEM_PATHS)

find_package(Boost  REQUIRED python numpy filesystem)
if(Boost_FOUND)
    message(STATUS "boost root is: ${Boost_ROOT}")
    message(STATUS "boost include path is : ${Boost_INCLUDE_DIR}")
    message(STATUS "boost library path is : ${Boost_LIBRARY_DIR}")
    message(STATUS "boost libraries is : ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})
else()
    message(WARNING "boost not found.")
endif()

target_link_libraries(run_mc Boost::filesystem)

target_link_libraries(run_mc Boost::python)
target_link_libraries(run_mc Boost::numpy)
#python
target_link_libraries(run_mc ${Python3_LIBRARIES})