cmake_minimum_required(VERSION 3.20)

project(polikit
    VERSION 0.4
    DESCRIPTION "Polikit: Polyhedral analyzing toolkits."
    LANGUAGES Fortran)

enable_language(Fortran)

option(DEBUG "Debug compile options" OFF)

if(DEBUG)
    add_compile_options(-fbounds-check)
    message(STATUS "Building in Debug mode")
else()
    add_compile_options(-O3)
    message(STATUS "Building in Release mode")
endif()

find_package(OpenMP)

# Get the number of CPU cores
if(NOT DEFINED OMP_THREAD_COUNT)
    cmake_host_system_information(RESULT CORE_COUNT QUERY NUMBER_OF_LOGICAL_CORES)
    set(OMP_THREAD_COUNT ${CORE_COUNT})
    message(STATUS "Detected CPU cores: ${CORE_COUNT}")
endif()

if(OpenMP_Fortran_FOUND)
    message(STATUS "Building with OpenMP on ${OMP_THREAD_COUNT} threads")

    # Set the number of threads as a compile definition
    add_compile_definitions(OMP_THREAD_COUNT=${OMP_THREAD_COUNT})

    add_compile_options(-fopenmp)
    add_link_options(-fopenmp)
    message(STATUS "Building with OMP")
endif()

file(GLOB SRC src/*.f90)

add_executable(polikit ${COMPILE_OPTIONS} ${SRC})
