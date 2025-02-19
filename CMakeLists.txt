cmake_minimum_required(VERSION 3.20)

project(polikit
    VERSION 0.3
    DESCRIPTION "Polikit: Polyhedral analyzing toolkits."
    LANGUAGES Fortran)

enable_language(Fortran)

add_subdirectory(src)
