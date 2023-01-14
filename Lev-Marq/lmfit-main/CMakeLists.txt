cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

option(LMFIT_CPPTEST "Build cpp demos" ON)

if(LMFIT_CPPTEST)
    project(lmfit VERSION 9.0.0 LANGUAGES C CXX)
else()
    project(lmfit VERSION 9.0.0 LANGUAGES C)
endif()

## Options.

option(FITTEST "Build with FitTest" OFF)
if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build as shared library" ON)
endif()
if(NOT DEFINED LIB_MAN)
    option(LIB_MAN "Build library man pages" ON)
endif()
option(WERROR "Treat warnings as errors" OFF)

## Compiler settings.

set(CMAKE_C_STANDARD 99)
if(LMFIT_CPPTEST)
    set(CMAKE_CXX_STANDARD 17)
endif()

if(WIN32)
    add_compile_options(/wd4267 /wd4710 /wd4820)
#    if(BUILD_SHARED_LIBS)
#        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
#        set(WINDOWS_EXPORT_ALL_SYMBOLS ON)
#    endif()
#    add_compile_options(-export-all-symbols)
else()
    option(PEDANTIC "Compile with pedantic warnings" ON)
    if(PEDANTIC)
        add_compile_options(-pedantic -Wall)
    endif()
endif()
add_compile_options(-O2)
if(WERROR)
    add_compile_options(-Werror)
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Install settings.

configure_file(lmfit.pc.in lmfit.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lmfit.pc DESTINATION lib/pkgconfig/)

include(CTest)

add_subdirectory(lib)
add_subdirectory(demo)
add_subdirectory(test)
if(LMFIT_CPPTEST)
    add_subdirectory(democpp)
endif()
if (LIB_MAN)
    add_subdirectory(man)
endif()
if (FITTEST)
    find_package(FitTest REQUIRED)
    add_subdirectory(fittest)
endif()
