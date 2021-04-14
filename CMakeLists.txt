cmake_minimum_required(VERSION 3.15)
project(XCS-IMG)

set(CMAKE_CXX_STANDARD 14)
SET(GCC_COVERAGE_COMPILE_FLAGS "-ggdb -g3 -rdynamic")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
include_directories(.)

add_executable(XCS-IMG
        classifier.cpp
        classifier.h
        codeFragment.cpp
        codeFragment.h
        configuration.h
        env.cpp
        env.h
        xcs.cpp
        xcs.h
        xcsMacros.cpp
        xcsMacros.h filter_list.cpp filter_list.h cf_list.h cf_list.cpp)
