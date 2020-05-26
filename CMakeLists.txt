cmake_minimum_required(VERSION 3.15)
project(RCFC_kb)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(RCFC_kb
        classifier.cpp
        classifier.h
        codeFragment.cpp
        codeFragment.h
        configuration.h
        env.cpp
        env.h
        xtime.cpp
        xtime.h
        xcs.cpp
        xcs.h
        xcsMacros.cpp
        xcsMacros.h filter_list.cpp filter_list.h)
