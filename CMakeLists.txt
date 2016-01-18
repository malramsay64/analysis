cmake_minimum_required(VERSION 3.3)
project(Analysis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "NO build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

add_subdirectory(test/lib/googletest)
include_directories(test/lib/googletest/include)

link_libraries(gtest pthread)

set(SOURCE_FILES
        src/functions.cpp
        src/functions.h
        src/Vector2d.cpp
        src/Vector2d.h
        )

set(TEST_FILES
        test/testVector2d.h test/testFunctions.h)

add_executable(analysis ${SOURCE_FILES})
add_executable(runTests test/runTests.cpp ${SOURCE_FILES} ${TEST_FILES})