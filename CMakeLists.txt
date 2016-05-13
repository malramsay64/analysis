cmake_minimum_required(VERSION 2.8)
project(Analysis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "NO build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()


set(SOURCE_FILES
        src/functions.cpp
        src/functions.h
        src/Vector2d.cpp
        src/Vector2d.h
        src/Particle.cpp
        src/Particle.h
        src/Molecule.cpp
        src/Molecule.h
        )
    FILE(GLOB SOURCE_FILES src/*)

set(TEST_FILES
        test/testVector2d.h
        test/testFunctions.h
        test/testParticle.h
        test/testMolecule.h)

if (NOT CMAKE_BUILD_TYPE STREQUAL "Release")

    # Include google test
    #add_subdirectory(test/lib/googletest)
    #include_directories(test/lib/googletest/include)
    #link_libraries(gtest pthread)

    #add_executable(runTests test/runTests.cpp ${SOURCE_FILES} ${TEST_FILES})

endif()

# Include stats
include_directories($ENV{HOME}/.local/include)
link_directories($ENV{HOME}/.local/lib)

link_libraries(stats)

add_executable(analysis ${SOURCE_FILES})

install(TARGETS analysis
        CONFIGURATIONS Release
        DESTINATION bin
        )
