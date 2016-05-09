cmake_minimum_required(VERSION 2.8)
project(Analysis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "NO build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(SOURCE_FILES
        src/functions.cpp   src/functions.h
        src/Particle.cpp    src/Particle.h
        src/Molecule.cpp    src/Molecule.h
        src/Frame.cpp       src/Frame.h
        src/neighbours.cpp  src/neighbours.h
        src/Vector.cpp      src/Vector.h
        src/Compute.cpp     src/Compute.h
        src/Search.cpp      src/Search.h
        src/Neighbour_List.cpp src/Neighbour_List.h
        src/input.cpp       src/input.h)

set(TEST_FILES
        test/testFunctions.h
        test/testParticle.h
        test/testMolecule.h
        test/testFrame.h
        test/testVector.h
        test/testCompute.h
        test/testSearch.h test/testInput.h test/testNeighbours.h)

# Include stats
include_directories($ENV{HOME}/.local/include)
link_directories($ENV{HOME}/.local/lib)

link_libraries(stats)

if (NOT CMAKE_BUILD_TYPE STREQUAL "Release")

    # Include google test
    add_subdirectory(lib/googletest)
    include_directories(lib/googletest/include)
    link_libraries(gtest pthread stats)

    add_executable(runTests test/runTests.cpp ${SOURCE_FILES} ${TEST_FILES})

endif()

add_executable(analysis src/main.cpp ${SOURCE_FILES})

install(TARGETS analysis CONFIGURATIONS Release DESTINATION bin)
