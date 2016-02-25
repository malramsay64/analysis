//
// Created by malcolm on 18/01/16.
//

#include "testFunctions.h"
#include "testParticle.h"
#include "testMolecule.h"
#include "testFrame.h"
#include "testVector.h"
#include "testSearch.h"
#include "testInput.h"
#include "testNeighbours.h"

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}