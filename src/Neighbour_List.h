//
// Created by malcolm on 24/02/16.
//

#ifndef ANALYSIS_NEIGHBOUR_LIST_H
#define ANALYSIS_NEIGHBOUR_LIST_H

#include "Molecule.h"
#include "Frame.h"

class Neighbour_List {
private:
    std::vector<std::vector<int>> neighs;
public:
    Neighbour_List() : neighs(std::vector<std::vector<int>>{}) {};

    Neighbour_List(const Frame &);
    void update(const Frame &);
};


#endif //ANALYSIS_NEIGHBOUR_LIST_H
