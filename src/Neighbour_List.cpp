//
// Created by malcolm on 24/02/16.
//

#include "Neighbour_List.h"

Neighbour_List::Neighbour_List(const Frame &frame) : neighs(std::vector<std::vector<int>>{frame.num_atoms()}) {
    for (auto atom1: frame.particles) {
        for (auto atom2: frame.particles){
            if (dist(atom1.pos, atom2.pos) < (atom1.radius + atom2.radius)){
                neighs[atom1.id-1].push_back(atom2.id-1);
                neighs[atom2.id-1].push_back(atom1.id-1);
            }
        }
    }
}

void Neighbour_List::update(const Frame &frame) {
    neighs = std::vector<std::vector<int>>{frame.num_atoms()};
    for (auto atom1: frame.particles) {
        for (auto atom2: frame.particles){
            if (dist(atom1.pos, atom2.pos) < (atom1.radius + atom2.radius)){
                neighs[atom1.id-1].push_back(atom2.id-1);
                neighs[atom2.id-1].push_back(atom1.id-1);
            }
        }
    }
}