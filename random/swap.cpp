//
//  swap.cpp
//  analysis
//
//  Created by Malcolm Ramsay on 27/04/2015.
//  Copyright (c) 2015 Malcolm Ramsay. All rights reserved.
//

#include "swap.h"

using namespace std;

int swap_neighbours(Frame * frame, int molindex){
    molecule * m1 = &frame->molecules.at(molindex);
    vector<int> m1_neigh = vector<int>();
    
    // Find all molecules in antiparallel configuration
    for (auto mol: m1->my_neighbours){
        if (order_type(m1, mol, frame) == 3 || order_type(m1, mol, frame) == 4){
             m1_neigh.push_back(mol->index());
        }
    }
    
    if (m1_neigh.size() == 0){
        return 0;
    }

    // Randomly select from list
    //srand(0);
    int r = rand() % m1_neigh.size();
    molecule * m2 = &frame->molecules.at(m1_neigh.at(r));
    vector<int> m1_atoms = {m1->atoms.at(0)->index(), m1->atoms.at(1)->index()};
    vector<int> m2_atoms = {m2->atoms.at(0)->index(), m2->atoms.at(1)->index()};
    
    // Do the swap bit
    m1->delete_mol_neighbours();
    m2->delete_mol_neighbours();
    
    m1->atoms = {&frame->particles.at(m1_atoms.at(0)), &frame->particles.at(m2_atoms.at(1))};
    m2->atoms = {&frame->particles.at(m2_atoms.at(0)), &frame->particles.at(m1_atoms.at(1))};
    
    
    // update molid
    frame->particles.at(m2_atoms.at(1)).molid = m1->id;
    frame->particles.at(m1_atoms.at(1)).molid = m2->id;
    
    // Set new particle neighbours
    for (auto &n: frame->particles.at(m2_atoms.at(1)).my_neighbours){
        if (n->id == frame->particles.at(m1_atoms.at(0)).id){
            n = &frame->particles.at(m2_atoms.at(0));
        }
    }
    for (auto &n: frame->particles.at(m2_atoms.at(0)).my_neighbours){
        if (n->id == frame->particles.at(m1_atoms.at(1)).id){
            n = &frame->particles.at(m2_atoms.at(1));
        }
    }
    for (auto &n: frame->particles.at(m1_atoms.at(1)).my_neighbours){
        if (n->id == frame->particles.at(m2_atoms.at(0)).id){
            n = &frame->particles.at(m1_atoms.at(0));
        }
    }
    for (auto &n: frame->particles.at(m1_atoms.at(0)).my_neighbours){
        if (n->id == frame->particles.at(m2_atoms.at(1)).id){
            n = &frame->particles.at(m1_atoms.at(1));
        }
    }

                        
    // Add new neighbours
    update_neighbours(m1, frame);
    update_neighbours(m2, frame);
    
    // Update mols
    m1->calc_COM();
    m2->calc_COM();
    
    m1->set_orientation(angle(m1, frame));
    m1->set_orientation(angle(m2, frame));
    
    return 1;
}

int update_neighbours(molecule * mol, Frame * frame){
    for (auto atom: mol->atoms){
        for (auto neigh: atom->my_neighbours){
            mol->add_neighbour(&frame->molecules.at(neigh->mol_index()));
        }
    }
    return 0;
}