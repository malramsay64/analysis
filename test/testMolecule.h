//
// Created by malcolm on 19/01/16.
//

#ifndef ANALYSIS_TESTMOLECULE_H
#define ANALYSIS_TESTMOLECULE_H

#include <gtest/gtest.h>
#include "../src/Molecule.h"

/* Create a molecule with:
 * COM = (1,1)
 * Mass = 3
 * 3 particles
 * 2 particle types
 */
Molecule create_molecule(){
    Molecule m1{};
    Particle *p1, *p2, *p3;
    p1 = new Particle{};
    p2 = new Particle{};
    p3 = new Particle{};
    p1->set_pos(Vector2d{2,1});
    p1->id = 1;
    p1->type = 2;
    p2->set_pos(Vector2d{1,1});
    p2->id = 2;
    p2->type = 1;
    p3->set_pos(Vector2d{0,1});
    p3->id = 3;
    p3->type = 2;
    m1.atoms.push_back(p1);
    m1.atoms.push_back(p2);
    m1.atoms.push_back(p3);
    return m1;
}

void move_particles(Molecule &m){
    for (auto &atom :m.atoms) {
        atom->pos += 0.1;
    }
}

void delete_molecule(Molecule &m){
    for (auto &a: m.atoms){
        delete a;
    }
}

TEST(Molecule, Initialisation){
    Molecule m{};
    EXPECT_EQ(0, m.contacts);
    EXPECT_EQ(0, m.type);
    EXPECT_EQ(0, m.id);
    EXPECT_EQ(0, m.get_rotation());
    EXPECT_EQ(0, m.atoms.size());
}

TEST(Molecule, AddNeighbour){
    Molecule m1{}, m2{}, m3{}, m4{};
    m1.id = 1;
    m2.id = 2;
    m3.id = 3;
    m4.id = 4;
    m1.add_neighbour(&m2);
    m1.add_neighbour(&m3);
    m1.add_neighbour(&m4);
    m1.add_neighbour(&m2);
    m1.add_neighbour(&m2);
    EXPECT_EQ(5, m1.num_contacts());
    EXPECT_EQ(3, m1.num_neighbours());
}

TEST(Molecule, Mass){
    Molecule m1 = create_molecule();
    EXPECT_DOUBLE_EQ(3,m1.get_mass());
    for (auto &atom :m1.atoms) {
        atom->mass = 2.5;
    }
    EXPECT_DOUBLE_EQ(7.5, m1.get_mass());
    delete_molecule(m1);
}

TEST(Molecule, COM){
    Molecule m1 = create_molecule();
    EXPECT_EQ(Vector2d(1,1), m1.get_COM());
    double pi2 = PI*2;
    double x,y, x_m, y_m;
    for (int i=1; i<32; i++){
        move_particles(m1);
        m1.update_COM();
        x = y = fmod(1.+i/10.,pi2);
        x_m = y_m = 1+i/10.;
        EXPECT_NEAR(x, m1.get_COM().x, 1e-14);
        EXPECT_NEAR(y, m1.get_COM().y, 1e-14);
        EXPECT_NEAR(x_m, m1.moved_COM().x, 1e-14);
        EXPECT_NEAR(y_m, m1.moved_COM().y, 1e-14);
    }
    delete_molecule(m1);
}

TEST(Molecule, CopyConstructor){
    Molecule m1 = create_molecule();
    Molecule m2{m1};
    move_particles(m1);
    EXPECT_EQ(Vector2d(1.1,1.1), m2.get_COM());
    m1.id = 1;
    m2.id = 2;
    EXPECT_EQ(1,m1.id);
    EXPECT_EQ(2,m2.id);
    m2.atoms = std::vector<Particle *>{};
    EXPECT_EQ(3, m1.num_particles());
    EXPECT_EQ(0, m2.num_particles());
    delete_molecule(m1);
}



#endif //ANALYSIS_TESTMOLECULE_H
