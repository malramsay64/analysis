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
class MoleculeTest : public testing::Test {
protected:
    Molecule m1{};
    Particle *p1, *p2, *p3;

    virtual void SetUp(){
        p1 = new Particle{};
        p2 = new Particle{};
        p3 = new Particle{};
        p1->set_pos(Vector{2, 1, 0});
        p1->id = 1;
        p1->type = 2;
        p2->set_pos(Vector{1, 1, 0});
        p2->id = 2;
        p2->type = 1;
        p3->set_pos(Vector{0, 1, 0});
        p3->id = 3;
        p3->type = 2;
        m1.atoms.push_back(p1);
        m1.atoms.push_back(p2);
        m1.atoms.push_back(p3);
    }
    virtual void TearDown(){
        delete p1;
        delete p2;
        delete p3;
    };

    void move_particles(double dist=0.1) {
        for (auto &atom :m1.atoms) {
            atom->pos[0] += dist;
            atom->pos[1] += dist;
        }
    }

    void rot_particles() {
        for (auto  &atom: m1.atoms){
            Vector r = (atom->pos - m1.get_COM());
            double theta = atan2(r)+0.1;
            atom -> pos = Vector{std::sin(theta), std::cos(theta)};
        }
    }
};

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

TEST_F(MoleculeTest, Mass){
    EXPECT_DOUBLE_EQ(3,m1.get_mass());
    for (auto &atom :m1.atoms) {
        atom->mass = 2.5;
    }
    EXPECT_DOUBLE_EQ(7.5, m1.get_mass());
}

TEST_F(MoleculeTest, COM) {
    EXPECT_DOUBLE_EQ(1, m1.get_COM()[0]);
    EXPECT_DOUBLE_EQ(1, m1.get_COM()[1]);
}

TEST_F(MoleculeTest, UpdateCOM){
    double pi2 = PI*2;
    double x,y, x_m, y_m;
    for (int i=1; i<32; i++){
        move_particles(0.1);
        m1.update_COM();
        x = y = fmod(1.+i/10.,pi2);
        x_m = y_m = 1+i/10.;
        // Testing update_com()
        EXPECT_NEAR(x, m1.get_COM()[0], 1e-14);
        EXPECT_NEAR(y, m1.get_COM()[1], 1e-14);
        EXPECT_NEAR(x_m, m1.moved_COM()[0], 1e-14);
        EXPECT_NEAR(y_m, m1.moved_COM()[1], 1e-14);
    }
}

TEST_F(MoleculeTest, CopyConstructor){
    Molecule m2{m1};
    move_particles();
    EXPECT_EQ(Vector({1.1, 1.1}), m2.get_COM());
    m1.id = 1;
    m2.id = 2;
    EXPECT_EQ(1,m1.id);
    EXPECT_EQ(2,m2.id);
    m2.atoms = std::vector<Particle *>{};
    EXPECT_EQ(3, m1.num_particles());
    EXPECT_EQ(0, m2.num_particles());
}

TEST(Molecule, DeleteNeighbours){
    /* TODO
     *
     * delete_neighbours()
     * delete_mol_neighbours()
     */

}

TEST(Molecule, Getters){
    /* TODO
     *
     * uniqe_contacts()
     * num_particles()
     * num_contacts()
     * num_neighbours()
     * get_mass()
     * get_large()
     * get_COM()
     * atom_pos(int)
     * max_pairing()
     * index()
     * get_rotation()
     * get_orient_vect()
     * get_neighbours()
     */
}


TEST(Molecule, BoolOperators){
    /* TODO
     *
     * ==
     * !=
     * <
     * <=
     * >
     * >=
     */
}



#endif //ANALYSIS_TESTMOLECULE_H
