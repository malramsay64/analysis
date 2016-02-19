//
// Created by malcolm on 28/01/16.
//

#ifndef ANALYSIS_COMPUTE_H
#define ANALYSIS_COMPUTE_H

#include "Frame.h"
#include <complex>

/*
 * Compute
 *
 * This class forms the basis for a number of compute functions.
 * The compute is performed on the const Frame with values being stored in an array as well as a cumulative statistic.
 *
 * The array is populated with data when the compute_array function is called.
 *
 * Eventually this file will be broken up into many smaller files to make everything easier to read, find and
 * understand. This will also include many of the functions that are either attached to the Particle, Molecule
 * or Frame classes to reduce their complexity.
 *
 */
class Compute {
protected:
    const Frame &frame;
public:
    std::vector<double> array;
    Stats::Stats stats{};

    Compute() : frame(Frame{}) {};
    Compute(const Frame &inFrame) : frame(inFrame) {};
    virtual void compute_array() {};
};


/* ComputeMol
 *
 * A subclass of compute which implements compute_array for molecules.
 * This initialises the size of array to the number of molecules in the frame, as well as an implementation
 * of the compute_array function which will be used by subclasses.
 *
 * The compute_array function makes called to the compute_single function, which computes the value in the
 * array for a single molecule. Each subclass of ComputeMol needs to have this implemented.
 *
 * The purpose of this class is to abstract the looping over all the molecules from calculating the value of
 * a single molecule. This allows for optimisation eg. in the form of multithreading, without having to
 * implement for every subclass.
 */
class ComputeMol : public Compute {
protected:
    virtual double compute_single(const Molecule &m) {};
public:
    ComputeMol(const Frame &inFrame) : Compute(inFrame), array(std::vector<double>(frame.num_mol(), 0.)){};

    void compute_array();
};


/*
 * ComputeAtom
 *
 * This is the complementary class to ComputeMol, through instead of computing values for every molecule,
 * the values are calculated for every particle. See ComputeMol for more information.
 */
class ComputeAtom : public Compute {
protected:
    virtual double compute_single(const Particle &m) {};
public:
    ComputeAtom(const Frame &inFrame) : Compute(inFrame), array(std::vector<double>(frame.num_atoms(), 0.)){};

    void compute_array();
};

/*
 * ComputeAngle
 *
 * Class to calculate the angle of every molecule. The angle is the value stored by each individual
 * molecule.
 */
class ComputeAngle : public ComputeMol {
public:
    ComputeAngle(const Frame &inFrame) : ComputeMol(inFrame) { }

    double compute_single(const Molecule &m);
};


/*
 * ComptuteCoord
 *
 * This class will calculate the coordination number of each particle.
 * There is currently no differentiation between particle types.
 * The default cutoff of 1.12 is the value at which the LJ potential reaches
 * the minimum in the potential well with sigma = 1
 */
class ComputeCoord : public ComputeAtom {
private:
    const double cutoff;
public:
    ComputeCoord::ComputeCoord(const Frame &inFrame, double cutoff = 1.12) : ComputeAtom(inFrame), cutoff(cutoff){}

    double compute_single(const Particle &);
};



class ComputeMolHexatic : public ComputeMol {
private:
    const int degree;
public:
    ComputeMolHexatic::ComputeMolHexatic(const Frame &inFrame, int d) : ComputeMol(inFrame), degree(d) { }

    double compute_single(const Molecule &m);
};

class ComputeCircleOrder: public ComputeAtom {
public:
    ComputeCircleOrder(const Frame &inFrame) : ComputeAtom(inFrame) {};

    double compute_single(const Particle &);
};

#endif //ANALYSIS_COMPUTE_H
