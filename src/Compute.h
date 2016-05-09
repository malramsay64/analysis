//
// Created by malcolm on 28/01/16.
//

#ifndef ANALYSIS_COMPUTE_H
#define ANALYSIS_COMPUTE_H

#include "Frame.h"
#include "Neighbour_List.h"
#include <complex>
#include <queue>
#include <set>

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
    Compute(const Frame &inFrame, size_t array_size): frame(inFrame), array(std::vector<double>(array_size, 0.)) {};
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
    ComputeMol(const Frame &inFrame) : Compute(inFrame, frame.num_mol()) {};

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
    ComputeAtom(const Frame &inFrame) : Compute(inFrame, frame.num_atoms()) {};

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
    ComputeCoord(const Frame &inFrame, double cutoff = 1.12) : ComputeAtom(inFrame), cutoff(cutoff){}

    double compute_single(const Particle &);
};

class ComputeMolHexatic : public ComputeMol {
private:
    const int degree;
public:
    ComputeMolHexatic(const Frame &inFrame, int d) : ComputeMol(inFrame), degree(d) { }

    double compute_single(const Molecule &m);
};

class ComputeCircleOrder: public ComputeAtom {
public:
    ComputeCircleOrder(const Frame &inFrame) : ComputeAtom(inFrame) {};

    double compute_single(const Particle &);
};

class ComputeOrientationOrder: public ComputeMol {
public:
    ComputeOrientationOrder(const Frame &inFrame) : ComputeMol(inFrame) {};

    double compute_single(const Molecule &);
};


class ComputeAtomDynamic: public ComputeAtom {
protected:
    Frame &movedFrame;
public:
    ComputeAtomDynamic(const Frame& initFrame, Frame& movedFrame) : ComputeAtom(initFrame), movedFrame(movedFrame) {};
    //ComputeAtomDynamic(const Frame& initFrame) : ComputeAtom(initFrame) { };

    virtual void setMoved(const Frame& frame) { movedFrame = frame; };
};

/* Computes the displacement between all particles between the inital frame
 * and the moved frame.
 */
class ComputeDisplacement: public ComputeAtomDynamic {
public:
    ComputeDisplacement(const Frame &initFrame, Frame& movedFrame) : ComputeAtomDynamic(initFrame, movedFrame) {};
    //ComputeDisplacement(const Frame &initFrame) : ComputeAtomDynamic(initFrame) {};

    double compute_single(const Particle&);
};

class ComputeOverlap: public ComputeAtomDynamic {
private:
    ComputeDisplacement c_disp;
    const double factor = 0.25;
    const double limit = 0.2;
    double overlap_func(double d){
        return std::exp(-(d*d)/factor);
    }
public:
    ComputeOverlap(const Frame &initFrame, Frame& movedFrame) : ComputeAtomDynamic(initFrame, movedFrame), c_disp(initFrame,movedFrame) {
        c_disp.compute_array();
    };

    void set_moved(const Frame& frame) {
        c_disp.setMoved(frame);
        c_disp.compute_array();
    }

    double compute_single(const Particle&);
};

class ComputeClustering : public ComputeAtomDynamic {
private:
    ComputeOverlap c_overlap;
    Neighbour_List neighs;
    int cluster = 0;
    std::set<int> traversed;

public:
    ComputeClustering(const Frame &initFrame, Frame& movedFrame, Neighbour_List &neighs) :
            ComputeAtomDynamic(initFrame, movedFrame), c_overlap(initFrame,movedFrame), neighs(neighs) {
        c_overlap.compute_array();
    };
    ComputeClustering(const Frame &initFrame, Frame& movedFrame) :
            ComputeAtomDynamic(initFrame, movedFrame), c_overlap(initFrame, movedFrame) {
        c_overlap.compute_array();
        neighs = Neighbour_List{initFrame};
    }

    void set_moved(const Frame& frame) {
        c_overlap.setMoved(frame);
        c_overlap.compute_array();
        cluster = 0;
    }

    double compute_single(const Particle&);
};

#endif //ANALYSIS_COMPUTE_H
