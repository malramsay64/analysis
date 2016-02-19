//
// Created by malcolm on 28/01/16.
//

#include "Compute.h"

void ComputeMol::compute_array() {
    for (auto const &m: frame.molecules){
        double val = compute_single(m);
        array[m.index()] = val;
        stats.push(val);
    }
}

void ComputeAtom::compute_array() {
    for (auto const &p: frame.particles){
        double val = compute_single(p);
        array[p.index()] = val;
        stats.push(val);
    }
}

double ComputeAngle::compute_single(const Molecule &m){
    return m.get_orientation();
}

double ComputeMolHexatic::compute_single(const Molecule &m) {
    double theta;
    std::complex<double> sum = std::complex<double>(0,0);
    std::complex<double> i = std::complex<double>(0,1);
    for (auto &m2: m.my_neighbours){
        theta = direction(m2.first->get_COM(), m.get_COM(), frame).angular()[1];
        sum += (1./m.num_neighbours())*exp(6*theta*i);
    }
    return abs(sum);
}

double ComputeCoord::compute_single(const Particle &p) {
    // TODO actually check the distance against the cutoff distance
    /* TODO search further than just neighbours
     * this will require searching the neighbour list
     */
    return p.numn();
}