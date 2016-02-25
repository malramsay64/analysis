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

double ComputeCircleOrder::compute_single(const Particle &p) {
    int small = 0, large = 0;
    for (auto &neigh: p.my_neighbours){
        if (neigh->type == 1){
            large++;
        }
        else{
            small++;
        }
    }
    // Ordering for d=1.637556 dimers
    // Not including bonded particle
    if (p.type == 1 && large == 3 && small == 3 ){
        return 1;
    }
    else if (p.type == 2 && large == 3 && small == 1 ){
        return 1;
    }
    return 0;
}

double ComputeOrientationOrder::compute_single(const Molecule &m) {
    Stats::Stats mean;
    for (auto m2: m.my_neighbours){
        mean.push(pow(dot_product(m.get_orient_vect(), m2.first->get_orient_vect()),2));
    }
    return mean.getMean();
}