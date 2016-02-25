#include "Frame.h"

Frame::Frame(){
    molecules = std::vector<Molecule>{};
    particles = std::vector<Particle>{};
    coloured = false;
    neighbours = false;
}

Frame::Frame(const Frame &frame){
    molecules = frame.molecules;
    particles = frame.particles;
    coloured = frame.coloured;
    neighbours = frame.neighbours;
    a = frame.a;
    b = frame.b;
    theta = frame.theta;
    timestep = frame.timestep;
    update_links();
}

Frame::Frame(const FrameVars &f) : Frame() {
    timestep = f.timestep;
    a = f.a;
    b = f.b;
    theta = f.theta;
    particles = std::vector<Particle>(f.num_atoms, Particle{});
};

void Frame::set_atoms(int a){
    particles.reserve(a);
}

size_t Frame::num_atoms() const{
    return particles.size();
}

void Frame::set_num_mol(int n){
    molecules.resize(n);
}

size_t Frame::num_mol() const{
    return molecules.size();
}

void Frame::add_particle(const Particle &p){
    particles.push_back(Particle{p});
}

void Frame::add_link(int molpos, int ppos){
    add_link(molpos, ppos, true);
}

void Frame::add_link(int molpos, int ppos, bool b_sort=true){
    Molecule *m = &molecules[molpos];
    m->atoms.push_back(&particles[ppos]);
    if (b_sort){
        sort(m->atoms.begin(), m->atoms.end());
    }
}

void Frame::update_links(){
    for (auto &m: molecules){
        // Links to particles
        for (auto &atom : m.atoms) {
            atom = &particles.at(atom->index());
        }
        // TODO Links to neighbours
    }
}

Molecule Frame::at(int i) const{
    return molecules.at(i);
}

void Frame::set_timestep(int t){
    timestep = t;
}

int Frame::get_timestep() const{
    return timestep;
}

double Frame::get_step_size() const {
    return step_size;
}

int Frame::get_time() const{
    return get_timestep()*get_step_size();
}

// Crystal Coordinates
void Frame::set_crys(double a, double b, double theta){
    this->a = a;
    this->b = b;
    this->theta = theta;
}

double Frame::get_a() const{
    return a;
}

double Frame::get_b() const{
    return b;
}

double Frame::get_theta() const{
    return theta;
}

double Frame::get_height() const{
    return b*sin(theta);
}

double Frame::get_tilt() const{
    return b*cos(theta);
}

double Frame::get_area() const{
    return get_a()*get_height();
}

double Frame::get_density() const{
    return num_atoms()/get_area();
}

double dist(const Vector &v1, const Vector &v2, const Frame &f){
    return direction(v1,v2,f).length();
}

Vector direction(const Vector &v1, const Vector &v2, const Frame &f){
    return cartesian(direction(v1,v2), f);
}

double angle(const Vector &v, const Frame &f){
    return atan2(cartesian(v,f));
}

Vector cartesian(const Vector &vi, const Frame &f){
    Vector v{vi};
    v /= (2*PI);
    v.r[0] = v.r[0]*f.get_a() + v.r[1]*f.get_b()*cos(f.get_theta());
    v.r[1] = v.r[1]*f.get_b()*sin(f.get_theta());
    return v;
}

Vector fractional(const Vector &vi, const Frame &f){
    Vector v{};
    v.r[0] = vi.r[0]*(1/f.get_a()) + vi.r[1]*(-cos(f.get_theta())/(f.get_a()*sin(f.get_theta())));
    v.r[1] = vi.r[1]/(f.get_b()*sin(f.get_theta()));
    v *= 2*PI;
    return v;
}

