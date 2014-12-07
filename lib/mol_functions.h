#include "particle.h"
#include "frame.h"
#include <list>

#ifndef MY_MOL_FUNCTIONS
#define MY_MOL_FUNCTIONS

template <class type> class dyn_queue{
    std::list<type *> q;
    std::list<int> depth;
    int push(type *, int);
  public:
    dyn_queue();
    dyn_queue(type *);
    type * pop();
    int remove(type *);
    int get_depth();
};    

class angle_list{
    std::vector<double> a;
    std::vector<int> count;
    std::vector<double> dist;
  public:
    angle_list();
    int push(double);
    int push(double, double);
    int print();
};

int get_colour(particle *p, Frame *frame);
int set_colour(Frame *frame);
int graph_colour(particle *p, Frame *frame);
double average_bonded_frac(particle *p, Frame * frame);
vect orientation(molecule *, Frame *);
double angle(molecule *, Frame *);
int colourise(Frame *frame);
angle_list like_me(Frame *frame);
double local_order(molecule * m, Frame * frame);

#endif
