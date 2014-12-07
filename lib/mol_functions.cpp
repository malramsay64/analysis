#include "mol_functions.h"

using namespace std;

static double deltaA = 5*PI/180;
static double deltaD = 0.15;

double average_bonded_frac(particle *p, Frame * frame){
    double mean = 0, stdev = 0, mprev, d, frac;
    vect v;
    vector<particle *>::iterator i;
    int num = 1;
    for (i = p->my_neighbours.begin(); i != p->my_neighbours.end(); i++){
        v = frame->direction(p->pos_vect(), (*i)->pos_vect());
        d = v.length();
        frac = d/(p->radius + (*i)->radius);
        mprev = mean;
        mean += (frac - mean)/(num);
        stdev += (frac-mean) * (frac - mprev);
        num++;
    }
    if (mean == 0){
        return 1;
    }
    return mean;
}

vect orientation(molecule *m, Frame *frame){
    vect com = m->COM();
    vect v = frame->direction(com, m->atoms.front()->pos_vect());
    if (m->nump() == 2 || (v.length() < EPS)){
         v = frame->direction(com, m->atoms.at(1)->pos_vect());
         v.orthogonalise();
     }
     v.normalise();
     return v;
}

double angle(molecule *m, Frame * frame){
    return atan2(orientation(m, frame)) + PI;
}

int get_colour(particle *p, Frame *frame){
    return frame->molecules.at(p->m_i()).get_colour();
}

int graph_colour(particle *p, Frame *frame){
    return frame->molecules.at(p->m_i()).graph_colour();
}

template <class type>
dyn_queue<type>::dyn_queue(){
    
}

template <class type>
dyn_queue<type>::dyn_queue(type * t){
    push(t,0);
}

template <class type>
int dyn_queue<type>::push(type * t, int d){
    if (!t->get_traversed()){
        t->traverse();
        q.push_back(t);
        depth.push_back(d);

    }
    return 0;
}

template <class type>
type * dyn_queue<type>::pop(){
    if (q.size()){
        type * t = q.front();
        q.pop_front();
        typename vector<type *>::iterator i;
        for (i = t->my_neighbours.begin(); i != t->my_neighbours.end(); i++){
            if ((*i)->type == t->type){
                push(*i, depth.front()+1);
                depth.pop_front();
            }
        }
        return t;
    }
    return 0;
}

template <class type>
int dyn_queue<type>::remove(type * t){
    q.remove(t);
    return 0;
}

template <class type>
int dyn_queue<type>::get_depth(){
    return depth.front();
}

angle_list::angle_list(){
    a = vector<double>(0,0); 
}


int angle_list::push(double d){
    for (int i = 0; i < a.size(); i++){
        if (abs(d-a.at(i)) < deltaA){
            a.at(i) += (d-a.at(i))/count.at(i);
            count.at(i)++;
            return i+1;
        }
    }
    a.push_back(d);
    count.push_back(1);
    return a.size();
}


int angle_list::push(double angle, double d){
    for (int i = 0; i < a.size(); i++){
        if (abs(angle-a.at(i)) < deltaA && abs(d-dist.at(i)) < deltaD){
            a.at(i) += (angle-a.at(i))/count.at(i);
            dist.at(i) += (d-dist.at(i))/count.at(i);
            count.at(i)++;
            return i+1;
        }
    }
    a.push_back(angle);
    dist.push_back(d);
    count.push_back(1);
    return a.size();
}

int angle_list::print(){
    if (dist.size() == a.size()){
        for (int i = 0; i < a.size(); i++){
            cout << a.at(i) << " " << dist.at(i) << endl;
        }
    }
    else{ 
        for (int i = 0; i < a.size(); i++){
            cout << a.at(i) << endl;
        }
    }
    return 0;
}

angle_list like_me(Frame *frame){
    vector<molecule>::iterator mol;
    angle_list angles;
    vect com1, com2, d;
    for (mol = frame->molecules.begin(); mol != frame->molecules.end(); mol++){
    dyn_queue<molecule> queue(&(*mol));
    queue.pop();
    molecule * m;
    com1 = (*mol).COM();
    while (queue.get_depth() < 4){
            m = queue.pop();
            if (m->get_colour() == (*mol).get_colour()){
                com2 = m->COM();
                d = direction(com1,com2);
                angles.push(atan2(d), d.length());          
            }
        }
    }
    return angles;
}


int colourise(Frame * frame){
    if (frame->get_coloured()){
        return 0;
    }
    angle_list a;
    vector<molecule>::iterator m;
    for (m = frame->molecules.begin(); m != frame->molecules.end(); m++){
        (*m).set_colour(a.push(angle(&(*m),frame)));
    }
    frame->set_coloured();
    return 0;
}

int graph_colourise(Frame *frame){
    // Frame already coloured
    if (frame->get_coloured()){
        return 0;
    }
    dyn_queue<particle> q = dyn_queue<particle>(&(frame->particles.front()));
    particle *p = q.pop();
    while (p){
        graph_colour(p, frame);
        p = q.pop();
    }
    frame->set_coloured();
    return 0;
}

double local_order(molecule * m, Frame * frame){
    vector<molecule *>::iterator i;
    double sum = 0;
    for (i = m->my_neighbours.begin(); i != m->my_neighbours.end(); i++){
        sum += pow(dot_product(orientation(m,frame),orientation((*i), frame)),2);
    }
    return sum/m->num_neighbours();
}
