//
// Created by malcolm on 19/02/16.
//

#ifndef ANALYSIS_SEARCH_H
#define ANALYSIS_SEARCH_H

#include <list>
#include <set>
#include "Molecule.h"

template <typename T>
class Search {
protected:
    std::list<std::pair<T *, int>> q;
    std::set<T *> traversed;

    Search() : q(std::list<std::pair<T*, int>>{}), traversed(std::set<T*>{}) {};
    Search(T *elem) : Search() { push(elem, 0); };
    Search(T &elem) : Search() { push(&elem, 0); };

public:
    virtual bool empty() { return q.empty(); };
    virtual std::size_t size() { return q.size(); };
    virtual T* front() { return q.front().first; };
    virtual int get_depth() { return q.front().second; };
    virtual void push(T* elem, int depth) = 0;
    virtual T* pop() {
        if (q.empty()){
            return nullptr;
        }
        T* f = front();
        q.pop_front();
        for (auto i: f->get_neighbours()){
            push(i, get_depth()+1);
        }
        return f;
    };
};

template <typename T>
class BFS : public Search<T> {
public:
    BFS() : Search<T>() {};
    BFS(T* elem) : Search<T>() { push(elem, 0); };
    BFS(T& elem) : Search<T>() { push(&elem, 0); };

    void push(T* elem, int depth) {
        if (this->traversed.count(elem)) return;
        this->q.push_back(std::pair<T*, int>(elem, depth));
        this->traversed.insert(elem);
    };
};

template <typename T>
class DFS : public Search<T> {
public:
    DFS() : Search<T>() {};
    DFS(T* elem) : Search<T>() { push(elem, 0); };
    DFS(T& elem) : Search<T>() { push(&elem, 0); };

    void push(T* elem, int depth) {
        if (this->traversed.count(elem)) return;
        this->q.push_front(std::pair<T*, int>(elem, depth));
        this->traversed.insert(elem);
    };
};

#endif //ANALYSIS_SEARCH_H
