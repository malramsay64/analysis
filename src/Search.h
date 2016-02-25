//
// Created by malcolm on 19/02/16.
//

#ifndef ANALYSIS_SEARCH_H
#define ANALYSIS_SEARCH_H

#include <list>
#include <set>
#include <climits>
#include "Molecule.h"

template <typename T>
class Search {
protected:
    std::list<std::pair<T *, int>> q;
    std::set<T *> traversed;
    size_t max_depth;

    Search() : q(std::list<std::pair<T*, int>>{}), traversed(std::set<T*>{}), max_depth(INT_MAX) {};
    Search(size_t depth) :  q(std::list<std::pair<T*, int>>{}), traversed(std::set<T*>{}), max_depth(depth) {}
    Search(T *elem) : Search() { push(elem, 0); };
    Search(T *elem, size_t depth) : Search(depth) {}

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
        int depth{get_depth()};
        q.pop_front();
        if (depth < max_depth ) {
            for (auto i: f->get_neighbours()) {
                push(i, depth + 1);
            }
        }
        return f;
    };
};

template <typename T>
class BFS : public Search<T> {
public:
    BFS() : Search<T>() {};
    BFS(size_t depth) : Search<T>(depth) {}
    BFS(T* elem) : Search<T>() { push(elem, 0); };
    BFS(T* elem, size_t depth) : Search<T>(depth) { push(elem, 0); };

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
    DFS(size_t depth) : Search<T>(depth) {};
    DFS(T* elem) : Search<T>() { push(elem, 0); };
    DFS(T* elem, size_t depth) : Search<T>(depth) { push(elem, 0); };

    void push(T* elem, int depth) {
        if (this->traversed.count(elem)) return;
        this->q.push_front(std::pair<T*, int>(elem, depth));
        this->traversed.insert(elem);
    };
};

#endif //ANALYSIS_SEARCH_H
