//
// Created by malcolm on 19/02/16.
//

#ifndef ANALYSIS_TESTSEARCH_H
#define ANALYSIS_TESTSEARCH_H

#include <gtest/gtest.h>
#include "../src/Search.h"
#include "../src/Particle.h"
#include "../src/Molecule.h"

class modelMol {
public:
    int id;
    std::vector<modelMol *> my_neighbours;
    modelMol() : my_neighbours(std::vector<modelMol *>{}), id(0) {};

    void add_neighbour(modelMol* m){
        my_neighbours.push_back(m);
    }
    std::vector<modelMol *> get_neighbours() { return std::vector<modelMol*>{my_neighbours}; };
};

template <typename T>
class SearchTest : public testing::Test {
public:
    std::vector<T> list;
    std::vector<T> list_recip;

    virtual void SetUp() {
        list = std::vector<T>(13,T{});
        for (auto i=0; i != list.size(); i++) { list[i].id = i; };
        list[0].add_neighbour(&list[1]); list[0].add_neighbour(&list[2]); list[0].add_neighbour(&list[3]); list[0].add_neighbour(&list[4]);
        list[1].add_neighbour(&list[5]);
        list[2].add_neighbour(&list[8]); list[2].add_neighbour(&list[9]); list[2].add_neighbour(&list[10]); list[2].add_neighbour(&list[11]);
        list[4].add_neighbour(&list[12]);
        list[5].add_neighbour(&list[6]);
        list[6].add_neighbour(&list[7]);

        list_recip = std::vector<T>(13,T{});
        for (auto i=0; i != list_recip.size(); i++) { list_recip[i].id = i; };
        list_recip[0].add_neighbour(&list_recip[1]); list_recip[0].add_neighbour(&list_recip[2]); list_recip[0].add_neighbour(&list_recip[3]); list_recip[0].add_neighbour(&list_recip[4]);
        list_recip[1].add_neighbour(&list_recip[5]);
        list_recip[2].add_neighbour(&list_recip[8]); list_recip[2].add_neighbour(&list_recip[9]); list_recip[2].add_neighbour(&list_recip[10]); list_recip[2].add_neighbour(&list_recip[11]);
        list_recip[4].add_neighbour(&list_recip[12]);
        list_recip[5].add_neighbour(&list_recip[6]);
        list_recip[6].add_neighbour(&list_recip[7]);

        list_recip[1].add_neighbour(&list_recip[0]); list_recip[2].add_neighbour(&list_recip[0]); list_recip[3].add_neighbour(&list_recip[0]); list_recip[4].add_neighbour(&list_recip[0]);
        list_recip[5].add_neighbour(&list_recip[1]);
        list_recip[8].add_neighbour(&list_recip[2]); list_recip[9].add_neighbour(&list_recip[2]); list_recip[10].add_neighbour(&list_recip[2]); list_recip[11].add_neighbour(&list_recip[2]);
        list_recip[12].add_neighbour(&list_recip[4]);
        list_recip[6].add_neighbour(&list_recip[5]);
        list_recip[7].add_neighbour(&list_recip[6]);
    }

    void check_BFS_Order(std::vector<int> &actual){
        std::vector<int> expected{0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 6, 7};
        for (auto i=0; i< expected.size(); i++){
            EXPECT_EQ(expected[i], actual[i]);
        }
    }

    void check_DFS_Order(std::vector<int> &actual){
        std::vector<int> expected{0, 4, 12, 3, 2, 11, 10, 9, 8, 1, 5, 6, 7};
        for (auto i=0; i< expected.size(); i++){
            EXPECT_EQ(expected[i], actual[i]);
        }
    }

    void check_Order(Search<T> &search){
        std::vector<int> actual{};
        while (!search.empty()){
            actual.push_back(search.front()->id);
            search.pop();
        }
    }
};

typedef testing::Types<modelMol, Molecule, Particle> SearchClasses;

TYPED_TEST_CASE(SearchTest, SearchClasses);

TYPED_TEST(SearchTest, BFSSingle){
    BFS<TypeParam> search = BFS<TypeParam>(this->list.front());
    std::vector<int> actual{};
    while (!search.empty()){
        actual.push_back(search.front()->id);
        search.pop();
    }
    this->check_BFS_Order(actual);
}

TYPED_TEST(SearchTest, BFSRecip){
    BFS<TypeParam> search = BFS<TypeParam>(this->list_recip.front());
    std::vector<int> actual{};
    while (!search.empty()){
        actual.push_back(search.front()->id);
        search.pop();
    }
    this->check_BFS_Order(actual);
}

TYPED_TEST(SearchTest, DFSSingle){
    DFS<TypeParam> search = DFS<TypeParam>(this->list.front());
    std::vector<int> actual{};
    while (!search.empty()){
        actual.push_back(search.front()->id);
        search.pop();
    }
    this->check_DFS_Order(actual);
}

TYPED_TEST(SearchTest, DFSRecip){
    DFS<TypeParam> search = DFS<TypeParam>(this->list_recip.front());
    std::vector<int> actual{};
    while (!search.empty()){
        actual.push_back(search.front()->id);
        search.pop();
    }
    this->check_DFS_Order(actual);
}
/*
TYPED_TEST(SearchTest, InterfcaceSingle){
    Search<TypeParam> search = DFS<TypeParam>(this->list.front());
    std::vector<int> actual{};
    while (!search.empty()){
        actual.push_back(search.front()->id);
        search.pop();
    }
    this->check_DFS_Order(actual);

    search = BFS<TypeParam>(this->list.front());
    actual = std::vector<int>{};
    while (!search.empty()){
        actual.push_back(search.front()->id);
        search.pop();
    }
    this->check_DFS_Order(actual);
}

TYPED_TEST(SearchTest, InterfcaceRecip){
    Search<TypeParam> search = DFS<TypeParam>(this->list_recip.front());
    std::vector<int> actual{};
    while (!search.empty()){
        actual.push_back(search.front()->id);
        search.pop();
    }
    this->check_DFS_Order(actual);

    search = BFS<TypeParam>(this->list_recip.front());
    actual = std::vector<int>{};
    while (!search.empty()){
        actual.push_back(search.front()->id);
        search.pop();
    }
    this->check_DFS_Order(actual);
}
*/
#endif //ANALYSIS_TESTSEARCH_H
