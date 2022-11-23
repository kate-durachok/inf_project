#pragma once

#include <cmath>

#include "Interaction_base.hpp"

#include <complex>
#include <iomanip>
#include <vector>

#include "Particle.hpp"
#include "Gas.hpp"


class Node{
private:
    static size_t id_total;
public:
    Vec_2d pos;
    double r;

    std::vector<Particle *> gas;
    std::vector< std::complex<double> > expansion_global;
    std::vector< std::complex<double> > expansion_local;

    Node *parent;
    Node *son[2][2];
    Node *neighbour[3][3];

    bool leaf;
    const size_t id;

    static size_t get_id(const Node *curr);

    Node(Vec_2d pos, double r, Node *parent, std::vector<Particle *> &gas, size_t expansion_degree);

    friend std::ostream& operator<<(std::ostream& os, const Node & rha);
};

class adaptive_tree{
private:
    void subdivide(Node *curr);
    void ctor_recursion (Node *curr);
    void dtor_recursion (Node *curr);
    std::vector<Node *> get_interaction_list(Node *curr);
    std::vector<std::pair<Node *, Node *> > get_direct_list(Node *curr);

    void print_recursion (Node *curr, std::ostream& os, size_t depth);
public:
    const size_t gas_size_max;
    const size_t expansion_degree;
    Node *root;

    adaptive_tree (std::vector<Particle *> &gas, size_t gas_size_max, size_t expansion_degree);
    ~adaptive_tree ();

    friend std::ostream& operator<<(std::ostream& os, adaptive_tree & rha);
};
