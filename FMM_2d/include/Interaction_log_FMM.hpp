#pragma once

#include <cmath>

#include <complex>
#include <iomanip>
#include <vector>

#include "Interaction_base.hpp"
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

    static bool cout_detailed;
    static bool cout_paraview;
    friend std::ostream& operator<<(std::ostream& os, const Node & rha);
};

class Adaptive_tree{
private:
    void subdivide(Node *curr);
    void ctor_recursion (Node *curr);
    void dtor_recursion (Node *curr);

    void print_recursion (Node *curr, std::ostream& os, size_t depth);
public:
    const size_t gas_size_max;
    const size_t expansion_degree;
    Node *root;

    Adaptive_tree (std::vector<Particle *> &gas, size_t gas_size_max, size_t expansion_degree);
    ~Adaptive_tree ();
    Adaptive_tree (Adaptive_tree const &) = delete;
    Adaptive_tree& operator=(Adaptive_tree const &) = delete;

    friend std::ostream& operator<<(std::ostream& os, Adaptive_tree & rha);
};

class binomial_coef_LUT{
private:
    double *table;
    unsigned int size;
public:
    double operator()(unsigned int n, unsigned int k);
    binomial_coef_LUT(unsigned int size);

    ~binomial_coef_LUT();
    binomial_coef_LUT (binomial_coef_LUT const &src);
    binomial_coef_LUT& operator=(binomial_coef_LUT const &src);
};

class Interaction_log_FMM: public Interaction_base{
private:
    static std::complex<double> vec_to_complex(Vec_2d vec);
    static std::vector< std::complex<double> > get_powers(std::complex<double> z, unsigned int max_pow);
    static std::vector<Node *> get_interaction_list(Node *curr);
    static std::vector<Node *> get_neighbour_list(Node *curr);

    void make_global_expansion      (Node* node);
    void translate_global_expansion (Node* parent, Node* son);   ///son to parent
    void global_to_local_expansion  (Node* source, Node* target);
    void translate_local_expansion  (Node* parent, Node* son);   ///parent to son

    void make_global_expansion_leafs         (Node* curr);
    void translate_global_expansion_internal (Node* curr);
    void global_to_local_expansion_total     (Node* curr);
    void translate_local_expansion_total     (Node* curr);

    void complete_expansions(Node* root);

    double calc_energy_pcl_node (Particle *pcl_1, Node *node);
    double calc_energy_pcl_pcl  (Particle *pcl_1, Particle *pcl_2);
    Vec_2d calc_force_pcl_node  (Particle *pcl_1, Node *node);
    Vec_2d calc_force_pcl_pcl   (Particle *pcl_1, Particle *pcl_2);

    double calc_energy_local_expansion (Node* curr);
    double calc_energy_neighbours      (Node* curr);
    void   calc_forces_local_expansion (Node* curr);
    void   calc_forces_neighbours      (Node* curr);

    size_t expansion_degree;
    binomial_coef_LUT binomial_coef;
    //void print_nodes(Node *curr, int level);
public:
    double factor;
    size_t gas_size_max;

    void   set_expansion_degree(size_t degree);
    size_t get_expansion_degree();

    Interaction_log_FMM (double factor, size_t level_max, size_t expansion_degree);
    double calc_energy (std::vector<Particle *> &gas);
    void   calc_force  (std::vector<Particle *> &gas);
    double calc        (std::vector<Particle *> &gas);
};


void print_tree_paraview_recursion (Node *curr, std::ostream& os);
void print_tree_paraview  (Adaptive_tree &tree, int i);
