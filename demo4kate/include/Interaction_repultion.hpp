#pragma once

#include <cmath>

#include "Interaction_base.hpp"
#include "Gas.hpp"

class collision_table{
private:
    struct smol_vec{
        int x, y;
    };
public:
    std::vector< std::vector<Particle*> > arr;
    size_t size_x, size_y;
    Vec_2d pos;
    double cell_size;
    double pcl_diam;

    std::vector<Particle*>& cell(Particle *pcl);
    std::vector<Particle*>& cell(unsigned int x, unsigned int y);
    collision_table(std::vector<Particle *> &gas, double cell_size, double pcl_diam);
    bool check_collision(Particle *pcl_1, Particle *pcl_2);
    std::vector<std::pair<Particle *, Particle *> > get_collisions();
    friend std::ostream& operator<<(std::ostream& os, const collision_table &rha);

};

struct Interaction_repulsion: public Interaction_base{
private:

public:
    double r_max;
    double factor;
    unsigned int degree_a; ///P(r) = (r.len() < r_max) * factor * (r^(-2*a) - 1)^b
    unsigned int degree_b; ///F(r) = (r.len() < r_max) * factor * ab*(r^(-2*a) - 1)^b-1 / r^(2*a+2) * r
    bool friction;
    double friction_k;

    Interaction_repulsion(double r_max, double factor);
    Interaction_repulsion(double r_max, double factor, unsigned int degree_a, unsigned int degree_b);
    double calc_energy_pcl_pcl (Particle *pcl_1, Particle *pcl_2);
    Vec_2d calc_force_pcl_pcl  (Particle *pcl_1, Particle *pcl_2);
    void calc_force   (std::vector<Particle *> &gas);
    double calc_energy(std::vector<Particle *> &gas);
    double calc       (std::vector<Particle *> &gas);
};
