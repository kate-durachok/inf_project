#pragma once

#include <cmath>

#include "Interaction_base.hpp"

class Interaction_trap: public Interaction_base{
private:
public:
    Vec_2d pos;
    double r_trap;         ///r = (pcl->pos-pos)
    double factor;         ///P(r) = (r.len() > r_trap) * factor * (r.len() - r_trap)^degree
    unsigned int degree;   ///F(r) = (r.len() > r_trap) * factor * degree * (r.len() - r_trap)^(degree-1)

    Interaction_trap(Vec_2d pos, double r_trap, double factor, unsigned int degree);
    void calc_force   (std::vector<Particle *> &gas);
    double calc_energy(std::vector<Particle *> &gas);
    double calc       (std::vector<Particle *> &gas);
};
