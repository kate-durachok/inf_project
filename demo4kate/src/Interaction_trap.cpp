#include "../include/Interaction_trap.hpp"

Interaction_trap::Interaction_trap(Vec_2d pos, double r_trap, double factor, unsigned int degree):
     pos(pos), r_trap(r_trap), factor(factor), degree(degree) { }

void Interaction_trap::calc_force(std::vector<Particle *> &gas){
    for(auto pcl : gas){
        Vec_2d r = pcl->pos - pos;
        if (r.sqr() > r_trap * r_trap){
            pcl->force -= (factor * degree * std::pow(r.len() - r_trap, degree-1) / r.len()) * r;
        }
    }
}

double Interaction_trap::calc_energy(std::vector<Particle *> &gas){
    double energy = 0;
    for(auto pcl : gas){
        Vec_2d r = pcl->pos - pos;
        if (r.sqr() > r_trap * r_trap){
            energy += factor * std::pow(r.len() - r_trap, degree);
        }
    }
    return energy;
}

double Interaction_trap::calc(std::vector<Particle *> &gas){
    double energy = 0;
    for(auto pcl : gas){
        Vec_2d r = pcl->pos - pos;
        if (r.sqr() > r_trap * r_trap){
            double r_len = r.len();
            double r_pow = std::pow(r_len - r_trap, degree-1);
            pcl->force -= (factor * degree * r_pow / r_len) * r;
            energy += factor * r_pow * (r_len - r_trap);
        }
    }
    return energy;
}
