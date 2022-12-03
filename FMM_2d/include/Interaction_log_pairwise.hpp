#pragma once

#include <cmath>

#include "Interaction_base.hpp"

class Interaction_log_pairwise: public Interaction_base{
private:

public:
    double factor;

    Interaction_log_pairwise (double factor);
    double calc_energy_pcl_pcl(Particle *pcl_1, Particle *pcl_2);
    Vec_2d calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2);
    double calc_energy (std::vector<Particle *> &gas);
    void   calc_force  (std::vector<Particle *> &gas);
    double calc        (std::vector<Particle *> &gas);
};
