#include "../include/Interaction_log_pairwise.hpp"

Interaction_log_pairwise::Interaction_log_pairwise (double factor): factor(factor) {}

double Interaction_log_pairwise::calc_energy_pcl_pcl(Particle *pcl_1, Particle *pcl_2){
    return -factor * ((pcl_1->charge * pcl_2->charge) / 2) * std::log((pcl_1->pos - pcl_2->pos).sqr());
}

Vec_2d Interaction_log_pairwise::calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    return factor * ((pcl_1->charge * pcl_2->charge) / (pcl_1->pos - pcl_2->pos).sqr()) * (pcl_1->pos - pcl_2->pos);
}

double Interaction_log_pairwise::calc_energy (std::vector<Particle *> &gas){
    double energy = 0;
    for (auto it_1 = gas.begin(); it_1 != gas.end(); ++it_1){
        Particle *pcl_1 = *it_1;
        for (auto it_2 = std::next(it_1); it_2 != gas.end(); ++it_2){
            Particle *pcl_2 = *it_2;
            energy += calc_energy_pcl_pcl(pcl_1, pcl_2);
        }
    }
    return energy;
}

void Interaction_log_pairwise::calc_force (std::vector<Particle *> &gas){
    for (auto it_1 = gas.begin(); it_1 != gas.end(); ++it_1){
        Particle *pcl_1 = *it_1;
        for (auto it_2 = std::next(it_1); it_2 != gas.end(); ++it_2){
            Particle *pcl_2 = *it_2;
            Vec_2d force = calc_force_pcl_pcl(pcl_1, pcl_2);
            pcl_1->force += force;
            pcl_2->force -= force;
        }
    }
}

double Interaction_log_pairwise::calc (std::vector<Particle *> &gas){
    calc_force(gas);
    return calc_energy(gas);
}
