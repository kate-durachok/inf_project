#pragma once

#include <vector>

#include "Particle.hpp"
#include "Stopwatch.hpp"

class Interaction_base{
private:

public:
    virtual double calc_energy (std::vector<Particle *> &gas) = 0;
    virtual void   calc_force  (std::vector<Particle *> &gas) = 0;
    virtual double calc        (std::vector<Particle *> &gas) = 0;

    virtual ~Interaction_base() = default;
};

class Interaction_empty: public Interaction_base{
private:

public:
    double calc_energy (std::vector<Particle *> &gas){
        double energy = 0;
        return energy;
    };
    void   calc_force (std::vector<Particle *> &gas){

    };
    double calc (std::vector<Particle *> &gas){
        double energy = 0;
        return energy;
    };
};
