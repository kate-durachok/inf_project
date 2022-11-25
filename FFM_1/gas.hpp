#pragma once

#include <iostream>
#include <fstream>

#include "Particle.hpp"

void get_bounding_box(std::vector<Particle *> &gas, Vec_2d &pos, Vec_2d &dim){
    double min_x = INFINITY, min_y = INFINITY, max_x = -INFINITY, max_y = -INFINITY;
    for (auto pcl:gas){
        if (min_x > pcl->pos.x){ min_x = pcl->pos.x; }
        if (max_x < pcl->pos.x){ max_x = pcl->pos.x; }
        if (min_y > pcl->pos.y){ min_y = pcl->pos.y; }
        if (max_y < pcl->pos.y){ max_y = pcl->pos.y; }
    }
    pos = Vec_2d( (min_x + max_x) / 2, (min_y + max_y) / 2);
    dim = Vec_2d( (max_x - min_x) / 2, (max_y - min_y) / 2);
}

void get_bounding_box(std::vector<Particle *> &gas, Vec_2d &pos, double &r){
    Vec_2d dim;
    get_bounding_box(gas, pos, dim);
    r = dim.x;
    if (r < dim.y){
        r = dim.y;
    }
}

void reset_force (std::vector<Particle *> &gas){
    for(auto pcl : gas){
        pcl->force = Vec_2d(0, 0);
    }
}

void force_to_acc (std::vector<Particle *> &gas){
    for(auto pcl : gas){
        pcl->acc = pcl->force / pcl->mass;
    }
}

void print_gas(std::vector<Particle *> &gas, std::ostream& os){
    os << "x y vx vy charge mass\n";
    for(auto pcl : gas){
        pcl->print_minimal(os);
    }
}

void print_frame(std::vector<Particle *> &gas, int i){
    std::ofstream f_out;
    f_out.open ("data/frames/frame" + std::to_string(i) + ".txt");
    print_gas(gas, f_out);
    f_out.close();
}

double calc_kinetic_energy(std::vector<Particle *> &gas){
    double energy = 0;
    for(auto pcl : gas){
        energy += pcl->vel.sqr() * pcl->mass;
    }
    return energy / 2;
}
