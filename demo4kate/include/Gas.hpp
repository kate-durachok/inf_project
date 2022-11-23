#pragma once

#include <iostream>
#include <fstream>
#include <vector>

#include "Particle.hpp"
#include "Interaction_base.hpp"

void get_bounding_box (std::vector<Particle *> &gas, Vec_2d &pos, Vec_2d &dim);
void get_bounding_box (std::vector<Particle *> &gas, Vec_2d &pos, double &r);

void reset_force  (std::vector<Particle *> &gas);
void force_to_acc (std::vector<Particle *> &gas);

///std::vector<Particle *> read_gas (std::istream& in);
///std::vector<Particle *> read_gas (const std::string &filename);
void print_gas   (std::vector<Particle *> &gas, std::ostream& os);
void print_frame (std::vector<Particle *> &gas, int i);

double calc_kinetic_energy   (std::vector<Particle *> &gas);
double calc_potential_energy (std::vector<Particle *> &gas, std::vector<Interaction_base *> &inter_arr);
void verlet (std::vector<Particle *> &gas, std::vector<Interaction_base *> &inter_arr, double dt);
