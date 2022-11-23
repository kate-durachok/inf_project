#include "../include/Interaction_repultion.hpp"

///collision_table

std::vector<Particle*>& collision_table::cell(Particle *pcl){
    Vec_2d rel_pos = pcl->pos - pos;
    return arr[ (int)(rel_pos.x / cell_size) + (int)(rel_pos.y / cell_size) * size_x ];
}

std::vector<Particle*>& collision_table::cell(unsigned int x, unsigned int y){
    return arr[x + y * size_x];
}

collision_table::collision_table(std::vector<Particle *> &gas, double cell_size, double pcl_diam):cell_size(cell_size), pcl_diam(pcl_diam) {
    Vec_2d pos, dim;
    get_bounding_box(gas, pos, dim);
    this->pos = pos - dim;
    size_x = (int)(2*dim.x/cell_size) + 1;
    size_y = (int)(2*dim.y/cell_size) + 1;
    arr = std::vector< std::vector<Particle*> > (size_x * size_y);

    for (auto pcl : gas){
        cell(pcl).push_back(pcl);
    }
}

bool collision_table::check_collision(Particle *pcl_1, Particle *pcl_2){
    return (pcl_1->pos - pcl_2->pos).sqr() < pcl_diam * pcl_diam;
}

std::vector<std::pair<Particle *, Particle *> > collision_table::get_collisions(){
    const size_t neighbour_count = 4;
    smol_vec begin     [neighbour_count] = { { 0,  1}, { 0,  0}, { 0,  0}, { 0,  0} };
    smol_vec end       [neighbour_count] = { {-1,  0}, {-1,  0}, {-1, -1}, { 0, -1} };
    smol_vec neighbour [neighbour_count] = { { 1, -1}, { 1,  0}, { 1,  1}, { 0,  1} };
    std::vector<std::pair<Particle *, Particle *> > collisions;
    for (size_t i_limit = 0; i_limit < neighbour_count; ++i_limit){
        for (size_t i_y = begin[i_limit].y; i_y < size_y + end[i_limit].y; ++i_y){
            for (size_t i_x = begin[i_limit].x; i_x < size_x + end[i_limit].x; ++i_x){
                for (auto pcl_1 : cell(i_x, i_y)){
                    for (auto pcl_2 : cell(i_x+neighbour[i_limit].x, i_y+neighbour[i_limit].y)){
                        if ( check_collision(pcl_1, pcl_2) ){
                            collisions.push_back(std::make_pair(pcl_1, pcl_2));
                        }
                    }
                }
            }
        }
    }
    for (size_t i_y = 0; i_y < size_y; ++i_y){
        for (size_t i_x = 0; i_x < size_x; ++i_x){
            for (auto it_1 = cell(i_x, i_y).begin(); it_1 != cell(i_x, i_y).end(); ++it_1){
                for (auto it_2 = std::next(it_1); it_2 != cell(i_x, i_y).end(); ++it_2){
                    if ( check_collision(*it_1, *it_2) ){
                        collisions.push_back(std::make_pair(*it_1, *it_2));
                    }
                }
            }
        }
    }
    return collisions;
}

std::ostream& operator<<(std::ostream& os, collision_table &rha){
    os << "coll_table: " << rha.size_x << " " << rha.size_y << " " << rha.pos << " " << rha.cell_size << ":\n";
    for (size_t i_y=0; i_y<rha.size_y; ++i_y){
        for (size_t i_x=0; i_x<rha.size_x; ++i_x){
            os << i_x << " " << i_y << ":\n";
            for (size_t i=0; i<rha.cell(i_x, i_y).size(); ++i){
                os << *(rha.cell(i_x, i_y)[i]);
            }
        }
    }
    return os;
}


///Interaction_repulsion

Interaction_repulsion::Interaction_repulsion(double r_max, double factor):
     r_max(r_max), factor(factor), degree_a(2), degree_b(2) { }
Interaction_repulsion::Interaction_repulsion(double r_max, double factor, unsigned int degree_a, unsigned int degree_b):
     r_max(r_max), factor(factor), degree_a (degree_a), degree_b(degree_b) { }

double Interaction_repulsion::calc_energy_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    Vec_2d r = pcl_1->pos - pcl_2->pos;
    double r_max_sqr = r_max * r_max;
    if (r.sqr() > r_max_sqr){
        return 0;
    }
    return factor * std::pow(std::pow(r_max_sqr/r.sqr(), degree_a) - 1, degree_b);
}

Vec_2d Interaction_repulsion::calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    Vec_2d r = pcl_2->pos - pcl_1->pos;
    double r_max_sqr = r_max * r_max;
    if (r.sqr() > r_max_sqr){
        return Vec_2d(0, 0);
    }
    return (-factor * 2.0*degree_a*degree_b * std::pow(std::pow(r_max_sqr/r.sqr(), degree_a) - 1, degree_b-1) * std::pow(r_max_sqr/r.sqr(), (degree_a+1)) / r_max_sqr) * r;
}

void Interaction_repulsion::calc_force(std::vector<Particle *> &gas){
    collision_table table(gas, r_max, r_max);
    std::vector<std::pair<Particle *, Particle *> > collisions = table.get_collisions();
    for (auto pair : collisions) {
        Vec_2d force = calc_force_pcl_pcl (pair.first, pair.second);
        pair.first->force  += force;
        pair.second->force -= force;
    }
}

double Interaction_repulsion::calc_energy(std::vector<Particle *> &gas){
    double energy = 0;
    collision_table table(gas, r_max, r_max);
    std::vector<std::pair<Particle *, Particle *> > collisions = table.get_collisions();
    for (auto pair : collisions) {
        energy += calc_energy_pcl_pcl(pair.first, pair.second);
    }
    return energy;
}

double Interaction_repulsion::calc(std::vector<Particle *> &gas) {
    double energy = 0;
    collision_table table(gas, r_max, r_max);
    std::vector<std::pair<Particle *, Particle *> > collisions = table.get_collisions();
    for (auto pair : collisions) {
        energy += calc_energy_pcl_pcl(pair.first, pair.second);
        Vec_2d force = calc_force_pcl_pcl (pair.first, pair.second);
        pair.first->force  += force;
        pair.second->force -= force;
    }
    return energy;
}

