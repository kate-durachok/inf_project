#include "../include/Interaction_short.hpp"

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

double Interaction_repulsion::calc_energy(std::vector<Particle *> &gas){
    double energy = 0;
    collision_table table(gas, r_max, r_max);
    std::vector<std::pair<Particle *, Particle *> > collisions = table.get_collisions();
    for (auto pair : collisions) {
        energy += calc_energy_pcl_pcl(pair.first, pair.second);
    }
    return energy;
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


///smooth_func

smooth_func::smooth_func(double r_0, double r_1, unsigned int degree):
    r_0(r_0), r_1(r_1), degree(degree) {}

double smooth_func::slow_pow(double x, unsigned int n){
    double pow = 1.0;
    for (unsigned int i=0; i<n; ++i){
        pow *= x;
    }
    return pow;
}

double smooth_func::operator()(double x){
    double x_normed = (x-r_0)/(r_1-r_0);
    if(x_normed < 0){
        return 0;
    }
    if(x_normed > 1){
        return 1;
    }
    double x_pow_0 = slow_pow(x_normed,     degree);
    double x_pow_1 = slow_pow(1 - x_normed, degree);
    return x_pow_0 / (x_pow_0 + x_pow_1);
}

double smooth_func::derivative(double x){
    double x_normed = (x-r_0)/(r_1-r_0);
    if(x_normed < 0 || x_normed > 1){
        return 0;
    }
    double x_pow_0 = slow_pow(x_normed,     degree-1);
    double x_pow_1 = slow_pow(1 - x_normed, degree-1);
    return (degree*x_pow_0*x_pow_1) / (slow_pow(x_normed*x_pow_0 + (1-x_normed)*x_pow_1, 2) * (r_1 - r_0));
}


///Interaction_6_12_smoothed

Interaction_6_12_smoothed::Interaction_6_12_smoothed(double eps, double sigma, double r_max, double r_unmod, unsigned int degree):
    eps(eps), sigma(sigma), smooth_f(r_max, r_unmod, degree) { }

double Interaction_6_12_smoothed::calc_energy_pcl_pcl_unmod (Particle *pcl_1, Particle *pcl_2){
    double dist_sqr = (pcl_2->pos - pcl_1->pos).sqr();

    double sigma_div_dist_2  = (sigma*sigma) / dist_sqr;
    double sigma_div_dist_6  = sigma_div_dist_2 * sigma_div_dist_2 * sigma_div_dist_2;
    double sigma_div_dist_12 = sigma_div_dist_6 * sigma_div_dist_6;

    return 4 * eps * (sigma_div_dist_12 - sigma_div_dist_6);
}

Vec_2d Interaction_6_12_smoothed::calc_force_pcl_pcl_unmod (Particle *pcl_1, Particle *pcl_2){
    double dist_sqr = (pcl_2->pos - pcl_1->pos).sqr();
    double sigma_sqr = sigma*sigma;

    double sigma_div_dist_2  = sigma_sqr/dist_sqr;
    double sigma_div_dist_4  = sigma_div_dist_2 * sigma_div_dist_2;
    double sigma_div_dist_8  = sigma_div_dist_4 * sigma_div_dist_4;
    double sigma_div_dist_14 = sigma_div_dist_2 * sigma_div_dist_4 * sigma_div_dist_8;

    return (-24*eps/sigma_sqr) * (2*sigma_div_dist_14 - sigma_div_dist_8) * (pcl_2->pos - pcl_1->pos);
}

double Interaction_6_12_smoothed::calc_energy_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    return calc_energy_pcl_pcl_unmod(pcl_1, pcl_2) * smooth_f((pcl_2->pos - pcl_1->pos).len());
}

Vec_2d Interaction_6_12_smoothed::calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    Vec_2d delta_pos = pcl_2->pos - pcl_1->pos;
    double dist = (delta_pos).len();
    return calc_force_pcl_pcl_unmod  (pcl_1, pcl_2) * smooth_f(dist) +
          (calc_energy_pcl_pcl_unmod (pcl_1, pcl_2) * smooth_f.derivative(dist)/dist) * delta_pos;

}

double Interaction_6_12_smoothed::calc_energy (std::vector<Particle *> &gas){
    double energy = 0;
    collision_table table(gas, smooth_f.r_0, smooth_f.r_0);
    std::vector<std::pair<Particle *, Particle *> > collisions = table.get_collisions();
    for (auto pair : collisions) {
        energy += calc_energy_pcl_pcl(pair.first, pair.second);
    }
    return energy;
}

void Interaction_6_12_smoothed::calc_force (std::vector<Particle *> &gas){
    collision_table table(gas, smooth_f.r_0, smooth_f.r_0);
    std::vector<std::pair<Particle *, Particle *> > collisions = table.get_collisions();
    for (auto pair : collisions) {
        Vec_2d force = calc_force_pcl_pcl (pair.first, pair.second);
        pair.first->force  += force;
        pair.second->force -= force;
    }
}

double Interaction_6_12_smoothed::calc (std::vector<Particle *> &gas){
    double energy = 0;
    collision_table table(gas, smooth_f.r_0, smooth_f.r_0);
    std::vector<std::pair<Particle *, Particle *> > collisions = table.get_collisions();
    for (auto pair : collisions) {
        energy      += calc_energy_pcl_pcl (pair.first, pair.second);
        Vec_2d force = calc_force_pcl_pcl (pair.first, pair.second);
        pair.first->force  += force;
        pair.second->force -= force;
    }
    return energy;
}
