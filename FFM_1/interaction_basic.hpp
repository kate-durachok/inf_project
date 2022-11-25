#pragma once

#include <cmath>
#include <vector>

#include "gas.hpp"

class Interaction{
private:

public:
    virtual double calc_energy (std::vector<Particle *> &gas) = 0;
    virtual void   calc_force  (std::vector<Particle *> &gas) = 0;
    virtual double calc (std::vector<Particle *> &gas) = 0;

    virtual ~Interaction() = default;
};

class Interaction_empty: public Interaction{
private:

public:
    double calc_energy_pcl_pcl(Particle *pcl_1, Particle *pcl_2){
        double energy = 0;
        return energy;
    };
    Vec_2d calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
        Vec_2d force;
        return force;
    };
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

    ~Interaction_empty() = default;
};

class Interaction_logarithmic_pairwise: public Interaction{
private:

public:
    double factor;

    Interaction_logarithmic_pairwise (double factor): factor(factor) {}

    double calc_energy_pcl_pcl(Particle *pcl_1, Particle *pcl_2){
        return -factor * ((pcl_1->charge * pcl_2->charge) / 2) * std::log((pcl_1->pos - pcl_2->pos).sqr());
    }

    Vec_2d calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
        return factor * ((pcl_1->charge * pcl_2->charge) / (pcl_1->pos - pcl_2->pos).sqr()) * (pcl_1->pos - pcl_2->pos);
    }

    double calc_energy (std::vector<Particle *> &gas){
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

    void calc_force (std::vector<Particle *> &gas){
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

    double calc (std::vector<Particle *> &gas){
        double energy = 0;
        for (auto it_1 = gas.begin(); it_1 != gas.end(); ++it_1){
            Particle *pcl_1 = *it_1;
            for (auto it_2 = std::next(it_1); it_2 != gas.end(); ++it_2){
                Particle *pcl_2 = *it_2;
                energy += calc_energy_pcl_pcl(pcl_1, pcl_2);
                Vec_2d force = calc_force_pcl_pcl(pcl_1, pcl_2);
                pcl_1->force += force;
                pcl_2->force -= force;
            }
        }
        return energy;
    }

    ~Interaction_logarithmic_pairwise() = default;
};

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

    std::vector<Particle*>& cell(Particle *pcl){
        Vec_2d rel_pos = pcl->pos - pos;
        return arr[ (int)(rel_pos.x / cell_size) + (int)(rel_pos.y / cell_size) * size_x ];
    }

    std::vector<Particle*>& cell(unsigned int x, unsigned int y){
        return arr[x + y * size_x];
    }

    collision_table(std::vector<Particle *> &gas, double cell_size, double pcl_diam):cell_size(cell_size), pcl_diam(pcl_diam) {
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

    bool check_collision(Particle *pcl_1, Particle *pcl_2){
        return (pcl_1->pos - pcl_2->pos).sqr() < pcl_diam * pcl_diam;
    }

    std::vector<std::pair<Particle *, Particle *> > get_collisions(){
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

    friend std::ostream& operator<<(std::ostream& os, collision_table &rha){
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
};

struct Interaction_repulsion: public Interaction{
private:

public:
    double r_max;
    double factor;
    unsigned int degree_a; ///P(r) = (r.len() < r_max) * factor * (r^(-2*a) - 1)^b
    unsigned int degree_b; ///F(r) = (r.len() < r_max) * factor * ab*(r^(-2*a) - 1)^b-1 / r^(2*a+2) * r
    bool friction;
    double friction_k;

    Interaction_repulsion(double r_max, double factor):
         r_max(r_max), factor(factor), degree_a(2), degree_b(2) { }
    Interaction_repulsion(double r_max, double factor, unsigned int degree_a, unsigned int degree_b):
         r_max(r_max), factor(factor), degree_a (degree_a), degree_b(degree_b) { }

    double calc_energy_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
        Vec_2d r = pcl_1->pos - pcl_2->pos;
        double r_max_sqr = r_max * r_max;
        if (r.sqr() > r_max_sqr){
            return 0;
        }
        return factor * std::pow(std::pow(r_max_sqr/r.sqr(), degree_a) - 1, degree_b);
    }

    Vec_2d calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
        Vec_2d r = pcl_2->pos - pcl_1->pos;
        double r_max_sqr = r_max * r_max;
        if (r.sqr() > r_max_sqr){
            return Vec_2d(0, 0);
        }
        return (-factor * 2.0*degree_a*degree_b * std::pow(std::pow(r_max_sqr/r.sqr(), degree_a) - 1, degree_b-1) * std::pow(r_max_sqr/r.sqr(), (degree_a+1)) / r_max_sqr) * r;
    }

    void calc_force(std::vector<Particle *> &gas){
        collision_table table(gas, r_max, r_max);
        std::vector<std::pair<Particle *, Particle *> > collisions = table.get_collisions();
        for (auto pair : collisions) {
            Vec_2d force = calc_force_pcl_pcl (pair.first, pair.second);
            pair.first->force  += force;
            pair.second->force -= force;
        }
    }

    double calc_energy(std::vector<Particle *> &gas){
        double energy = 0;
        collision_table table(gas, r_max, r_max);
        std::vector<std::pair<Particle *, Particle *> > collisions = table.get_collisions();
        for (auto pair : collisions) {
            energy += calc_energy_pcl_pcl(pair.first, pair.second);
        }
        return energy;
    }

    double calc(std::vector<Particle *> &gas) {
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
};

class Interaction_trap: public Interaction{
private:
public:
    Vec_2d pos;
    double r_trap;         ///r = (pcl->pos-pos)
    double factor;         ///P(r) = (r.len() > r_trap) * factor * (r.len() - r_trap)^degree
    unsigned int degree;   ///F(r) = (r.len() > r_trap) * factor * degree * (r.len() - r_trap)^(degree-1)

    Interaction_trap(Vec_2d pos, double r_trap, double factor, unsigned int degree):
         pos(pos), r_trap(r_trap), factor(factor), degree(degree) { }

    void calc_force(std::vector<Particle *> &gas){
        for(auto pcl : gas){
            Vec_2d r = pcl->pos - pos;
            if (r.sqr() > r_trap * r_trap){
                pcl->force -= (factor * degree * std::pow(r.len() - r_trap, degree-1) / r.len()) * r;
            }
        }
    }

    double calc_energy(std::vector<Particle *> &gas){
        double energy = 0;
        for(auto pcl : gas){
            Vec_2d r = pcl->pos - pos;
            if (r.sqr() > r_trap * r_trap){
                energy += factor * std::pow(r.len() - r_trap, degree);
            }
        }
        return energy;
    }

    double calc(std::vector<Particle *> &gas){
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
};
