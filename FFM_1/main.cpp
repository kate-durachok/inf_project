#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>
#include <set>
#include <complex>
#include <climits>
#include <string>
#include <fstream>

#include "vec.hpp"
#include "Particle.hpp"
#include "interaction_basic.hpp"
#include "gas.hpp"

double log_interaction_coef = 1;

unsigned long long f_rand (){
    static unsigned long long seed = 100;
    static unsigned long long prev_random = seed;
    prev_random ^= prev_random << 21;
    prev_random ^= prev_random >> 35;
    prev_random ^= prev_random << 4;
    return prev_random;
}

double f_rand_double(double min, double max){
    return (1.0 * f_rand()/ULLONG_MAX) * (max-min) + min;
}

unsigned long long get_nanoseconds(){
    static auto t_0 = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_0).count();
}

double potential_energy_pcl_pcl(Particle *pcl_1, Particle *pcl_2){
    return -log_interaction_coef * ((pcl_1->charge * pcl_2->charge) / 2) * std::log((pcl_1->pos - pcl_2->pos).sqr());
}

Vec_2d calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    return log_interaction_coef * ((pcl_1->charge * pcl_2->charge) / (pcl_1->pos - pcl_2->pos).sqr()) * (pcl_1->pos - pcl_2->pos);
}

std::complex<double> vec_to_complex(Vec_2d vec){
    return std::complex<double> (vec.x, vec.y);
}

std::vector< std::complex<double> > get_powers(std::complex<double> z, unsigned int pow){
    std::vector< std::complex<double> > powers(pow+1);
    std::complex<double> z_pow = 1.0;
    for (size_t i=0; i<=pow; ++i){
        powers[i] = z_pow;
        z_pow*=z;
    }
    return powers;
}

class binomial_coef_LUT{
private:
    double *table;
    unsigned int size;
public:

    double operator()(unsigned int n, unsigned int k){
        return table[n + k * size];
    }
    binomial_coef_LUT(unsigned int size): size(size){
        table = new double [size * size];
        for (unsigned int n=0; n<size; ++n){
            table[n + 0 * size] = 1;
        }
        for (unsigned int k=1; k<size; ++k){
            for (unsigned int n=k; n<size; ++n){
                table[n + k * size] = table[n + (k-1) * size] * (n-k+1)/(k);
            }
        }
    }

    ~binomial_coef_LUT(){
         delete[] table;
    };
    binomial_coef_LUT (binomial_coef_LUT const &src): size(src.size){
        table = new double [size * size];
        for (unsigned int i=0; i < size * size; ++i){
            table[i] = src.table[i];
        }
    };
    binomial_coef_LUT& operator=(binomial_coef_LUT const &src){
        if (this == &src){
            return *this;
        }
        binomial_coef_LUT tmp(src);
        std::swap(this->size,  tmp.size);
        std::swap(this->table, tmp.table);
        return *this;
    };
} binomial_coef(100 * 2);

struct tree_coordinates{
    int x, y;
    int level;

    tree_coordinates() = default;
    tree_coordinates(int x, int y, int level):x(x), y(y), level(level){}

    unsigned int upper_area(){
        return ((1<<(level*2))-1)/3;
    }
    unsigned int index(){
        return x + y*(1<<level) + upper_area();
    }
    tree_coordinates son(unsigned int i_x, unsigned int i_y){
        return tree_coordinates(2*x+i_x, 2*y+i_y, level+1);
    }
    tree_coordinates parent(){
        return tree_coordinates(x/2, y/2, level-1);
    }
    tree_coordinates neighbour(int i_x, int i_y){
        if ( (x+i_x>=0) && (x+i_x<(1<<level)) && (y+i_y>=0) && (y+i_y<(1<<level)) ){
            return tree_coordinates(x+i_x, y+i_y, level);
        }else{
            return tree_coordinates(0, 0, -1);
        }
    }
};

struct node{
    Vec_2d pos;
    double r;
    std::vector<Particle *> gas;

    std::vector< std::complex<double> > expansion_coef;
    std::vector< std::complex<double> > expansion_local;

    void calc_expansion(size_t expansion_degree){
        expansion_coef.resize(expansion_degree + 1);
        for (auto pcl : gas){
            expansion_coef[0] += pcl->charge;
            std::vector< std::complex<double> > z_power = get_powers(vec_to_complex (pcl->pos - pos), expansion_degree);
            for (size_t k=1; k<=expansion_degree; ++k){
                expansion_coef[k] += (double)pcl->charge * z_power[k];
            }
        }
        for (size_t k=1; k<=expansion_degree; ++k){
            expansion_coef[k] *= -1.0/k;
        }

    }

    void upward_translation(std::vector<node> &tree, tree_coordinates curr){
        expansion_coef = std::vector< std::complex<double> >(tree[curr.son(0, 0).index()].expansion_coef.size());
        for (int i_x=0; i_x<2; ++i_x){
            for (int i_y=0; i_y<2; ++i_y){
                node *node_p = &tree[curr.son(i_x, i_y).index()];

                expansion_coef[0] += node_p->expansion_coef[0];

                std::vector< std::complex<double> > z_power = get_powers(vec_to_complex(node_p->pos - pos), expansion_coef.size());
                for (size_t i=1; i<expansion_coef.size(); ++i){
                    expansion_coef[i] -= node_p->expansion_coef[0] * z_power[i] / ((double)i);
                    for (size_t k=1; k<=i; ++k){
                        expansion_coef[i] += node_p->expansion_coef[k] * z_power[i-k] * binomial_coef(i-1, k-1);
                    }
                }
            }
        }
    }

    node(){};
    node(Vec_2d pos, double r): pos(pos), r(r) {};
    node(Vec_2d pos, double r, std::vector<Particle *> &gas): pos(pos), r(r), gas(gas) {};
    node(Vec_2d pos, double r, std::vector<Particle *> &gas, size_t expansion_degree): pos(pos), r(r), gas(gas) {calc_expansion(expansion_degree);};

    friend std::ostream& operator<<(std::ostream& os, const node & rha){
        os << rha.pos << " " << rha.r << " " << rha.gas.size();
        os << "\n\n";
        for (auto coef : rha.expansion_coef){
            os << coef << " ";
        }
        os << "\n";
        for (auto coef : rha.expansion_local){
            os << coef << " ";
        }
        os << "\n";
        return os;
    };
};

void print_tree(std::vector<node> &tree, tree_coordinates curr, int max_level){
    if (curr.level > max_level){
        return;
    }
    for (int i=0; i<curr.level; ++i){
        std::cout << " ";
    }
    node *curr_p = &tree[curr.index()];
    std::cout << "|-" << *curr_p << "\n";

    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            print_tree(tree, curr.son(i_x, i_y), max_level);
        }
    }
}

void buildtree_recursion_translation(std::vector<node> &tree, tree_coordinates curr, size_t expansion_degree, int max_level){
    node *curr_p = &tree[curr.index()];
    if(curr.level == max_level){
        curr_p->calc_expansion(expansion_degree);
        return;
    }

    std::vector<Particle *> gas_son [2][2];
    for (auto pcl : curr_p->gas){
        int i_x=0, i_y=0;
        if (pcl->pos.x > curr_p->pos.x) { i_x=1; }
        if (pcl->pos.y > curr_p->pos.y) { i_y=1; }
        gas_son[i_x][i_y].push_back(pcl);
    }
    for (int i_x=0; i_x<2; ++i_x){
        for (int i_y=0; i_y<2; ++i_y){
            double r_new = curr_p->r / 2;
            Vec_2d delta_pos = Vec_2d(r_new * (2*i_x-1), r_new * (2*i_y-1));

            tree_coordinates son = curr.son(i_x, i_y);

            tree[son.index()] = node(curr_p->pos + delta_pos, r_new, gas_son[i_x][i_y]);
            buildtree_recursion_translation(tree, son, expansion_degree, max_level);
        }
    }
    curr_p->upward_translation(tree, curr);
}

std::vector<node> new_tree(std::vector<Particle *> &gas, size_t expansion_degree, int max_level){
    Vec_2d pos; double r;
    get_bounding_box(gas, pos, r);
    tree_coordinates end(0, 0, max_level+1);

    std::vector<node> new_tree( end.upper_area() );

    new_tree[0] = node(pos, r, gas);
    buildtree_recursion_translation(new_tree, tree_coordinates(0, 0, 0), expansion_degree, max_level);

    return new_tree;
}

std::vector<node *> get_interaction_list (std::vector<node> &tree, tree_coordinates curr){
    std::vector<node *> inter_list;
    tree_coordinates parent = curr.parent();
    for(int i_x=-1; i_x<=1; ++i_x){
        for(int i_y=-1; i_y<=1; ++i_y){
            tree_coordinates parent_neighbour = parent.neighbour(i_x, i_y);
            if (parent_neighbour.level != -1){
                for(int i_x_son=0; i_x_son<2; ++i_x_son){
                    for(int i_y_son=0; i_y_son<2; ++i_y_son){
                        tree_coordinates son = parent_neighbour.son(i_x_son, i_y_son);
                        if ( (son.x - curr.x > 1) || (son.x - curr.x < -1) || (son.y - curr.y > 1) || (son.y - curr.y < -1) ){
                            inter_list.push_back(&tree[son.index()]);
                        }
                    }
                }
            }
        }
    }
    return inter_list;
}

std::vector<node *> get_neighbour_list (std::vector<node> &tree, tree_coordinates curr){
    std::vector<node *> neighbour_list;
    for(int i_x=-1; i_x<=1; ++i_x){
        for(int i_y=-1; i_y<=1; ++i_y){
            tree_coordinates neighbour = curr.neighbour(i_x, i_y);
            if (neighbour.level != -1){
                neighbour_list.push_back(&tree[neighbour.index()]);
            }
        }
    }
    return neighbour_list;
}

std::vector< std::complex<double> > get_local_expansion(std::vector< std::complex<double> > &multipole_expansion, std::complex<double> z_0){
    std::vector< std::complex<double> > local_expansion (multipole_expansion.size());

    local_expansion[0] += multipole_expansion[0] * std::log(-z_0);

    std::vector< std::complex<double> > z_rev_power = get_powers((double)-1.0 / z_0, multipole_expansion.size() - 1);
    for (size_t k=1; k<multipole_expansion.size(); ++k){
        local_expansion[0] += multipole_expansion[k] * z_rev_power[k];
    }
    for (size_t l=1; l<multipole_expansion.size(); ++l){
        local_expansion[l] -= multipole_expansion[0] / ((double)l);
        for (size_t k=1; k<multipole_expansion.size(); ++k){
            local_expansion[l] += multipole_expansion[k] * z_rev_power[k] * binomial_coef(l+k-1, k-1);
        }
        local_expansion[l] *= z_rev_power[l] * ((double)1.0-2*(l%2));
    }

    return local_expansion;
}

void calc_local_samelevel(std::vector<node> &tree, tree_coordinates curr, int max_level){
    node *curr_p = &tree[curr.index()];
    std::vector<node *> inter_list = get_interaction_list(tree, curr);
    curr_p->expansion_local.resize(curr_p->expansion_coef.size());
    for (auto node_p : inter_list){
        std::complex<double> z_0 = vec_to_complex(node_p->pos - curr_p->pos);
        std::vector< std::complex<double> > node_p_local = get_local_expansion(node_p->expansion_coef, z_0);
        for(size_t i=0; i<curr_p->expansion_local.size(); ++i){
            curr_p->expansion_local[i] += node_p_local[i];
        }
    }
    if(curr.level == max_level){
        return;
    }
    for (int i_x=0; i_x<2; ++i_x){
        for (int i_y=0; i_y<2; ++i_y){
            calc_local_samelevel(tree, curr.son(i_x, i_y), max_level);
        }
    }
}

std::vector< std::complex<double> > translate_local_expansion(std::vector< std::complex<double> > &local_expansion, std::complex<double> z_0){
    std::vector< std::complex<double> > translated_expansion = local_expansion;
    size_t degree = translated_expansion.size()-1;
    for (size_t j=0; j<=degree-1; ++j){
        for (size_t k=degree-j-1; k<=degree-1; ++k){
            translated_expansion[k] -= z_0 * translated_expansion[k+1];
        }
    }
    return translated_expansion;
}

void calc_local_total(std::vector<node> &tree, tree_coordinates curr, int max_level){
    if(curr.level == max_level){
        return;
    }
    node *curr_p = &tree[curr.index()];
    for (int i_x=0; i_x<2; ++i_x){
        for (int i_y=0; i_y<2; ++i_y){
            tree_coordinates son = curr.son(i_x, i_y);
            node *son_p = &tree[son.index()];
            std::complex<double> z_0 = vec_to_complex(son_p->pos - curr_p->pos);
            std::vector< std::complex<double> > translated_expansion = translate_local_expansion(curr_p->expansion_local, (double)-1.0 * z_0);
            for (size_t i=0; i<curr_p->expansion_local.size(); ++i){
                son_p->expansion_local[i] += translated_expansion[i];
            }
            calc_local_total(tree, son, max_level);
        }
    }
}

double calc_potential_energy_pcl_node(Particle *pcl, node *node_p){
    std::complex<double> z = vec_to_complex(pcl->pos - node_p->pos);
    std::complex<double> potential_complex (0, 0);
    std::complex<double> z_pow = 1.0;
    for(size_t k=0; k<node_p->expansion_local.size(); ++k){
        potential_complex += node_p->expansion_local[k] * z_pow;
        z_pow *= z;
    }
    return -log_interaction_coef * pcl->charge * potential_complex.real();
}

double calc_potential_energy_FMM(std::vector<node> &tree, tree_coordinates curr, int max_level){
    double energy = 0;
    if (curr.level == max_level){
        node *curr_p = &tree[curr.index()];
        for (auto pcl : curr_p->gas){
            energy += calc_potential_energy_pcl_node(pcl, curr_p);
        }
        std::vector<node *> neighbour_list = get_neighbour_list(tree, curr);
        for (auto pcl_1 : curr_p->gas){
            for (auto node_p : neighbour_list){
                for (auto pcl_2 : node_p->gas){
                    if (pcl_1 != pcl_2){
                        energy += potential_energy_pcl_pcl(pcl_1, pcl_2);
                    }
                }
            }
        }
    }else{
        for (int i_x=0; i_x<2; ++i_x){
            for (int i_y=0; i_y<2; ++i_y){
                energy += calc_potential_energy_FMM(tree, curr.son(i_x, i_y), max_level);
            }
        }
    }
    return energy;
}

double calc_potential_energy_multipole(std::vector<Particle *> &gas, size_t expansion_degree){
    int max_level = 5;//std::log(gas.size() / (0.5 * expansion_degree) ) / std::log(4);
    if (max_level < 0) { max_level = 0; }
    //std::cout << "n = " << gas.size() << "   depth = " << max_level << "   pcl_per_leaf = " << 1.0 * gas.size() / (1<<(max_level*2)) << "\n";
    std::vector<node> tree = new_tree(gas, expansion_degree, max_level);
    calc_local_samelevel (tree, tree_coordinates(0, 0, 0), max_level);
    calc_local_total     (tree, tree_coordinates(0, 0, 0), max_level);
    //print_tree(tree, tree_coordinates(0, 0, 0), max_level);
    double energy = calc_potential_energy_FMM(tree, tree_coordinates(0, 0, 0), max_level) / 2;
    return energy;
}

Vec_2d calc_force_pcl_node(Particle *pcl, node *node_p){
    std::complex<double> z = vec_to_complex(pcl->pos - node_p->pos);
    std::complex<double> force_complex (0, 0);
    std::complex<double> z_pow = 1.0;
    for(size_t k=1; k<node_p->expansion_local.size(); ++k){
        force_complex += node_p->expansion_local[k] * (double)k * z_pow;
        z_pow *= z;
    }
    return log_interaction_coef * pcl->charge * Vec_2d(force_complex.real(), -force_complex.imag());
}

void calc_acceleration_FMM(std::vector<node> &tree, tree_coordinates curr, int max_level){
    if (curr.level == max_level){
        node *curr_p = &tree[curr.index()];
        for (auto pcl : curr_p->gas){
            pcl->force  += calc_force_pcl_node(pcl, curr_p);
        }
        std::vector<node *> neighbour_list = get_neighbour_list(tree, curr);
        for (auto pcl_1 : curr_p->gas){
            for (auto node_p : neighbour_list){
                for (auto pcl_2 : node_p->gas){
                    if (pcl_1 != pcl_2){
                        pcl_1->force  += calc_force_pcl_pcl(pcl_1, pcl_2);
                    }
                }
            }
        }
        for (auto pcl : curr_p->gas){
            pcl->acc /= pcl->mass;
        }
    }else{
        for (int i_x=0; i_x<2; ++i_x){
            for (int i_y=0; i_y<2; ++i_y){
                calc_acceleration_FMM(tree, curr.son(i_x, i_y), max_level);
            }
        }
    }
}

void calc_acceleration_multipole(std::vector<Particle *> &gas, size_t expansion_degree){
    int max_level = std::log(gas.size() / (0.5 * expansion_degree) ) / std::log(4);
    if (max_level < 0) { max_level = 0; }
    std::cout << "n = " << gas.size() << "   depth = " << max_level << "   pcl_per_leaf = " << 1.0 * gas.size() / (1<<(max_level*2)) << "\n";
    std::vector<node> tree = new_tree(gas, expansion_degree, max_level);
    calc_local_samelevel  (tree, tree_coordinates(0, 0, 0), max_level);
    calc_local_total      (tree, tree_coordinates(0, 0, 0), max_level);
    calc_acceleration_FMM (tree, tree_coordinates(0, 0, 0), max_level);
}

double calc_multipole(std::vector<Particle *> &gas, size_t expansion_degree){
    int max_level = std::log(gas.size() / (0.5 * expansion_degree) ) / std::log(4);
    if (max_level < 0) { max_level = 0; }
    //std::cout << "n = " << gas.size() << "   depth = " << max_level << "   pcl_per_leaf = " << 1.0 * gas.size() / (1<<(max_level*2)) << "\n";
    std::vector<node> tree = new_tree(gas, expansion_degree, max_level);
    calc_local_samelevel  (tree, tree_coordinates(0, 0, 0), max_level);
    calc_local_total      (tree, tree_coordinates(0, 0, 0), max_level);
    calc_acceleration_FMM (tree, tree_coordinates(0, 0, 0), max_level);
    double energy = calc_potential_energy_FMM(tree, tree_coordinates(0, 0, 0), max_level) / 2;
    return energy;
}

double calc_potential_energy_pairwise (std::vector<Particle *> &gas){
    double energy = 0;
    for (auto it_1 = gas.begin(); it_1 != gas.end(); ++it_1){
        Particle *pcl_1 = *it_1;
        for (auto it_2 = std::next(it_1); it_2 != gas.end(); ++it_2){
            Particle *pcl_2 = *it_2;
            energy += potential_energy_pcl_pcl(pcl_1, pcl_2);
        }
    }
    return energy;
}

void calc_force_pairwise (std::vector<Particle *> &gas){
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

double verlet_init(std::vector<Particle *> &gas){
    size_t expansion_degree = 20;
    double energy = calc_multipole(gas, expansion_degree);
    return energy;
}

double verlet(std::vector<Particle *> &gas, double dt){
    for(auto pcl : gas){
        pcl->pos += pcl->vel * dt + pcl->acc * (dt*dt/2);
        pcl->acc_prev = pcl->acc;
    }
    size_t expansion_degree = 10; ///

    Interaction_logarithmic_pairwise inter_1(-1.0);
    Interaction_repulsion inter_2(1.0, 0.1, 1, 2);
    Interaction_trap      inter_3(Vec_2d(0, 0), 60.0, 1.0, 10);

    reset_force(gas);
    //for (auto pcl:gas){ pcl->force -= pcl->vel * 0.05; }

    double energy = 0;
    unsigned long long t_0 = get_nanoseconds();
    //energy += inter_1.calc(gas);
    energy += calc_multipole(gas, expansion_degree);
    unsigned long long t_1 = get_nanoseconds();
    energy += inter_2.calc(gas);
    unsigned long long t_2 = get_nanoseconds();
    energy += inter_3.calc(gas);
    unsigned long long t_3 = get_nanoseconds();
//    printf("-----------------\n");
//    printf("FMM:   %10llu\n", (t_1-t_0)/1000);
//    printf("rep:   %10llu\n", (t_2-t_1)/1000);
//    printf("trap:  %10llu\n", (t_3-t_2)/1000);

    force_to_acc(gas);
    for(auto pcl : gas){
        pcl->vel += (pcl->acc + pcl->acc_prev) * (dt/2);
    }
    return energy;
}

int main()
{
    std::cout.precision(3);
    std::cout << std::fixed;
    int n = 1<<10;
    std::vector<Particle *> gas;

//    gas.push_back(new Particle(Vec_2d(0.0, 0), Vec_2d(  0, 0), 1, 1.0));
//    gas.push_back(new Particle(Vec_2d(0.5, 0), Vec_2d(-10, 0), 1, 1.0));
    {
        double dist = 1.0;
        int width = 10;
        for (int i=0; i<n/4; ++i){
            for (int k=0; k<4; ++k){
                int x_side = (2*(int)(k/2%2)-1);
                int y_side = (2*(int)(k%2)-1);
                Vec_2d pos(dist*(i/width+0.5) * x_side, dist*(i%width+0.5 ) * y_side);
                Vec_2d vel(0.0, 0.0);
                double mass = 1.0;
                double charge = 0.05 * y_side;
                gas.push_back(new Particle(pos, vel, mass, charge));
            }
        }
    }

    print_frame(gas, 0);
    verlet_init(gas);
    double dt = 1.0;
    int subframe_amm = 20;
    for (int i=1; i<100; ++i){
        double potential = 0;
        for (int k=0; k<subframe_amm; ++k){
            potential = verlet(gas, dt/subframe_amm);
        }
        std::cout << "iter_n:" << i << "\n";
        std::cout << calc_kinetic_energy(gas) + potential << "\n";
        std::cout << "------------" << "\n";
        print_frame(gas, i);
    }

    for (auto pcl : gas){
        delete pcl;
    }

    return 0;
}
