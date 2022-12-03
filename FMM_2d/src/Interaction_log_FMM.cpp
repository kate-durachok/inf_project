#include "../include/Interaction_log_FMM.hpp"


///Node

size_t Node::id_total = 0;
bool Node::cout_paraview = true;
bool Node::cout_detailed = false;

size_t Node::get_id(const Node *curr){
    if (curr){
        return curr->id;
    }
    return 0;
}

Node::Node(Vec_2d pos, double r, Node *parent, std::vector<Particle *> &gas, size_t expansion_degree):
    pos(pos), r(r), gas(gas), expansion_global(expansion_degree+1), expansion_local(expansion_degree+1), parent(parent), son{{},{}}, neighbour{{},{},{}}, leaf(true), id(++id_total)
{
    neighbour[1][1] = this;
};

std::ostream& operator<<(std::ostream& os, const Node & rha){
    if (Node::cout_paraview){
        os << rha.pos.x << " " << rha.pos.y << " " << rha.r << " " << rha.gas.size() << "\n";
        return os;
    }
    size_t (*get_id)(const Node *) = Node::get_id;

    os << rha.pos << " " << rha.r << " " << rha.gas.size() << " ";
    if (Node::cout_detailed){
        os << "  " << get_id(&rha) << "  (" << get_id(rha.parent) << ", {" << get_id(rha.son[0][0]) << ", " << get_id(rha.son[1][0]) << ", " << get_id(rha.son[0][1]) << ", " << get_id(rha.son[1][1]) << "})";
    }
    os << "\n";
    if (Node::cout_detailed){
        const size_t num_w = 2;
        os <<  "                |" << std::setw(num_w) << get_id(rha.neighbour[0][0]) << ", " << std::setw(num_w) << get_id(rha.neighbour[1][0]) << ", " << std::setw(num_w) << get_id(rha.neighbour[2][0]) << "|\n";
        os <<  "                |" << std::setw(num_w) << get_id(rha.neighbour[0][1]) << ", " << std::setw(num_w) << get_id(rha.neighbour[1][1]) << ", " << std::setw(num_w) << get_id(rha.neighbour[2][1]) << "|\n";
        os <<  "                |" << std::setw(num_w) << get_id(rha.neighbour[0][2]) << ", " << std::setw(num_w) << get_id(rha.neighbour[1][2]) << ", " << std::setw(num_w) << get_id(rha.neighbour[2][2]) << "|\n";
    }
    return os;
}

///Tree

void Adaptive_tree::subdivide(Node *curr){
    curr->leaf = false;
    std::vector<Particle *> gas_son [2][2];
    for (auto pcl : curr->gas){
        int i_x=0, i_y=0;
        if (pcl->pos.x > curr->pos.x) { i_x=1; }
        if (pcl->pos.y > curr->pos.y) { i_y=1; }
        gas_son[i_x][i_y].push_back(pcl);
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            double new_r = curr->r / 2;
            Vec_2d new_pos = curr->pos + Vec_2d(new_r * (2*i_x-1), new_r * (2*i_y-1));

            curr->son[i_x][i_y] = new Node(new_pos, new_r, curr, gas_son[i_x][i_y], expansion_degree);

            for (int n_y=0; n_y<3; ++n_y){
                for (int n_x=0; n_x<3; ++n_x){
                    int n_x_curr = (i_x + n_x + 1)/2;
                    int n_y_curr = (i_y + n_y + 1)/2;
                    Node *curr_neighbour = curr->neighbour[n_x_curr][n_y_curr];
                    if(curr_neighbour){
                        int i_x_neig = (i_x + n_x + 1)%2;
                        int i_y_neig = (i_y + n_y + 1)%2;
                        Node *son_neighbour = curr_neighbour->son[i_x_neig][i_y_neig];
                        if (son_neighbour){
                            curr->son[i_x][i_y]->neighbour[n_x][n_y] = son_neighbour;
                            son_neighbour->neighbour[2-n_x][2-n_y] = curr->son[i_x][i_y];
                        }
                    }
                }
            }
        }
    }
}

void Adaptive_tree::ctor_recursion(Node *curr){
    if(curr->gas.size() <= gas_size_max){
        return;
    }
    subdivide(curr);
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            ctor_recursion(curr->son[i_x][i_y]);
        }
    }
}

void Adaptive_tree::dtor_recursion(Node *curr){
    if (curr == NULL){
        return;
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            dtor_recursion(curr->son[i_x][i_y]);
        }
    }
    delete curr;
}

void Adaptive_tree::print_recursion (Node *curr, std::ostream& os, size_t depth){
    //os << depth << ")";
    for(size_t i=0; i<depth; ++i){
        os << "   ";
    }
    os << "|-" << *curr;
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            if (curr->son[i_x][i_y] != NULL){
                print_recursion(curr->son[i_x][i_y], os, depth+1);
            }
        }
    }
};

Adaptive_tree::Adaptive_tree (std::vector<Particle *> &gas, size_t gas_size_max, size_t expansion_degree):
    gas_size_max(gas_size_max), expansion_degree(expansion_degree)
{
    Vec_2d pos;
    double r;
    get_bounding_box(gas, pos, r);
    root = new Node(pos, r, NULL, gas, expansion_degree);
    ctor_recursion(root);
}

Adaptive_tree::~Adaptive_tree (){
    dtor_recursion(root);
}

std::ostream& operator<<(std::ostream& os, Adaptive_tree & rha){
    os << rha.gas_size_max << "\n";
    rha.print_recursion(rha.root, os, 0);
    return os;
}

///binomial_coef_LUT

double binomial_coef_LUT::operator()(unsigned int n, unsigned int k){
    return table[n + k * size];
}

binomial_coef_LUT::binomial_coef_LUT(unsigned int size): size(size){
    table = new double [size * size + 1];
    for (unsigned int n=0; n<size; ++n){
        table[n + 0 * size] = 1;
    }
    for (unsigned int k=1; k<size; ++k){
        for (unsigned int n=k; n<size; ++n){
            table[n + k * size] = table[n + (k-1) * size] * (n-k+1)/(k);
        }
    }
}

binomial_coef_LUT::~binomial_coef_LUT(){
     delete[] table;
};

binomial_coef_LUT::binomial_coef_LUT(binomial_coef_LUT const &src): size(src.size){
    table = new double [size * size + 1];
    for (unsigned int i=0; i < size * size; ++i){
        table[i] = src.table[i];
    }
};

binomial_coef_LUT& binomial_coef_LUT::operator=(binomial_coef_LUT const &src){
    if (this == &src){
        return *this;
    }
    binomial_coef_LUT tmp(src);
    std::swap(this->size,  tmp.size);
    std::swap(this->table, tmp.table);
    return *this;
};

///FMM implementation

std::complex<double> Interaction_log_FMM::vec_to_complex(Vec_2d vec){
    return std::complex<double> (vec.x, vec.y);
}

std::vector< std::complex<double> > Interaction_log_FMM::get_powers(std::complex<double> z, unsigned int max_pow){
    std::vector< std::complex<double> > powers(max_pow+1);
    powers[0] = 1.0;
    for (size_t i=1; i<=max_pow; ++i){
        powers[i] = powers[i-1] * z;
    }
    return powers;
}

std::vector<Node *> Interaction_log_FMM::get_interaction_list(Node *curr){
    std::vector<Node *> interaction_list;
    if (curr->parent != NULL){
        for (int n_y=0; n_y<3; ++n_y){
            for (int n_x=0; n_x<3; ++n_x){
                Node *elder_neighbour = curr->parent->neighbour[n_x][n_y];
                if (elder_neighbour != NULL){
                    for (int i_y=0; i_y<2; ++i_y){
                        for (int i_x=0; i_x<2; ++i_x){
                            Node *tmp = elder_neighbour->son[i_x][i_y];
                            if(tmp != NULL){
                                Vec_2d rel_pos = tmp->pos - curr->pos;
                                double neig_r = 3*curr->r;
                                if (rel_pos.x > neig_r || rel_pos.x < -neig_r || rel_pos.y > neig_r || rel_pos.y < -neig_r){
                                    interaction_list.push_back(tmp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return interaction_list;
}

std::vector<Node *> Interaction_log_FMM::get_neighbour_list(Node *curr){
    std::vector<Node *> neighbour_list;
    for (int n_y=0; n_y<3; ++n_y){
        for (int n_x=0; n_x<3; ++n_x){
            Node *neighbour = curr->neighbour[n_x][n_y];
            if (neighbour != NULL){
                neighbour_list.push_back(neighbour);
            }
        }
    }
    return neighbour_list;
}


void Interaction_log_FMM::make_global_expansion (Node* node){
    for (auto pcl : node->gas){
        node->expansion_global[0] += pcl->charge;
        std::vector< std::complex<double> > z_power = get_powers(vec_to_complex (pcl->pos - node->pos), expansion_degree);
        for (size_t k=1; k<=expansion_degree; ++k){
            node->expansion_global[k] += pcl->charge * z_power[k];
        }
    }
    for (size_t k=1; k<=expansion_degree; ++k){
        node->expansion_global[k] *= -1.0/k;
    }
}

void Interaction_log_FMM::translate_global_expansion (Node* parent, Node* son){
    std::vector< std::complex<double> > z_power = get_powers(vec_to_complex(son->pos - parent->pos), expansion_degree);
    parent->expansion_global[0] += son->expansion_global[0];
    for (size_t i=1; i<=expansion_degree; ++i){
        parent->expansion_global[i] -= son->expansion_global[0] * z_power[i] / ((double)i);
        for (size_t k=1; k<=i; ++k){
            parent->expansion_global[i] += son->expansion_global[k] * z_power[i-k] * binomial_coef(i-1, k-1);
        }
    }
}

void Interaction_log_FMM::global_to_local_expansion (Node* source, Node* target){
    std::complex<double> z = vec_to_complex(source->pos - target->pos);
    std::vector< std::complex<double> > z_rev_power = get_powers(-1.0 / z, expansion_degree);
    std::vector< std::complex<double> > expansion_local(expansion_degree+1);

    expansion_local[0] += source->expansion_global[0] * std::log(-z);
    for (size_t k=1; k<expansion_degree; ++k){
        expansion_local[0] += source->expansion_global[k] * z_rev_power[k];
    }
    for (size_t l=1; l<expansion_degree; ++l){
        expansion_local[l] -= source->expansion_global[0] / ((double)l);
        for (size_t k=1; k<expansion_degree; ++k){
            expansion_local[l] += source->expansion_global[k] * z_rev_power[k] * binomial_coef(l+k-1, k-1);
        }
        expansion_local[l] *= z_rev_power[l] * (1.0-2*(l%2));
    }

    for (size_t i=0; i<=expansion_degree; ++i){
        target->expansion_local[i] += expansion_local[i];
    }
}

void Interaction_log_FMM::translate_local_expansion (Node* parent, Node* son){
    std::complex<double> z = vec_to_complex(son->pos - parent->pos);
    std::vector< std::complex<double> > translated_expansion = parent->expansion_local;
    for (size_t j=0; j<=expansion_degree-1; ++j){
        for (size_t k=expansion_degree-j-1; k<=expansion_degree-1; ++k){
            translated_expansion[k] += z * translated_expansion[k+1];
        }
    }
    for (size_t i=0; i<=expansion_degree; ++i){
        son->expansion_local[i] += translated_expansion[i];
    }
}


void Interaction_log_FMM::make_global_expansion_leafs (Node* curr){
    if (curr->leaf){
        make_global_expansion(curr);
        return;
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            make_global_expansion_leafs(curr->son[i_x][i_y]);
        }
    }
}

void Interaction_log_FMM::translate_global_expansion_internal (Node* curr){
    if (curr->leaf){
        return;
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            translate_global_expansion_internal(curr->son[i_x][i_y]);
            translate_global_expansion(curr, curr->son[i_x][i_y]);
        }
    }
}

void Interaction_log_FMM::global_to_local_expansion_total (Node* curr){
    std::vector<Node *> interaction_list = get_interaction_list(curr);
    for (auto source : interaction_list){
        global_to_local_expansion(source, curr);
    }
    if (curr->leaf){
        return;
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            global_to_local_expansion_total(curr->son[i_x][i_y]);
        }
    }
}

void Interaction_log_FMM::translate_local_expansion_total (Node* curr){
    if (curr->leaf){
        return;
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            translate_local_expansion(curr, curr->son[i_x][i_y]);
            translate_local_expansion_total(curr->son[i_x][i_y]);
        }
    }
}


void Interaction_log_FMM::complete_expansions(Node* root){
    make_global_expansion_leafs         (root);
    translate_global_expansion_internal (root);
    global_to_local_expansion_total     (root);
    translate_local_expansion_total     (root);
}


double Interaction_log_FMM::calc_energy_pcl_node (Particle *pcl, Node *node){
    std::complex<double> z = vec_to_complex(pcl->pos - node->pos);
    std::complex<double> potential_complex (0, 0);
    std::complex<double> z_pow = 1.0;
    for(size_t k=0; k<=expansion_degree; ++k){
        potential_complex += node->expansion_local[k] * z_pow;
        z_pow *= z;
    }
    return -factor * pcl->charge * potential_complex.real();
}

double Interaction_log_FMM::calc_energy_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    return -factor * ((pcl_1->charge * pcl_2->charge) / 2) * std::log((pcl_1->pos - pcl_2->pos).sqr());
}

Vec_2d Interaction_log_FMM::calc_force_pcl_node (Particle *pcl, Node *node){
    std::complex<double> z = vec_to_complex(pcl->pos - node->pos);
    std::complex<double> force_complex (0, 0);
    std::complex<double> z_pow = 1.0;
    for(size_t k=1; k<=expansion_degree; ++k){
        force_complex += node->expansion_local[k] * (double)k * z_pow;
        z_pow *= z;
    }
    return factor * pcl->charge * Vec_2d(force_complex.real(), -force_complex.imag());
}

Vec_2d Interaction_log_FMM::calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    return factor * ((pcl_1->charge * pcl_2->charge) / (pcl_1->pos - pcl_2->pos).sqr()) * (pcl_1->pos - pcl_2->pos);
}


double Interaction_log_FMM::calc_energy_local_expansion (Node* curr){
    if (curr->leaf){
        double energy = 0;
        for (auto pcl : curr->gas){
            energy += calc_energy_pcl_node(pcl, curr);
        }
        return energy;
    }
    double energy = 0;
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            energy+=calc_energy_local_expansion(curr->son[i_x][i_y]);
        }
    }
    return energy;
}

double Interaction_log_FMM::calc_energy_neighbours (Node* curr){
    double energy = 0;
    std::vector<Node *> neighbour_list = get_neighbour_list(curr);
    for (auto node : neighbour_list){
        if(curr->leaf || node->leaf){
            for (auto pcl_1 : curr->gas){
                for (auto pcl_2 : node->gas){
                    if (pcl_1 != pcl_2){
                        energy += calc_energy_pcl_pcl(pcl_1, pcl_2);
                    }
                }
            }
        }
    }
    if (curr->leaf){
        return energy;
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            energy+=calc_energy_neighbours(curr->son[i_x][i_y]);
        }
    }
    return energy;
}

void Interaction_log_FMM::calc_forces_local_expansion (Node* curr){
    if (curr->leaf){
        for (auto pcl : curr->gas){
            pcl->force += calc_force_pcl_node(pcl, curr);
        }
        return;
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            calc_forces_local_expansion(curr->son[i_x][i_y]);
        }
    }
}

void Interaction_log_FMM::calc_forces_neighbours (Node* curr){
    std::vector<Node *> neighbour_list = get_neighbour_list(curr);
    for (auto node : neighbour_list){
        if(curr->leaf || node->leaf){
            for (auto pcl_1 : curr->gas){
                for (auto pcl_2 : node->gas){
                    if (pcl_1 != pcl_2){
                        pcl_1->force  += calc_force_pcl_pcl(pcl_1, pcl_2);
                    }
                }
            }
        }
    }
    if (curr->leaf){
        return;
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            calc_forces_neighbours(curr->son[i_x][i_y]);
        }
    }
}

/**
void Interaction_log_FMM::print_nodes (Node* curr, int level){
    for(size_t i=0; i<level; ++i){
        std::cout << "   ";
    }
    std::cout << "|-" << *curr;
    std::cout << "\n";
    for (auto coef : curr->expansion_global){
        std::cout << coef << " ";
    }
    std::cout << "\n";
    for (auto coef : curr->expansion_local){
        std::cout << coef << " ";
    }
    std::cout << "\n";

    if (curr->leaf){
        return;
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            print_nodes(curr->son[i_x][i_y], level+1);
        }
    }
}
**/

void Interaction_log_FMM::set_expansion_degree(size_t degree){
    binomial_coef = binomial_coef_LUT(degree * 2);
    expansion_degree = degree;
}

size_t Interaction_log_FMM::get_expansion_degree(){
    return expansion_degree;
}

Interaction_log_FMM::Interaction_log_FMM (double factor, size_t gas_size_max, size_t expansion_degree):
    expansion_degree(expansion_degree), binomial_coef(expansion_degree*2), factor(factor), gas_size_max(gas_size_max) {};

double Interaction_log_FMM::calc_energy(std::vector<Particle *> &gas){
    double energy = 0;
    Adaptive_tree tree(gas, gas_size_max, expansion_degree);

    complete_expansions(tree.root);

    //print_nodes(tree.root, 0);

    energy += calc_energy_local_expansion (tree.root);
    energy += calc_energy_neighbours      (tree.root);
    energy /= 2;

    return energy;
}

void Interaction_log_FMM::calc_force (std::vector<Particle *> &gas){
    Adaptive_tree tree(gas, gas_size_max, expansion_degree);

    complete_expansions(tree.root);

    calc_forces_local_expansion (tree.root);
    calc_forces_neighbours      (tree.root);
}

double Interaction_log_FMM::calc (std::vector<Particle *> &gas){
    double energy = 0;
    Adaptive_tree tree(gas, gas_size_max, expansion_degree);

    complete_expansions(tree.root);

    energy += calc_energy_local_expansion (tree.root);
    energy += calc_energy_neighbours      (tree.root);
    energy /= 2;
    calc_forces_local_expansion (tree.root);
    calc_forces_neighbours      (tree.root);

    return energy;
}

/// 4Paraview

void print_tree_paraview_recursion(Node *curr, std::ostream& os){
    os << *curr;
    if (curr->leaf){
        return;
    }
    for (int i_y=0; i_y<2; ++i_y){
        for (int i_x=0; i_x<2; ++i_x){
            print_tree_paraview_recursion(curr->son[i_x][i_y], os);
        }
    }
}

void print_tree_paraview(Adaptive_tree &tree, int i){
    std::ofstream f_out;
    f_out.open ("data/frames/tree" + std::to_string(i) + ".txt");
    Node::cout_paraview = true;
    f_out << "x y r gas_size\n";
    print_tree_paraview_recursion(tree.root, f_out);
    f_out.close();
}
