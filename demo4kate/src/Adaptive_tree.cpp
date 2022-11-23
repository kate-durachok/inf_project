#include "../include/Adaptive_tree.hpp"


///Node

size_t Node::id_total = 0;

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
    size_t (*get_id)(const Node *) = Node::get_id;

    os << rha.pos << " " << rha.r << " " << rha.gas.size() << " ";
    os << "  " << get_id(&rha) << "  (" << get_id(rha.parent) << ", {" << get_id(rha.son[0][0]) << ", " << get_id(rha.son[1][0]) << ", " << get_id(rha.son[0][1]) << ", " << get_id(rha.son[1][1]) << "})";
    os << "\n";
    const size_t num_w = 2;
    os <<  "                |" << std::setw(num_w) << get_id(rha.neighbour[0][0]) << ", " << std::setw(num_w) << get_id(rha.neighbour[1][0]) << ", " << std::setw(num_w) << get_id(rha.neighbour[2][0]) << "|\n";
    os <<  "                |" << std::setw(num_w) << get_id(rha.neighbour[0][1]) << ", " << std::setw(num_w) << get_id(rha.neighbour[1][1]) << ", " << std::setw(num_w) << get_id(rha.neighbour[2][1]) << "|\n";
    os <<  "                |" << std::setw(num_w) << get_id(rha.neighbour[0][2]) << ", " << std::setw(num_w) << get_id(rha.neighbour[1][2]) << ", " << std::setw(num_w) << get_id(rha.neighbour[2][2]) << "|\n";
    return os;
}

///Tree

void adaptive_tree::subdivide(Node *curr){
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

void adaptive_tree::ctor_recursion(Node *curr){
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

void adaptive_tree::dtor_recursion(Node *curr){
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

std::vector<Node *> adaptive_tree::get_interaction_list(Node *curr){
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

void adaptive_tree::print_recursion (Node *curr, std::ostream& os, size_t depth){
    os << depth << ")";
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


adaptive_tree::adaptive_tree (std::vector<Particle *> &gas, size_t gas_size_max, size_t expansion_degree):
    gas_size_max(gas_size_max), expansion_degree(expansion_degree)
{
    Vec_2d pos;
    double r;
    get_bounding_box(gas, pos, r);
    root = new Node(pos, r, NULL, gas, expansion_degree);
    ctor_recursion(root);
}

adaptive_tree::~adaptive_tree (){
    dtor_recursion(root);
}

std::ostream& operator<<(std::ostream& os, adaptive_tree & rha){
    os << rha.gas_size_max << "\n";
    rha.print_recursion(rha.root, os, 0);
    return os;
}
