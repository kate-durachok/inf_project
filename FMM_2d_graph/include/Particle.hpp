#pragma once

#include "Vec_2d.hpp"

class Particle{
private:
//    static size_t id_total = 1;
//    size_t id;
public:
    Vec_2d pos;
    Vec_2d vel;
    Vec_2d acc, acc_prev;
    Vec_2d force;

    double mass;
    double charge;

    Particle(Vec_2d pos, Vec_2d vel, double mass, double charge): pos(pos), vel(vel), mass(mass), charge(charge) {};
//    Particle(Vec_2d pos, Vec_2d vel, double mas, double charge) pos(pos), vel(vel), mas(mas), charge(charge), id(id_total) {++id_total};
//    Particle(Particle const &src) Particle(src.pos, src.vel, src.mas, src.charge) {};
//    Particle& operator=(Particle const &src)
//
//    size_t get_id(){
//        return id;
//    }

    friend std::ostream& operator <<(std::ostream& os, const Particle & rha){
        return os << "pos:" << rha.pos << "   vel:" << rha.vel << "   acc:" << rha.acc << "   mass:" << rha.mass << "   charge:" << rha.charge << "\n";
    };
    std::ostream& print_minimal(std::ostream& os){
        return os << pos.x << " " << pos.y << " " << vel.x << " " << vel.y << " " << mass << " " << charge << "\n";
    }
};
