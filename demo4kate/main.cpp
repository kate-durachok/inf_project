#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>
#include <set>
#include <complex>
#include <climits>
#include <string>
#include <fstream>
#include <iomanip>

#include "include/Gas.hpp"
#include "include/Interaction_log_pairwise.hpp"
#include "include/Interaction_repultion.hpp"
#include "include/Interaction_trap.hpp"
#include "include/Adaptive_tree.hpp"

double log_interaction_coef = 1;

unsigned long long f_rand (){
    static unsigned long long seed = 1758750;
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

int main()
{
    std::cout.precision(3);
    std::cout << std::fixed;

    ///создаём газ
    int n = 12;
    std::vector<Particle *> gas;

    //gas.push_back(new Particle(Vec_2d(0.0, 0.0), Vec_2d(0, 0), 1.0, 1.0));
    //gas.push_back(new Particle(Vec_2d(8.0, 8.0), Vec_2d(0, 0), 1.0, 1.0));

    for (int i=0; i<n; ++i){
        Vec_2d pos(f_rand_double(0, 8), f_rand_double(0, 8));
        Vec_2d vel(0.0, 0.0);
        double mass = 1.0;
        double charge = 0.05;
        gas.push_back(new Particle(pos, vel, mass, charge));
    }

    ///это нинада
//    adaptive_tree tree(gas, 3, 0);
//    std::cout << tree;

    ///Объявление взаимодействий, которые будут действовать на частицы
    std::vector<Interaction_base *> inter_arr;
    Interaction_log_pairwise  inter_1(1.0);
    Interaction_repulsion     inter_2(1.0, 0.1, 1, 2);
    Interaction_trap          inter_3(Vec_2d(0, 0), 25.0, 1.0, 10);
    inter_arr.push_back(&inter_1);
    inter_arr.push_back(&inter_2);
    inter_arr.push_back(&inter_3);



    //print_frame(gas, 0); <- пишет все координаты частиц в файл для сторонней обработке
    ///основной цикл
    double dt = 1.0;
    int subframe_amm = 10;
    for (int i=1; i<=100; ++i){
        for (int k=0; k<subframe_amm; ++k){
            verlet(gas, inter_arr, dt/subframe_amm);  // тут происходит движение частиц
        }
        double potential = calc_potential_energy(gas, inter_arr);
        std::cout << "iter_n:" << i << "\n";
        std::cout << calc_kinetic_energy(gas) << "   " << potential << "\n";
        std::cout << calc_kinetic_energy(gas) + potential << "\n";
        std::cout << "------------" << "\n";
        //print_frame(gas, i);

        /// в этой точке можно брать координаты газа и рисовать.
    }

    for (auto pcl : gas){
        delete pcl;
    }

    return 0;
}
