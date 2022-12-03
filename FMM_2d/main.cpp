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
#include "include/Interaction_log_FMM.hpp"
#include "include/Interaction_short.hpp"
#include "include/Interaction_external.hpp"
#include "include/Stopwatch.hpp"

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

int main()
{
    std::vector<Particle *> gas;
    //read_gas(gas, "data/saved_frames/circle_1k.txt");
    {
        int size_x = 50;
        int size_y = 40;
        double dist = 1.12;
        for (int i_y=0; i_y<size_y; ++i_y){
            for (int i_x=0; i_x<size_x; ++i_x){
                double x = (i_x - size_x/2) * 0.866;
                double y = (i_y + i_x%2*0.5 - size_y/2);
                Vec_2d pos(dist * x, dist * y);
                Vec_2d vel(f_rand_double(-0.5, 0.5), f_rand_double(-0.5, 0.5));
                double mass = 1.0;
                double charge = 0.01 * (2*(i_x%2) - 1);
                gas.push_back(new Particle(pos, vel, mass, charge));
            }
        }
    }
    {
        Vec_2d momentum(0, 0);
        double mass = 0;
        for (auto pcl : gas){
            momentum += pcl->vel * pcl->mass;
            mass += pcl->mass;
        }
        for (auto pcl : gas){
            pcl->vel -= momentum/mass;
        }
    }

    Interaction_log_pairwise  inter_pair (-1.0);
    int expansion_degree = 5;
    int gas_size_max = 10;
    Interaction_log_FMM       inter_1 (1.0, gas_size_max, expansion_degree);
    Interaction_repulsion     inter_2 (1.0, 0.1, 1, 2);
    Interaction_6_12_smoothed inter_3 (0.01, 1.0, 2.0, 1.5, 2);
    Interaction_trap          inter_4 (Vec_2d(0, 0), 60.0, 1.0, 3);
    Interaction_uni_field     inter_5 (Vec_2d(0.1, 0));

    std::vector<Interaction_base *> inter_arr;
    //inter_arr.push_back(&inter_pair);
    inter_arr.push_back(&inter_1);
    //inter_arr.push_back(&inter_2);
    inter_arr.push_back(&inter_3);
    inter_arr.push_back(&inter_4);
    inter_arr.push_back(&inter_5);

    print_frame(gas, 0);
    double dt = 5.00;
    int subframe_amm = 70;
    for (int i=1; i<=1000; ++i){
        unsigned long long t_0 = Stopwatch::get_nanoseconds();
        for (int k=0; k<subframe_amm; ++k){
            verlet(gas, inter_arr, dt/subframe_amm);
            for (auto pcl : gas){ pcl->vel -= pcl->vel * (0.002 * dt/subframe_amm); }
        }
        unsigned long long t_1 = Stopwatch::get_nanoseconds();
        double ener_p = calc_potential_energy(gas, inter_arr);
        double ener_k = calc_kinetic_energy(gas);

        std::cout << "------------------------\n";
        std::cout.precision(3);
        std::cout << std::fixed;
        std::cout << "iter_n:" << i << "\n";
        printf("%6.3E   %6.3E     %6.7E\n", ener_k, ener_p, ener_k + ener_p);
        std::cout << (t_1-t_0)*1.0E-9 << "\n";
        std::cout << (t_1-t_0)*1.0E-3/(subframe_amm * gas.size()) << "\n";
        std::cout << stopwatch;
        std::cout << "------------------------\n\n";

        print_frame(gas, i);

        Adaptive_tree tree(gas, gas_size_max, 0);
        print_tree_paraview(tree, i);
    }

    for (auto pcl : gas){
        delete pcl;
    }

    return 0;
}
