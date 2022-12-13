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
#include <SFML/Graphics.hpp>

#include "include/Gas.hpp"
#include "include/Interaction_log_pairwise.hpp"
#include "include/Interaction_log_FMM.hpp"
#include "include/Interaction_short.hpp"
#include "include/Interaction_external.hpp"
#include "include/Stopwatch.hpp"

unsigned long long f_rand ()
{
    static unsigned long long seed = 100;
    static unsigned long long prev_random = seed;
    prev_random ^= prev_random << 21;
    prev_random ^= prev_random >> 35;
    prev_random ^= prev_random << 4;
    return prev_random;
}

double f_rand_double(double _min, double _max)
{
    return (1.0 * f_rand()/ULLONG_MAX) * (_max - _min) + _min;
}


void updateParticle(std::vector<Particle *> &gas, std::vector<Interaction_base *> &inter_arr, double dt, int subframe_amm)
{
        for (int k=0; k<subframe_amm; ++k){
            verlet(gas, inter_arr, dt/subframe_amm);
            for (auto pcl : gas){ pcl->vel -= pcl->vel * (0.002 * dt/subframe_amm); }
        }
}


int velTocol(double vel)
{
    return std::min((int)std::round((255 - 100) / (0.8 - 0) * vel + 100), 255);
}

bool check_pcl_vs_gas(std::vector<Particle *> &gas, Vec_2d pos, double PARTICLE_RADIUS){
    for (auto pcl : gas){
        if ( (pcl->pos - pos).sqr() < PARTICLE_RADIUS*PARTICLE_RADIUS){
            return false;
        }
    }
    return true;
}

void make_some_gas(std::vector<Particle *> &gas, std::vector<sf::CircleShape *> &particles, Vec_2d pos_0, int size_x, int size_y, double PARTICLE_RADIUS){
    double dist = 2*PARTICLE_RADIUS*1.5;
    for (int i_y=0; i_y<size_y; ++i_y){
        for (int i_x=0; i_x<size_x; ++i_x){
            double x = (i_x - (size_x-1)/2.0);
            double y = (i_y - (size_y-1)/2.0);
            Vec_2d pos = Vec_2d(x, y) * dist + pos_0;
            if (check_pcl_vs_gas(gas, pos, 2*PARTICLE_RADIUS)){
                Vec_2d vel(f_rand_double(-0.5, 0.5), f_rand_double(-0.5, 0.5));
                double mass = 1.0;
                double charge = 0.02 * (1.0*(f_rand()%3) - 1);
                gas.push_back(new Particle(pos, vel, mass, charge));
                particles.push_back(new sf::CircleShape(PARTICLE_RADIUS));
                particles[i_x + i_y * size_y]->setPosition(pos.x, pos.y);
                sf::Color particleColor(0, 0, 0);
                if (charge < 0)
                    particleColor.b = velTocol(vel.len());
                else if (charge > 0)
                    particleColor.r = velTocol(vel.len());
                else
                    particleColor.g = velTocol(vel.len());
                particles[i_x + i_y * size_y]->setFillColor(particleColor);
            }
        }
    }
}

int main()
{

//    sf::Font font;
//    font.loadFromFile("ITCFranklinGothic.ttf");
//    sf::Text text;
//    text.setColor(sf::Color::White);
//    text.setCharacterSize(24);

    std::vector<Particle *> gas;
    std::vector<sf::CircleShape *> particles;
    const unsigned int windowW = 1200,
                        windowH = 1200;
    const unsigned int PARTICLE_RADIUS = 5;
    sf::RenderWindow window(sf::VideoMode(windowW, windowH), "I hate this!", sf::Style::Close);
    //read_gas(gas, "data/saved_frames/circle_1k.txt");
    {
        int size_x = 10;
        int size_y = 10;
        Vec_2d pos(0, 0);
        make_some_gas(gas, particles, pos, size_x, size_y, PARTICLE_RADIUS);
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

    const float TRAP_RADIUS = 300.0;
    sf::CircleShape *trapCircle = new sf::CircleShape(TRAP_RADIUS);
    Interaction_log_pairwise  inter_pair (-1.0);
    int expansion_degree = 5;
    int gas_size_max = 10;
    Interaction_log_FMM       inter_1 (1.0, gas_size_max, expansion_degree);
    Interaction_repulsion     inter_2 (1.0, 0.1, 1, 2);
    Interaction_6_12_smoothed inter_3 (0.02, 2*PARTICLE_RADIUS, 2*PARTICLE_RADIUS*2.0, 2*PARTICLE_RADIUS*1.5, 2);
    Interaction_trap          inter_4 (Vec_2d(0, 0), TRAP_RADIUS, 1.0, 3);
    Interaction_uni_field     inter_5 (Vec_2d(0.05 , 0));

    std::vector<Interaction_base *> inter_arr;
    inter_arr.push_back(&inter_pair);
    inter_arr.push_back(&inter_1);
    //inter_arr.push_back(&inter_2);
    inter_arr.push_back(&inter_3);
    inter_arr.push_back(&inter_4);
    inter_arr.push_back(&inter_5);

    //print_frame(gas, 0);
    sf::Color backgroundColor = sf::Color::White;
    trapCircle->setPosition(windowW / 2 - trapCircle->getRadius(), windowH / 2 - trapCircle->getRadius());
    trapCircle->setFillColor(backgroundColor);
    trapCircle->setOutlineColor(sf::Color::Black);
    trapCircle->setOutlineThickness(PARTICLE_RADIUS);
    trapCircle->setPointCount(300);
    bool updFlag = true;
    double dt = 1;
    int subframe_amm = 5;
    float scaleFactor = 1.0;

    for (int i=1; i<=30000; ++i)
    {
        unsigned long long t_0 = Stopwatch::get_nanoseconds();
        if(updFlag)
        {
            updateParticle(gas, inter_arr, dt, subframe_amm);
            for(int j = 0; j < gas.size(); ++j)
            {
                particles[j]->setPosition(scaleFactor * (gas[j]->pos.x - PARTICLE_RADIUS) + windowW / 2,
                                          scaleFactor * (gas[j]->pos.y - PARTICLE_RADIUS) + windowH / 2);
                sf::Color tmpColor = particles[j]->getFillColor();

                if (gas[j]->charge > 0){
                    tmpColor.r = velTocol(gas[j]->vel.len());
                    tmpColor.b = 0;
                    tmpColor.g = 0;
                }
                else if (gas[j]->charge < 0)
                {
                    tmpColor.r = 0;
                    tmpColor.b = velTocol(gas[j]->vel.len());
                    tmpColor.g = 0;
                }
                else
                {
                    tmpColor.r = 0;
                    tmpColor.b = 0;
                    tmpColor.g = velTocol(gas[j]->vel.len());;
                }
                particles[j]->setFillColor(tmpColor);
            }
        }
        else
        {
            --i;
        }

//        double ener_p = calc_potential_energy(gas, inter_arr);
//        double ener_k = calc_kinetic_energy(gas);



//        std::cout << "------------------------\n";
//        std::cout.precision(3);
//        std::cout << std::fixed;
//        std::cout << "iter_n:" << i << "\n";
//        printf("%6.3E   %6.3E     %6.7E\n", ener_k, ener_p, ener_k + ener_p);
//        std::cout << (t_1-t_0)*1.0E-9 << "\n";
//        std::cout << (t_1-t_0)*1.0E-3/(subframe_amm * gas.size()) << "\n";
//        std::cout << stopwatch;
//        std::cout << "------------------------\n\n";


        sf::Event event;
        while (window.pollEvent(event))
        {
            // Close window: exit
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == sf::Event::KeyPressed)
            {
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space))
                {
                    updFlag = !updFlag;
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
                {
                    scaleFactor *= 1.25;
                    trapCircle->setRadius(scaleFactor * TRAP_RADIUS);
                    trapCircle->setPosition(windowW / 2 - trapCircle->getRadius(),
                                            windowH / 2 - trapCircle->getRadius());
                    for(int j = 0; j < particles.size(); ++j)
                    {
                        particles[j]->setRadius(scaleFactor * PARTICLE_RADIUS);
                    }
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
                {
                    scaleFactor /= 1.25;
                    trapCircle->setRadius(scaleFactor * TRAP_RADIUS);
                    trapCircle->setPosition(windowW / 2 - trapCircle->getRadius(),
                                            windowH / 2 - trapCircle->getRadius());
                    for(int j = 0; j < particles.size(); ++j)
                    {
                        particles[j]->setRadius(scaleFactor * PARTICLE_RADIUS);
                    }
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::P))
                {
                    double x = sf::Mouse::getPosition(window).x - PARTICLE_RADIUS * scaleFactor;
                    double y = sf::Mouse::getPosition(window).y - PARTICLE_RADIUS * scaleFactor;
                    Vec_2d pos((x - windowW / 2) / scaleFactor,
                               (y - windowH / 2) / scaleFactor);
                    int size_x = 4;
                    int size_y = 4;
                    make_some_gas(gas, particles, pos, size_x, size_y, PARTICLE_RADIUS);
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::N))
                {
                    double x = sf::Mouse::getPosition(window).x - PARTICLE_RADIUS * scaleFactor;
                    double y = sf::Mouse::getPosition(window).y - PARTICLE_RADIUS * scaleFactor;
                    Vec_2d pos((x - windowW / 2) / scaleFactor,
                               (y - windowH / 2) / scaleFactor);
                    int size_x = 1;
                    int size_y = 1;
                    make_some_gas(gas, particles, pos, size_x, size_y, PARTICLE_RADIUS);
                }
            }
        }

        window.clear(sf::Color::Black);
        window.draw(*trapCircle);
        for(auto particle: particles)
        {
            window.draw(*particle);
        }

//        window.draw(text);
        unsigned long long t_1 = Stopwatch::get_nanoseconds();
        window.setTitle(std::to_string(1e9/(t_1 - t_0)));
        window.display();
    }

    for (auto pcl : gas){
        delete pcl;
    }

    return 0;
}
