#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>

class Stopwatch{
private:
    struct Note{
        std::string name;
        unsigned int counter;
        unsigned long long duration;

        Note(std::string& name, unsigned long long duration);
    };
    std::vector<Note> notes;
public:
    static unsigned long long get_nanoseconds();
    void propose_note(std::string& name, unsigned long long duration);
    friend std::ostream& operator<<(std::ostream& out, Stopwatch& stopwatch);
};
extern Stopwatch stopwatch;

class Stopwatch_call{
private:
    Stopwatch *stopwatch_p;
    std::string name;
    unsigned long long t_0;
public:
    Stopwatch_call(Stopwatch& stopwatch, std::string name);
    ~Stopwatch_call();
};
