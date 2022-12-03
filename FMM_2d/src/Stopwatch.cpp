#include "../include/Stopwatch.hpp"

Stopwatch::Note::Note(std::string& name, unsigned long long duration):name(name), counter(1), duration(duration) { }
unsigned long long Stopwatch::get_nanoseconds(){
    static auto t_0 = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_0).count();
}

void Stopwatch::propose_note(std::string& name, unsigned long long duration){
    for(auto &note : notes){
        if (name == note.name){
            ++note.counter;
            note.duration += duration;
            return;
        }
    }
    notes.push_back(Note(name, duration));
}

std::ostream& operator<<(std::ostream& out, Stopwatch& stopwatch){
    if (!stopwatch.notes.size()){
        return out;
    }
    out << "Stopwatch measurements:\n";
    size_t length_max = 0;
    for (auto note : stopwatch.notes){
        if (length_max < note.name.length()){
            length_max = note.name.length();
        }
    }
    for (auto note : stopwatch.notes){
        out << "|   " << note.name + ":";
        for (size_t i=0; i<length_max-note.name.size()+1; ++i){
            out << " ";
        }
        std::cout.precision(1);
        std::cout << std::fixed;
        out << std::setw(5);
        if (note.duration > 1E9) {
            out << note.duration * 1E-9 << " sec";
        } else if (note.duration > 1E6) {
            out << note.duration * 1E-6 << " ms ";
        } else {
            out << note.duration * 1E-3 << " mcs";
        }

        out << "   n=" << note.counter << "\n";
    }
    stopwatch.notes.clear();
    return out;
}
Stopwatch stopwatch;


Stopwatch_call::Stopwatch_call(Stopwatch& stopwatch, std::string name):stopwatch_p(&stopwatch), name(name){
    t_0 = Stopwatch::get_nanoseconds();
}

Stopwatch_call::~Stopwatch_call(){
    unsigned long long t_1 = Stopwatch::get_nanoseconds();
    stopwatch_p->propose_note(name, t_1 - t_0);
}

