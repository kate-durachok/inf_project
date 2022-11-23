#pragma once

#include <cmath>
#include <iostream>

class Vec_2d{
private:

public:
    double x, y;
    double &operator[](size_t);

    Vec_2d(double, double);
    Vec_2d();

    Vec_2d(Vec_2d const &) = default;
    ~Vec_2d() = default;
    Vec_2d& operator=(Vec_2d const &) = default;

    Vec_2d operator+(Vec_2d const &) const;
    Vec_2d operator-() const;
    Vec_2d operator-(Vec_2d const &) const;
    Vec_2d& operator+=(Vec_2d const &);
    Vec_2d& operator-=(Vec_2d const &);

    Vec_2d operator*(double const &) const;
    friend Vec_2d operator*(double const &, Vec_2d const &);
    Vec_2d operator/(double const &) const;
    Vec_2d& operator*=(double const &);
    Vec_2d& operator/=(double const &);

    double operator*(Vec_2d const &) const;
    double sqr() const;
    double len() const;

    friend std::ostream& operator<<(std::ostream& os, const Vec_2d& rha);
};
