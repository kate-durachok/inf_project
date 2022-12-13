#include "../include/Vec_2d.hpp"

double& Vec_2d::operator[](size_t ind){
    if (ind == 0) {return x;}
    if (ind == 1) {return y;}
    return x;
}

Vec_2d::Vec_2d(double x, double y):x(x), y(y){}

Vec_2d::Vec_2d():Vec_2d(0, 0){}

Vec_2d Vec_2d::operator+(Vec_2d const &rha) const{
    return Vec_2d(x + rha.x, y + rha.y);
}

Vec_2d Vec_2d::operator-() const{
    return Vec_2d(-x, -y);
}

Vec_2d Vec_2d::operator-(Vec_2d const &rha) const{
    return *this + (-rha);
}

Vec_2d& Vec_2d::operator+=(Vec_2d const &rha){
    *this = *this + rha;
    return *this;
}

Vec_2d& Vec_2d::operator-=(Vec_2d const &rha){
    *this = *this - rha;
    return *this;
}


Vec_2d Vec_2d::operator*(double const &k) const{
    return Vec_2d(k * this->x, k * this->y);
}

Vec_2d operator*(double const &k, Vec_2d const &rha){
    return rha * k;
}

Vec_2d Vec_2d::operator/(double const &k) const{
    return *this * (1/k);
}

Vec_2d& Vec_2d::operator*=(double const &k){
    *this = *this * k;
    return *this;
}

Vec_2d& Vec_2d::operator/=(double const &k){
    *this = *this / k;
    return *this;
}


double Vec_2d::operator*(Vec_2d const &rha) const{
    return x*rha.x + y*rha.y;
}

double Vec_2d::sqr() const{
    return (*this) * (*this);
}

double Vec_2d::len() const{
    return std::sqrt(this->sqr());
}

std::ostream& operator<<(std::ostream &os, const Vec_2d &rha) {
    os << "(" << rha.x << ", " << rha.y << ")";
    return os;
}
