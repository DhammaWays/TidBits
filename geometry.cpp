// Experiment with geometry
// Lekhraj

#include <limits>
#include <cmath>
#include <iostream>

#include "geometry.h"

// Implementation of Vector class member functions
Vec::Vec(const double _x, const double _y, const double _z) : x(_x), y(_y), z(_z) {}

Vec::Vec(const std::array<double, 3>& v) : x(v[0]), y(v[1]), z(v[2]) {}
Vec::Vec(const std::array<double, 2>& v) : x(v[0]), y(v[1]), z(0) {}
Vec::Vec(const Vec& dirV, double mag) : x(mag*dirV.x), y(mag*dirV.y), z(mag*dirV.z) {}

// Copy constructor
Vec::Vec(const Vec &mv) : x(mv.x), y(mv.y), z(mv.z) {}

// Assignment operator
Vec &Vec::operator=(const Vec &mv) {
  if (this == &mv) // Self assignment
    return (*this);
	
  x = mv.x;
  y = mv.y;
  z = mv.z;
  
  return (*this);
}

// Destructor
Vec::~Vec(void) {}


Vec::operator std::array<double, 3>() const {
  return (std::array<double, 3>({{x, y, z}}));
}

double Vec::operator[](std::size_t i) const {
  if (i >= 3)
    throw(std::logic_error("[] out of range"));
	
  return ((std::array<double, 3>)(*this))[i];
}

double &Vec::operator[](std::size_t i) {
  if (i >= 3)
    throw(std::logic_error("[] out of range"));
	
  return ((std::array<double, 3>)(*this))[i];
}

// Compare: v1 == v2
bool Vec::operator==(const Vec &mv) const {
  if (x != mv.x || y != mv.y || z != mv.z )
        return (false);

  // Two vector are same: v1 == v2		
  return (true);
}

bool Vec::operator!=(const Vec &mv) const { return !(*this == mv); }

// In-place Add: V += v1
Vec &Vec::operator+=(const Vec &mv) {
  x += mv.x;
  y += mv.y;
  z += mv.z;
  
  return (*this);
}

// In-place Subtract: V -= v1
Vec &Vec::operator-=(const Vec &mv) {
  x -= mv.x;
  y -= mv.y;
  z -= mv.z;
  
  return (*this);
}

// In-place: V *= a
Vec &Vec::operator*=(double scale) {
  x *= scale;
  y *= scale;
  z *= scale;
  
  return (*this);
}

// In-place V /= a
Vec &Vec::operator/=(double a) {
  if (a == 0)
    throw(std::logic_error("Divide by zero!"));
	
  x /= a;
  y /= a;
  z /= a;
  
  return (*this);
}

// v1 + V
Vec Vec::operator+(Vec mv) const {
  mv += *this;
  
  return (mv);
}

// V-v1
Vec Vec::operator-(const Vec &mv) const {
  Vec tmp(*this);
  tmp -= mv;
  
  return (tmp);
}

// 
Vec Vec::operator*(double scale) const {
  Vec tmp(*this);
  tmp *= scale;
  
  return (tmp);
}

// V*a
Vec operator*(double scale, const Vec &mv) {
  Vec tmp(mv);
  tmp *= scale;
  
  return (tmp);
}

// V/a
Vec Vec::operator/(double a) const {
  Vec tmp(*this);
  tmp /= a;
  
  return (tmp);
}

double Vec::mag(void) const {
  double s = x*x + y*y + z*z;
  
  return (std::sqrt(s));
}

double Vec::magsqr(void) const {
  double s = x*x + y*y + z*z;
  
  return (s);
}

Vec Vec::dir(void) const {
  return this->unit();
}

Vec Vec::unit(void) const {
  Vec tmp(*this);
  
  return (tmp/mag());
}

// Dot (Scalar Product): V.v1
// (x1, y1, z1).(x2, y2, z2) => x1*x2 + y1*y2 + z1*z2
double Vec::dot(const Vec &mv) const {
  return (x * mv.x + y * mv.y + z * mv.z);
}

// Dot (Scalar Product): V.v1
double Vec::operator*(const Vec &mv) const {  
  return dot(mv);
}

// Cross Product (Vector Product): V x v1
// 	| i  j  k  |
//	| x1 y1 z1 | => (y1*z2 - z1*y2)i + (z1*x2 - x1*z2)j + (x1*y2 - y1*x2)k
//	| x2 y2 z2 |
Vec Vec::cross(const Vec &mv) const {
  return Vec(y * mv.z - z * mv.y, z * mv.x - x * mv.z, x * mv.y - y * mv.x);
}

// Does this vector make an acute angle (<=90 deg) with given vector?
bool Vec::isAcute(const Vec& v) const {
	return this->dot(v) > 0;
}

// Does this vector make a right angle (90 deg) with given vector?
bool Vec::isRight(const Vec& v) const {
	return this->dot(v) == 0;
}

// Does given vector go in anticlock wise direction from our vector?
// Treat collinear vector as anticlockwise.
bool Vec::isAntiClock(const Vec& v) const {
	return (x * v.y - y * v.x) >= 0; // this->cross(v).z >= 0
}

// Project his vector on to given vector
double Vec::project(const Vec& v) const {
	return this->dot(v)/v.mag();
}

// Handy print
std::ostream &operator<<(std::ostream &o, const Vec &mv) {
  return (o << "[" << mv.x << ", " << mv.y << ", " << mv.z << "]");
}

