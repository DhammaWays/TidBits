// Experiment with geometry
// Lekhraj

#include <limits>
#include <cmath>
#include <iostream>

#include "geometry.h"
#include "numbers.h"

// Implementation of Vector class member functions
Vec::Vec(const double _x, const double _y, const double _z) : x(_x), y(_y), z(_z) {}

Vec::Vec(const std::array<double, 3>& v) : x(v[0]), y(v[1]), z(v[2]) {}
Vec::Vec(const std::array<double, 2>& v) : x(v[0]), y(v[1]), z(0) {}
Vec::Vec(const Vec& dirV, double mag) : x(mag*dirV.x), y(mag*dirV.y), z(mag*dirV.z) {}
Vec::Vec(const Point& p0, const Point& p1) : x(p1.x - p0.x), y(p1.y - p0.y), z(p1.z - p0.z) {}
  

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

// Shapes

// Do two bound boxes intersect with each other?
bool doBoundsOverlap(const BBox& b1, const BBox& b2) {
	// If either of upper corner is less then other's lower corner we do not overlap
	return !( b1.ur.x < b2.ll.x || b1.ur.y < b2.ll.y || b1.ur.z < b2.ll.z ||
			  b2.ur.x < b1.ll.x || b2.ur.y < b1.ll.y || b2.ur.z < b1.ll.z );
}

// Handy print
std::ostream &operator<<(std::ostream &o, const Shape &s) {
  return (s.print(o));
}

// Line

Line::Line(const Point& p0, const Point& p1) : first(p0), second(p1) {}

BBox Line::boundBox(void) const {
	BBox bound;
	bound.ll.x = MIN(first.x, second.x);
	bound.ll.y = MIN(first.y, second.y);
	bound.ll.z = MIN(first.z, second.z);
	bound.ur.x = MAX(first.x, second.x);
	bound.ur.y = MAX(first.y, second.y);
	bound.ur.z = MAX(first.z, second.z);
	
	return bound;
}

// Line to Line intersection
bool Line::doesIntersect(const Shape& s) const {
	// Line to Line intersection Algo:
	// 1) Do NOT intersect if their bound boxes do not overlap
	// 2) Do intersect if they are collinear (lie along the same line) 
	//    Case when their ends do not overlap will get rejected with bound box check
	// 3) Do NOT intersect if both ends of lines do not lie on same side w.r.t each other
	// 4) Rest should intersect?
	
	const Line* const pL2 = dynamic_cast<const Line* const>(&s);
	assert( pL2 );
	if( !pL2 ) {
		return false;
		}

	BBox bound1 = this->boundBox();
	BBox bound2= pL2->boundBox();
	
	if( !doBoundsOverlap(bound1, bound2) )
		return false; // No intersection:  L1 and L2 are not within each other bounds
		
	Vec v1f2f(first, pL2->first), v1f2s(first, pL2->second);
	if( v1f2f.cross(v1f2s).magsqr() == 0 ) // Both the lines are collinear
		return true; // Intersection: They must be overlapping since boundbox do overlap
		
	Vec v1f1s(first, second); 
	if( v1f1s.isAntiClock(v1f2f) == v1f1s.isAntiClock(v1f2s) ) // L2 end points are on same side of L1
		return false; // No intersection: L2 is completely on one side of L1 
		
	Vec v2f2s(pL2->first, pL2->second), v2f1f(pL2->first, first), v2f1s(pL2->first, second);
	if( v2f2s.isAntiClock(v2f1f) == v2f2s.isAntiClock(v2f1s) ) // L1 end points are on same side of L1
		return false; // No intersection: L1 is completely on one side of L2 

	// If we have reached so far, lines must be intersecting.
	// Their bound box overlap, they are not collinear and their end points lie on opposite of each other.				
	return true;
}

bool Line::isInside(const Point& p) const {
	return false;
}

double Line::area(void) const {
	return 0;
}
	
bool Line::isConvex(void) const {
	return true;
} 

std::ostream& Line::print(std::ostream &o) const {
  return (o << "Line : P0[" << first.x << ", " << first.y << ", " << first.z << "]" <<
  		  " to P1[" << second.x << ", " << second.y << ", " << second.z << "] ");
}
