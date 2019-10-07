#ifndef MY_GEOMETRY_H
#define MY_GEOMETRY_H

// Few experiments with geometry
// Lekhraj

#include <array>

class Vec {
public:
  // Constructor creating constant vector, also default constructor
  Vec(const double x=0, const double y=0, const double z=0); // for 2D, pass first 2 values
  Vec(const std::array<double, 3>& v); // 3D
  Vec(const std::array<double, 2>& v); // 2D
  Vec(const Vec& dirV, double mag);
  
  // Copy constructor, computational cost O(n)
  Vec(const Vec &mv);
  
  // Move constructor, computational cost O(1)
  Vec(Vec &&mv) = default;
  
  // Assignment operator, computational cost O(n)
  Vec &operator=(const Vec &mv);
  
  // Move assignment operator, computational cost O(1)
  Vec &operator=(Vec &&mv) = default;
  
  // Destructor
  virtual ~Vec(void);
  
  // Type conversion to STL array
  operator std::array<double, 3>() const;

  // Access operators: rvalue \& lvalue, with range check
  double operator[](std::size_t i) const;
  double &operator[](std::size_t i);
  
  // Comparison operators
  bool operator==(const Vec &mv) const;
  bool operator!=(const Vec &mv) const; 
  
  // Overloaded arithmetic operations
  // In place vector addition: x += y;
  Vec &operator+=(const Vec &mv);
  
  // In place vector subtraction: x-= y;
  Vec &operator-=(const Vec &mv);
  
  // In place scalar multiplication: x *= a;
  Vec &operator*=(double alpha);
  
  // In place scalar division: x /= a;
  Vec &operator/=(double alpha);
  
  // Vector addition
  Vec operator+(Vec mv) const;
  
  // Vector subtraction
  Vec operator-(const Vec &mv) const;
  
  // Scalar multiplication from right and left: x = a*y; x = y*a
  Vec operator*(double alpha) const;
  friend Vec operator*(double alpha, const Vec &);
  
  // Scalar divsion: x = y/a;
  Vec operator/(double alpha) const;
  
  // Magnitude and direction (unit vector)
  double mag(void) const;
  double magsqr(void) const;
  Vec dir(void) const;
   
  // Unit Vector
  Vec unit(void) const;
  
  // Dot (inner) product : scalar
  double operator*(const Vec &) const;
  double dot(const Vec &) const;
  
  // Cross product : vector
  Vec cross(const Vec &) const;
  
  // Do two vector make acute angle (< 90 deg)?
  bool isAcute(const Vec&) const;
  
  // Do two vectors make an right angle (90 deg)?
  bool isRight(const Vec&) const;
  
  // Does given vector fall on anticlockwise direction from our vector?
  bool isAntiClock(const Vec&) const;
  
  // Project on another vector
  double project(const Vec&) const;  
  
  // Output function
  friend std::ostream &operator<<(std::ostream &, const Vec &mv);

public: // Leaving data member spublic so that x,y,z components can be accessed easily
  double x, y, z; // Data: Vector components
};


#endif
