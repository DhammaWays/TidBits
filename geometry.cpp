// Experiment with geometry
// Lekhraj

#include <limits>
#include <cmath>
#include <iostream>

#include "geometry.h"
#include "numbers.h"


// Handy print
std::ostream& printPoint(std::ostream &o, const Point& p) {
  o << "[" << p.x << ", " << p.y << ", " << p.z << "]";
  return (o);
}

std::ostream& operator<<(std::ostream& o, const Point& p)
{
  o << "[" << p.x << ", " << p.y << ", " << p.z << "]";
  return (o);
}

// Vector
 
// Implementation of Vector class member functions
Vec::Vec(const double _x, const double _y, const double _z) : x(_x), y(_y), z(_z) {}

Vec::Vec(const std::array<double, 3>& v) : x(v[0]), y(v[1]), z(v[2]) {}
Vec::Vec(const std::array<double, 2>& v) : x(v[0]), y(v[1]), z(0) {}
Vec::Vec(const Vec& dirV, double mag) : x(mag*dirV.x), y(mag*dirV.y), z(mag*dirV.z) {}
Vec::Vec(const Point& p0, const Point& p1) : x(p1.x - p0.x), y(p1.y - p0.y), z(p1.z - p0.z) {}
Vec::Vec(const Point& p) : x(p.x), y(p.y), z(p.z) {}  

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

// Does this vector make an acute angle (< 90 deg) with given vector?
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
  return (o << Point({mv.x, mv.y, mv.z}));
}

// Shapes

// Do two bound boxes intersect with each other?
bool doBoundsOverlap(const BBox& b1, const BBox& b2) {
	// If either of upper corner is less then other's lower corner we do not overlap
	return !( b1.ur.x < b2.ll.x || b1.ur.y < b2.ll.y || b1.ur.z < b2.ll.z ||
			  b2.ur.x < b1.ll.x || b2.ur.y < b1.ll.y || b2.ur.z < b1.ll.z );
}

// Do two lines share an end point?
bool haveCommonEndPoint(const Line* pL1, const Line* pL2) {
	return ((pL1->first) == (pL2->first) || (pL1->first) == (pL2->second) || 
		    (pL1->second) == (pL2->first) || (pL1->second) == (pL2->second));
}

// Is the point inside bound box?
bool isPointInsideBoundBox(const BBox& b, const Point& p) {
	return !( p.x < b.ll.x || p.y < b.ll.y || p.z < b.ll.z ||
	          p.x > b.ur.x || p.y > b.ur.y || p.z > b.ur.z );
}

// Does an horizontal ray in +X axis direction intersect given line
// This test is much faster then general line to line intersection
bool doesRayIntersectLine(const Line& ray, const Line line) {
	BBox bound = line.boundBox();
	
	// For +X axis horizontal ray to interscet a line, the line must
	// lie to left of line and ray vertical position must be 
	// within vertical bounds of the line.
	if( bound.ll.y <= ray.first.y && ray.first.y <= bound.ur.y ) {
		if( ray.first.x <= bound.ll.x ) // intersection: ray starting point is completely left of line bounds
			return true;
		else if ( ray.first.x > bound.ur.x ) // no interscetion: ray starting point is completely right of line bounds
			return false;
		else { // Need to check if ray starting point lies on left of line
			Line increasingYline(line);			
			if( line.first.y > line.second.y ) { // Need to swap end points to make our line in increasing Y order
				increasingYline.first = line.second;
				increasingYline.second = line.first;
			}
			
			Vec lineV(increasingYline.first, increasingYline.second);
			Vec rayV(increasingYline.first, ray.first);
			return lineV.isAntiClock(rayV);			
		}
	}
	else
		return false;
	
}

// Handy print
std::ostream &operator<<(std::ostream &o, const Shape &s) {
  return (s.print(o));
}

// Line

Line::Line(const Point& p0, const Point& p1) : first(p0), second(p1) {}

BBox Line::boundBox(void) const {
	BBox bound;
	
	/*
	bound.ll.x = MIN(first.x, second.x);
	bound.ll.y = MIN(first.y, second.y);
	bound.ll.z = MIN(first.z, second.z);
	bound.ur.x = MAX(first.x, second.x);
	bound.ur.y = MAX(first.y, second.y);
	bound.ur.z = MAX(first.z, second.z);
	*/
		
	// MINMAX along x-axis
	if( first.x <= second. x ) {
		bound.ll.x = first.x;
		bound.ur.x = second.x;
	}
	else {
		bound.ll.x = second.x;
		bound.ur.x = first.x;
	}
	
	// MINMAX along y-axis
	if( first.y <= second.y ) {
		bound.ll.y = first.y;
		bound.ur.y = second.y;
	}
	else {
		bound.ll.y = second.y;
		bound.ur.y = first.y;
	}
	
	// MINMAX along z-axis
	if( first.z <= second.z ) {
		bound.ll.z = first.z;
		bound.ur.z = second.z;
	}
	else {
		bound.ll.z = second.z;
		bound.ur.z = first.z;
	}		
	
	return bound;
}

// Line to Line intersection
bool Line::doesIntersect(const Shape& s) const {
	// Line to Line intersection Algo:
	// 1) Do NOT intersect if their bound boxes do not overlap
	// 2) Do intersect if they are collinear (lie along the same line) or their one of end point is on other line
	//    Case when their ends do not overlap will get rejected with bound box check
	// 3) Do NOT intersect if both ends of lines do not lie on same side w.r.t each other
	// 4) Rest should intersect?
	
	BBox bound1 = this->boundBox();
	BBox bound2= s.boundBox();
	
	if( !doBoundsOverlap(bound1, bound2) )
		return false; // No intersection:  S1 and S2 are not within each other bounds
	
	// Need to branch based on what type of shape is coming in.
	
	// Not using dynamic_cats a sit is usually much slower. typeid is faster.	
	//const Line* const pL2 = dynamic_cast<const Line* const>(&s);
	if( typeid(s) == typeid(Polygon) ) { // Line to Polygon intersection
		const Polygon* const pP = static_cast<const Polygon* const>(&s);
		assert(pP); 
		
		// Check if any of the line segment point is inside the polygon.
		if( pP->isInside(first) || pP->isInside(second) ) {
			return true;
		}
		
		// Do proper intersection check
		for(int i=0; i < pP->mVertices.size(); i++) {
			Line l1(pP->mVertices[i], pP->mVertices[(i+1) % pP->mVertices.size()]);
			if( this->doesIntersect(l1) )
				return true;
		}
		
		// Checked al llines of polygon, no intersections found!
		return false;		
	}
	else if ( typeid(s) != typeid(Line) ) {
		return false;
	}
	
	// Line to Line intersection
	const Line* const pL2 = static_cast<const Line* const>(&s);
	assert( pL2 );
	if( !pL2 ) {
		return false;
		}
		
	// Edge cases : All of them are intersections
	// 1) Two lines collinear (overlap)
	// 2) Two lines share one of end points
	// 3) One of end point of one line lies within other line	
	Vec v1f2f(first, pL2->first), v1f2s(first, pL2->second);
	if( v1f2f.cross(v1f2s).magsqr() == 0 || haveCommonEndPoint(this, pL2) || 
	    this->isInside(pL2->first) || this->isInside(pL2->second) )
		return true; // Intersection: They must be overlapping since boundbox do overlap
	
	// Edge cases handled, we no wjust check if any of line lies completely on side of each other	
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
	// Check if given point is within line segment
	if( !isPointInsideBoundBox(this->boundBox(), p) )
		return false;
	
	// Check if point is along the line segment.
	// If point vector from line's one end is collinear with line,
	// point must be with line segment since point is already
	// in bounds of line.
	Vec vL(first, second), vP(first, p);
	return (vL.cross(vP) == 0);	
}

double Line::area(void) const {
	return 0;
}
	
bool Line::isConvex(void) const {
	return true;
} 

std::ostream& Line::print(std::ostream &o) const {
  return (o << "Line : {" << first << ", " << second << "}");
}

// Polygon

Polygon::Polygon(const std::vector<Point>& vertices) : mVertices(vertices) {}

BBox Polygon::boundBox(void) const {
	BBox bound = {{0, 0, 0}, {0, 0, 0}};
	
	if( mVertices.size() < 1 )
		return bound;
		
	bound.ll = mVertices[0];
	bound.ur = mVertices[0];
	
	for(int i=1; i < mVertices.size(); i++) {
		if( bound.ll.x > mVertices[i].x ) bound.ll.x = mVertices[i].x;
		if( bound.ll.y > mVertices[i].y ) bound.ll.y = mVertices[i].y;
		if( bound.ll.z > mVertices[i].z ) bound.ll.z = mVertices[i].z;
		if( bound.ur.x < mVertices[i].x ) bound.ur.x = mVertices[i].x;
		if( bound.ur.y < mVertices[i].y ) bound.ur.y = mVertices[i].y;
		if( bound.ur.z < mVertices[i].z ) bound.ur.z = mVertices[i].z;
	}
		
	return bound;
}

// Check if given two polygons intersect
// It wil lonly return true if two of them actually intersect (i.e. at least one line intersects)
// Overlap: It will return false in case one polygon is contained in other!

bool doPolygonsIntersect(const Polygon& poly1, const Polygon& poly2) {
	// Check if line segment of each polygon intersect with each other
	// Need to be careful to pick th elast connecting last point to first point.
	// (i+1) % N returns 0 for last point, allowing us to pick last line.
	
	const int N = poly1.mVertices.size();
	const int M = poly2.mVertices.size();	
	for(int i=0; i < N; i++)	{
		Line l1(poly1.mVertices[i], poly1.mVertices[(i+1) % N]);
		for(int j=0; j < M; j++) {
			Line l2(poly2.mVertices[j], poly2.mVertices[(j+1) % M]);
			if( l1.doesIntersect(l2) )
				return true;
		}	
	}

	// If we have reached so far, polygons are not intersecting
	// We ignore the case when one polygon is completely inside other
	return false;
	
} 

// Checks if any vertex of first polygon lies within second polygon
bool doesPoly1VertexLieInPoly2(const Polygon& poly1, const Polygon& poly2) {
	// Check if any point of first polygon lies inside second one
	
	for(int i=0; i < poly1.mVertices.size(); i++)	{
		if( poly2.isInside(poly1.mVertices[i]) ) 
			return true;
	}

	// If we have reached so far, no vertex of first polygon lies within second one!
	return false;	
} 

// Polygon to Shape intersection
bool Polygon::doesIntersect(const Shape& s) const {
	// Polygon to Polygon intersection Algo:
	// 1) Do NOT intersect if their bound boxes do not overlap
	// 2) Check if any edge of first shape intersects any edge of second shape
	//    We short circuit with "Intersect" as soon as we find the first intersection!
	// 3) If we have not found any intersection amongtheir edges, these shapes do not intersect.
	
	
	// Need to branch based on what kind of shape is coming in
	
	// Using typeid instead of dynamic cast as it is significantly faster.
	//const Polygon* const pS2 = dynamic_cast<const Polygon* const>(&s);	
	if( typeid(s) == typeid(Line) ) { // Polygon to line intersection
		const Line* const pL = static_cast<const Line* const>(&s);
		assert(pL);
	
		return pL->doesIntersect(*this);
	}
	else if( typeid(s) != typeid(Polygon) ) {
		return false;
	}
	
	// Polygon to Polygon intersection
	const Polygon* const pS2 = static_cast<const Polygon* const>(&s);
	assert( pS2 );
	if( !pS2 ) {
		return false;
		}
				
	// Bound check for quick rejection
	BBox bound1 = this->boundBox();
	BBox bound2= pS2->boundBox();
	
	if( !doBoundsOverlap(bound1, bound2) )
		return false; // No intersection:  S1 and S2 are not within each other bounds

	// Check for intersections including complete overlap
		
	if( mVertices.size() < 3 || pS2->mVertices.size() < 3 ) // Expecting at least a triangle
		return false;
		
	// We check for overlap, trying to see if one polygon vertex lies in another one
	// In case first one does not work out, we try other way around.
	if( doesPoly1VertexLieInPoly2(*this, *pS2) || doesPoly1VertexLieInPoly2(*pS2, *this) )
		return true;
		
	// Let us see if any line of one polygon intersects another one
	return doPolygonsIntersect(*this, *pS2);	
}

// Edge case: when ray hits either a vertex or collinear with an entire edge
// For an intersecting edge whose at least one vertex is on th horizontal ray (i.e. verteex edge),
// check if other vetex is below the ray. For crossing algorithm we will only count
// edges which are below the ray.
bool doesVertexEdgeLiesBelowRay(const Line& edge, const Line& ray) {
	if( edge.first.y == ray.first.y )
		return edge.second.y < ray.first.y; // count only if other vertex is below ray
	else {
		if( edge.second.y == ray.first.y )
			return edge.first.y < ray.first.y; // count omly if other vertex is below ray
		else
			return true; // count the intersection as no vertex is in path of horizontal ray
	}	
}

bool Polygon::isInside(const Point& p) const {
	BBox bound = this->boundBox();
	if( !isPointInsideBoundBox(bound, p) )
		return false;
		
	// Inside polygon (N sides) test using crossing method 
	const int N = mVertices.size();
	
	if( N < 3 ) // expecting at least a triangle
		return false;
			
	// Cast an horizontal ray in +X direction from point and count number of intersections
	// Number of intersections ODD => point is inside
	
	Line ray(p, Point({bound.ur.x+2, p.y, p.z}));
	int nIntersections = 0;
	
	for(int i=0; i < N; i++) {
		Line edge(mVertices[i], mVertices[(i+1) % N]);
		// We can use fast ray to line intersection instead of ray.doesIntersect(edge)
		if( doesRayIntersectLine(ray, edge) && doesVertexEdgeLiesBelowRay(edge, ray) )
			nIntersections++;
	}
	
	return isOdd(nIntersections);
}

double Polygon::area(void) const {
	// Area of a general 3D planar polygon having normal unit vector n'
	// Signed Area = 1/2 * (Normal . (SUM of cross product of its vertices))
	//
	// Area = 1/2 * (n'.SUM(Vi * Vi+1)), where i=0, n-1, and Vi is "i"th vertex as vector
	// Need to be careful for last vertex as its next vertex needs to be 0th vertex. We
	// can use "i%N" to take care of this.
	//
	// Positive sign of area means, polygon vertices are ordered anticlockwise.
	// Actual area is absolute value of this signed area.
	//
	// For normal 2D polygon, we can use faster computation:
	// Area = 1/2 * SUM( xi * (yi+1 - yi-1)), where i=1, n, and xi and yi belog to Vi vertex
	// To avoid using "i%N" to take care of last vertex, we can either temporarily insert
	// vertex v0 at the end or just do the last vetex computation outside the loop.
	
	double area = 0;
	const int N = mVertices.size();
	
	if( N < 3 ) // Expecting at least a triangle
		return 0;
	
	// Find sum of cross product of polygon vertices
	Vec v(0, 0, 0), normal(0, 0, 1);
	for(int i=0; i < N; i++) {
		Vec v0(mVertices[i]), v1(mVertices[(i+1) % N]);
		v += v0.cross(v1);			
	}
	
	// Calculate unit normal vector to the plane of planar polygon
	// We could have easily computed normal by taking cross product 
	// of first two edges, but these could be collinear!
	for(int i=0; i < N; i++) {
		Vec e0(mVertices[i], mVertices[(i+1) % N]), e1(mVertices[(i+1) % N], mVertices[(i+2) % N]);
		Vec n(e0.cross(e1));
		
		// We have the normal as soon as we find two consecutive edges which are not collinear
		if( n.magsqr() != 0 ) {
			normal = n.unit();
			break;
		}
	}
	
	// Area = 1/2 * Normal . SUM(cross product of its vertices)
	area = 0.5 * (normal.dot(v));
	
	// Return absolute value
	// We could have returned signed area, which can be used to find orientation of polygon
	// positive sign means polygon vertices are ordered anticlockwise.	
	
	return ABS(area);
}
	
bool Polygon::isConvex(void) const {
	// Assuming given polygon is not self-intersecting, we go round the polygon edges
	// seeing if winding (turn) is consistent as we move from one edge to another edge. If we
	// that any of these turns is not consistent with others, we have a concavity in the
	// polygon, ther epolygonis concave. If we find all turns are consitent with each other
	// our polygon is convex.
	const int N = mVertices.size();
	
	// Triangles are always convex!
	if( N < 4 ) {
		return true;
	}
		
	// Check turns at each corner of polygon
	// Need to be careful to pick the last edge and last turn. 
	// (i+1) % N returns 0 for last point, allowing us to pick last edge and first edge
	// for last turn.
	bool turnAntiClock = Vec(mVertices[0], mVertices[1]).isAntiClock(Vec(mVertices[1], mVertices[2]));
	for(int i=1; i < N; i++)	{
		Vec e1(mVertices[i], mVertices[(i+1) % N]);
		Vec e2(mVertices[(i+1) % N], mVertices[(i+2) % N]);		
		if( e1.isAntiClock(e2) != turnAntiClock ) // Check if this turn is same as other turns
			return false;	// Found first turn which is not consistent, polygon is concave
	}
	
	// All turns are consistent (go in same direction), our polygon is convex
	return true;
} 

std::ostream& Polygon::print(std::ostream &o) const {
  o << "Polygon: {";
  for(int i=0; i < mVertices.size(); i++) {
  	o << mVertices[i] << ", ";
  }
  o << "}";
  
  return (o);
}



