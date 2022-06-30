// Experiment with geometry
// Lekhraj

#include <limits>
#include <cmath>
#include <iostream>
#include <list>
#include <algorithm>
#include <numeric>

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

// Absolute value
Vec Vec::abs(void) const {
	Vec tmp(*this);
	tmp.x = ABS(tmp.x); tmp.y = ABS(tmp.y); tmp.z = ABS(tmp.z);
	
	return tmp;
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

// Simple minmax bound
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

// Check if given line and polygon intersect
// It wil lonly return true if two of them actually intersect (i.e. at least one line intersects)
// Overlap: It will return false in case line is contained in polygon!

bool doesLinePolygonIntersect(const Line& line, const Polygon& poly) {
	// Check if line segment intersects any of polygon edges. 
	// Need to be careful to pick the last connecting last point to first point.
	// (i+1) % N returns 0 for last point, allowing us to pick last edge.

	for(int i=0; i < poly.mVertices.size(); i++) {
		Line l1(poly.mVertices[i], poly.mVertices[(i+1) % poly.mVertices.size()]);
		if( line.doesIntersect(l1) )
			return true;
	}
	
	// Checked all edges of polygon, no intersection found!
	return false;
}		
	
// Line to Line intersection: Overlaps and concident end points are treated as intersections
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
		// We could have checked mid-point of line segment as well if check on end-points fail!
		if( pP->isInside(first) || pP->isInside(second) ) {
			return true;
		}
		
		// Do proper intersection check as line may be crossing polygon
		return doesLinePolygonIntersect(*this, *pP);
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

// Is given point within line segment? End points are treated as within.
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

double Line::Area(void) const {
	return 0;
}

double Line::signedArea(void) const {
	return 0;
}
	
bool Line::isConvex(void) const {
	return true;
} 

// Comparison operator, returns true if line end points are same in any order
bool Line::operator==(const Shape &rhsShape) const {
	
	if( typeid(rhsShape) != typeid(Line) ) 
			return false;
	
	// Line to Line comparsion
	const Line* const pL2 = static_cast<const Line* const>(&rhsShape);
	assert( pL2 );
	if( !pL2 ) {
		return false;
		}
	
	// Lines are same if their end points are same in any order
	return ((first == pL2->first && second == pL2->second) || (second == pL2->first && first == pL2->second));			
}

// Check if two lines are congruent (i.e. same geometrically)
bool Line::isCongruent(const Shape &rhsShape) const {
	// Two lines will be congruent if they are of same length! As they can be easily transformed into another
	// with affine transformation (translation and rotation).
	
	if( typeid(rhsShape) != typeid(Line) ) 
			return false;
	
	// Line to Line comparsion
	const Line* const pL2 = static_cast<const Line* const>(&rhsShape);
	assert( pL2 );
	if( !pL2 ) {
		return false;
		}
	
	// Lines are congruent if their lengths are same. 
	return (Vec(first, second).magsqr() == Vec(pL2->first, pL2->second).magsqr());			
}

std::ostream& Line::print(std::ostream &o) const {
  return (o << "Line : {" << first << ", " << second << "}");
}

// Polygon

// Return left bottommost vertex
int getBottomMostVertex(const std::vector<Point>& vertices) {
	int iBottom = 0;
	for(int i=1; i < vertices.size(); i++) { // in case of same bottom, pick leftmost
		if( (vertices[i].y < vertices[iBottom].y) || 
		    (vertices[i].y == vertices[iBottom].y && vertices[i].x < vertices[iBottom].x) ) 
			iBottom = i;
	}
	
	return iBottom;
}

// Return true if p1 comes first (before p2) when going from x-axis in anti-clockwise direction
// 
int compareAngle(const Point& p0, const Point& p1, const Point& p2) {
	// We are checking if edge p0-p1 is to the left of edge p0-p2 in anticlockwise direction
	Vec e1(p0, p1), e2(p0, p2);
	double angleMeasure = e1.cross(e2).z;
	if( angleMeasure == 0 ) // In case of same angle, vertex closer to bottommost point is picked
		return e1.magsqr() <= e2.magsqr();
		
	return angleMeasure >= 0; // vertex which is left to other (i.e. form anticlock turn w.r.t other vertex edge)
}

// Polygon vertices can be given either in anti-clockwise order or random order (default is CCW order)
// For random vertices order, a simple closed polygon is constructed in increasing angle from leftmost bottommost vertex
Polygon::Polygon(const std::vector<Point>& vertices, const bool randomOrder) : mVertices(vertices) {
	if( randomOrder ) {
		// Construct a simple closed path polygon 
		// Algo:
		// 1) Find the bottommost point
		// 2) For every other vertex, calculate its angle from bottommost vertex w.r.t to x-axis
		// 3) Connect the vertices in the increasing order of their angle (natural anti-clock order)
		//    -For same angle pick them in increasing order of their distance from bottom most vertex
		//	  -We do not need actual angle, just measure of orientation (like cross product w.r.t bottommost point) would do
		
		// Fidn leftmost bottom most vertex		
		int iBottom = getBottomMostVertex(mVertices);
		Point p0(mVertices[iBottom]);
		std::swap(mVertices[0], mVertices[iBottom]);
		
		// Sort the rest of vertex array in order of their angle from x-axis
		std::sort(mVertices.begin()+1, mVertices.end(), 
		          [&p0](const Point& p1, const Point& p2){ return compareAngle(p0, p1, p2);});
		
	}
}

typedef std::pair<Point, Point> EDGE;

std::list<EDGE>::iterator findEDGE(std::list<EDGE> vList, EDGE& edge) {
	for(auto it=vList.begin(); it != vList.end(); it++) {
		if( (*it) == edge )
			return it;
	}
	
	return vList.end();
}

// Construct the shape from list of connected triangles (triangles are in anti-clockwise order)
Polygon::Polygon(const std::vector<std::tuple<Point,Point,Point>>& vecTri) {
	// Algo to construct shape from list of connected triangles
	// 1) Construct a list of exterior edges from list of triangles (remove all interior edges)
	//    For each edge in triangle:
	//			a) Remove the edge if already exists (interior edge)
	//			b) Insert the edge in right place (linking with previous edge)
	// 2) Spit out the shape by traversing the edge list
	
	// We are using list here but could have used map for better performance of find operation.
	
	// Make the exterior edge list
	Point v[3];
    EDGE pairPoints;
	std::list<EDGE> vList;
	for(int i=0; i < vecTri.size(); i++) {
		std::tie(v[0], v[1], v[2]) = vecTri[i];
		for( int j = 0; j < 3; j++ ) {
		    pairPoints = std::make_pair(v[(j+1)%3], v[j]);
			auto it = std::find(vList.begin(), vList.end(), pairPoints);
			//auto it = findEDGE(vList, pairPoints);
			if( it != vList.end() ) { // interior edge found
				vList.erase(it); // remove the edge
			}
			else { // Insert the edge in right place 
				auto itInsert = std::find_if(vList.begin(), vList.end(), [&v, j](EDGE& e){return e.second == v[j];});
				if( itInsert != vList.end() ) 
					vList.insert(std::next(itInsert), std::make_pair(v[j], v[(j+1)%3]));
				else
					vList.push_back(std::make_pair(v[j], v[(j+1)%3]));			
			}			
		}		
	}
	
	// Put elements in connected order if they are not already
	for(auto it = vList.begin(); it != vList.end(); it++) {
		auto itNext = std::next(it);
		if( it != vList.end() && !((*it).second == (*itNext).first) ) {
			auto itInsert = std::find_if(it, vList.end(), [&it](EDGE& e){return e.first == (*it).second;});
			if( itInsert != vList.end() ) { // Move the element to its right connected position
				vList.insert(itNext, (*itInsert));
				vList.erase(itInsert);
			}
		}
	}
	
	// Populate our polygon vertices
	// There could be multiple loops if some interior triangles are missing! It can be used to detect
	// holes in shape!
	// Here we return the first loop, which works when there are no missing interior triangles.
	for(auto& e: vList) {
		mVertices.push_back(e.first);
		if( e.second == vList.front().first ) // reached back to start of loop
			break;
	}
	
	// Done	
}

// Simple minmax bound box

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
	for(int i=0; i < N; i++)	{
		Line l1(poly1.mVertices[i], poly1.mVertices[(i+1) % N]);
		if( doesLinePolygonIntersect(l1, poly2) )
				return true;	
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

// Polygon to Shape intersection, we do check for full overlap (one inside another) here!
bool Polygon::doesIntersect(const Shape& s) const {
	// Polygon to Polygon intersection Algo:
	// 1) Do NOT intersect if their bound boxes do not overlap
	// 2) Check for overlaps pure edge to edge intersection will miss out the case when one is inside another!)
	// 3) Check if any edge of first shape intersects any edge of second shape
	//    We short circuit with "Intersect" as soon as we find the first intersection!
	// 4) If we have not found any intersection among their edges, these shapes do not intersect.
	
	
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

// Is given point inside the polygon? Boundary is treated as inside.

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

// Signed Area: +ive for anti-clockwise ordered polygon, -ive for clockwise ordered polygon
double Polygon::signedArea(void) const {
	// Signed area of a general 3D planar polygon having normal unit vector n'
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
			// We want positive unit normal to get proper sign otherwise sign of "v" and sign of "normal"
			// cancel out, not giving us true signed area!
			normal = n.unit().abs();
			
			break;
		}
	}
	
	// Area = 1/2 * Normal . SUM(cross product of its vertices)
	area = 0.5 * (normal.dot(v));
	
	// Returning sign which can be used to find orientation of polygon
	// positive sign means polygon vertices are ordered anticlockwise.	
	
	return area;
}

// Return absolute value of signed area
double Polygon::Area(void) const {
	return ABS(signedArea());
}

// Is the convex (i.e. no cavity in it)?	
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

// Triangulation

struct VertexStatus {
	int idx; // index into polygon vertices
	bool isConcave; // Is the corner at this vertex concave?
};

// Check if a corner v0-v1-v2 (edge v0v1 is first edge, edge v1-v2 is second edge in anti-clockwise order)
// is concave (i.e. repersents a cavity in polygon).
//
// Assuming polygon vertices are ordered in anti-clockwise (otherwise check will be reverse), we check
// if first edge vector to second edg vector goes in clockwise direction. For regular corner (convex),
// it should be anticlockwise direction (i fgiven polygon vertices are anticlockwise oredered).

bool isCornerConcave(const Point& v0, const Point& v1, const Point& v2) {
	return !(Vec(v0, v1).isAntiClock(Vec(v1, v2)));
}

// Circular walking of list
std::list<VertexStatus>::iterator prevElem(std::list<VertexStatus>& list, std::list<VertexStatus>::iterator it) {
	if( it == list.begin() )
		return std::prev(list.end());
	else
		return std::prev(it); 
}

std::list<VertexStatus>::iterator nextElem(std::list<VertexStatus>& list, std::list<VertexStatus>::iterator it) {
	if( it == std::prev(list.end()) || it == list.end() )
		return list.begin();
	else
		return std::next(it); 
}

bool nextConcaveCorner(const std::list<VertexStatus>& list, std::list<VertexStatus>::iterator& it) {
	while( it != list.end() && !(*it).isConcave ) {
		it++;
	}
	
	return (it != list.end() && (*it).isConcave);
}

bool isConcaveVertex(const std::list<VertexStatus>& vList, const int v) {
	for( auto it = vList.begin(); it != vList.end(); it++ ) {
		if( (*it).idx == v )
			return (*it).isConcave;
	}
	
	return false;	
}

// Check if this a valid triangle, therefore no polygon vertex is inside given triangle
bool validTriangle(const Polygon& poly, const std::list<VertexStatus>& vList,
                   const int v0, const int v1, const int v2) {
	Polygon tri({poly.mVertices[v0], poly.mVertices[v1], poly.mVertices[v2]});
	for(int i=0; i < poly.mVertices.size(); i++) {
		if( i != v0  && i != v1 && i != v2 && isConcaveVertex(vList, i) )  { // Do not want to check with the triangle itself!
			if( tri.isInside(poly.mVertices[i]) ) 
				return false;
		}
	}
	
	// No vertex found inside given triangle
	return true;			
	
}

void updateConcavity(const Polygon& poly, std::list<VertexStatus>& vList, std::list<VertexStatus>::iterator& it) {
	// Update concavity status of our current concave vertex
	auto itPrev1 = prevElem(vList, it);
	auto itNext1 = nextElem(vList, it);
	(*it).isConcave = isCornerConcave(poly.mVertices[(*itPrev1).idx], poly.mVertices[(*it).idx], 
									  poly.mVertices[(*itNext1).idx]);	
}

// Generate triangles from the shape (tessellate/triangulate shape into triangles)
int Polygon::genTriangles(std::vector<std::tuple<int,int,int>>& vecTri) const {
	// Tesslate any arbitary polygon (concave or convex) into triangles.
	// Assumption: Polygon vertices are given in anticlockwise order and it si not self-intersecting.
	// Our approach:
	// 1) Classify all corners (vertex) as either having concavity or not
	// 2) While there is an concave corner left, do
	//    	Chop-off concave corner:
	//			Spit out a triangle with its previous and next vertex if possible
	//				Only attempt to spit out if it is a valid triangle
	//				Remove the spitted out previous and next vertex
	//		Update concavity status of current vertex and new adjacent vertices (after previous adjacent vertex removals)
	//		If current vertex is no longer concave, move to next concave corner
	// 3) After we have removed all concavity, we will be left with convex polygon
	// 4) Spit out fanned triangles for the remaining convex polygon
	
	const int N = mVertices.size();
	vecTri.clear();
	
	if( N < 3 ) {
		return 0; 
	}
	else if( N == 3 ) { // Do not need to triangulate an already triangle!
		vecTri.push_back(std::make_tuple(0, 1, 2));
		return 1;
	}
	else if( this->isConvex() ) { // Pure convex shape
		// Create a fan of triangles with first vertex
		for(int i = 1; i < N-1; i++) {
			vecTri.push_back(std::make_tuple(0, i, i+1));
		}
		
		return vecTri.size();		
	}
	
	// Initialize each vertex with its concavity info
	std::list<VertexStatus> vList;
	VertexStatus vS;
	vS.idx = 0;
	vS.isConcave = isCornerConcave(mVertices[N-1], mVertices[0], mVertices[1]);
	vList.push_back(vS);
	for(int i=1; i < N; i++) {
		vS.idx = i;
		vS.isConcave = isCornerConcave(mVertices[i-1], mVertices[i], mVertices[(i+1)%N]);
		vList.push_back(vS);
		}
		
	// Start chopping of concave corners
	auto it = vList.begin();
	auto itPrev1=it, itNext1=it, itPrev2=it, itNext2=it;
	int i, j, k;
	while( nextConcaveCorner(vList, it) && vecTri.size() < N )  {
	
		// Spit out a triangle with pervious vertex if we can
		itPrev1 = prevElem(vList, it);
		itPrev2 = prevElem(vList, itPrev1);
		if( !((*itPrev1).isConcave) ) {
			i = (*itPrev2).idx; j = (*itPrev1).idx; k = (*it).idx;
			if( validTriangle(*this, vList, i, j, k) ) {
				vecTri.push_back(std::make_tuple(i, j, k));	
						
				// Remove the previous vertex
				vList.erase(itPrev1);
			}			
		}
		
		// Spit out a triangle with next vertex if we can
		itNext1 = nextElem(vList, it);
		itNext2 = nextElem(vList, itNext1);
		if( !((*itNext1).isConcave) ) {
			i = (*it).idx; j = (*itNext1).idx; k = (*itNext2).idx;
			if( validTriangle(*this, vList, i, j, k) ) {		
				vecTri.push_back(std::make_tuple(i, j, k));
			
				// Remove the next vertex
				vList.erase(itNext1);	
			}		
		}
		
		// Update concavity status of our current concave vertex and its two neigbhors
		updateConcavity(*this, vList, it);
		
		itPrev1 = prevElem(vList, it);
		updateConcavity(*this, vList, itPrev1);
				
		itNext1 = nextElem(vList, it);
		updateConcavity(*this, vList, itNext1);
		
		// Degeneracy check
		if( itPrev1 == itNext1 ) 
			break;
	}
	
	// We should only have convex polygon remaining now
	// Just create a fan of triangles starting from first vertex
	
	auto it0 = vList.begin();
	auto it1 = std::next(it0);
	for( auto it = it1; it != vList.end(); it++ ) {
	    it1 = std::next(it);
		if( it1 != vList.end() ) {
			assert( !((*it0).isConcave) && !((*it).isConcave) && !((*it1).isConcave) );
			vecTri.push_back(std::make_tuple((*it0).idx, (*it).idx, (*it1).idx));
			}		
	}
			

	// Done with triangulating our polygon, phew!
	return vecTri.size();		
}  

// Comparison operator, returns true if polygon is exactly same shape (0..k..n == k..n..0)
bool Polygon::operator==(const Shape &rhsShape) const {
	
	if( typeid(rhsShape) != typeid(Polygon) ) 
			return false;
	
	// Polygon to Polygon comparsion
	const Polygon* const pS2 = static_cast<const Polygon* const>(&rhsShape);
	assert( pS2 );
	if( !pS2 ) {
		return false;
		}
	
	// Trivial case of exact same vertices
	if( mVertices == pS2->mVertices )
		return true;
	else if( mVertices.size() != pS2->mVertices.size() ) // Even their vertex count do not match!
		return false;
		
	// Need to check if order of vertices in both the polygons is basically same
	// We first check where does first vertex of second polygon lie. Then we cyclically
	// keep checking until we reach back to same vertex.
	auto it = std::find(mVertices.begin(), mVertices.end(), pS2->mVertices[0]);
	if( it == mVertices.end() ) // Did not even find match for first vertex of second polygon
		return false;
		
	for( int i=1; i < pS2->mVertices.size(); i++ ) {
		// Advance first polygon vertices in circular fashion
		it++;
		if( it == mVertices.end() )
			it = mVertices.begin();
		
		// Keep comparing until there is match	
		if( !((*it) == pS2->mVertices[i]) ) // mismatch, bail out, polygons are not same!
			return false;		
	}
	// At this state, Next element in first polygon should be same as first element of second polygon
	assert( mVertices[((it - mVertices.begin()) + 1)%mVertices.size()] == pS2->mVertices[0] ); // back to original vertex	
		
	return true;		
}

// Generate corner measure of given polygon
// Our corner measure should be proxy to comparing two edges length and their interior angle
void genCornerPoly(const Polygon& poly, std::vector<Point>& cornerPoly) {
	const int N = poly.mVertices.size();
	Point cornerMeasure;
	cornerPoly.clear();
	for(int i=1; i <= N; i++) {
		Vec e1(poly.mVertices[i%N], poly.mVertices[(i+1) % N]);
		Vec e2(poly.mVertices[i%N], poly.mVertices[(i-1) % N]);
		// Our corner measure is cross product between two edges and their individual magnitude
		// We are storing edge magnitude in sorterd order to make it work whichever order they come in
		cornerMeasure.z = e1.cross(e2).z; // e1.dot(e2) is not effective when angle is 90 deg between them!
		cornerMeasure.x = e1.magsqr();
		cornerMeasure.y = e2.magsqr();
		if( cornerMeasure.x > cornerMeasure.y ) // swap them, keep minimum first
			std::swap(cornerMeasure.x, cornerMeasure.y);
			
		cornerPoly.push_back( cornerMeasure );			
	}	
}

// Compare two corner measures in either anti-clock wise (+1) or clock-wise (-1) direction
bool compareCyclic(const std::vector<Point>& cornerPoly1, const std::vector<Point>& cornerPoly2, const int iDir = 1) {
	if( cornerPoly1 == cornerPoly2 ) // trivial case when they are same!
		return true;
		
	// Find the first element of second poly in first 
	auto it = std::find(cornerPoly1.begin(), cornerPoly1.end(), cornerPoly2[0]);
	if( it == cornerPoly1.end() ) // could not even find match for first measure anywhere in second!
		return false;
	
	// We may have duplicates in poly1, keep looping until we exhaust the poly1
	bool bNext = true;
	int iFirst = 1;
	if( iDir < 0 )
		iFirst = cornerPoly2.size()-1;
	
		
	for( ; bNext && it != cornerPoly1.end(); it = std::find(++it, cornerPoly1.end(), cornerPoly2[0]) ) {			
		// Start comparing cyclically from next element
		bNext = false;	
		for(int i=iFirst, j = (it - cornerPoly1.begin())+1; (i < cornerPoly2.size() && i > 0) ; i += iDir, j++ ) {
			// Walk cyclically in given direction
			if( !(cornerPoly2[i] == cornerPoly1[j%cornerPoly1.size()]) ) {
				bNext = true;
				break; 
			}
		}
	}
	
	// Did we match?
	return !bNext;
}

// Check if two shapes are congruent (i.e. same geometrically)
bool Polygon::isCongruent(const Shape &rhsShape) const {
	// Two polygons are congruent if:
	// a) Same number of points (i.e. same number of edges)
	// b) Their sides are same size
	// c) Their inetrior angles are same
	
	// For our purpose we do not need to calculate exact size or angle.
	// We can use some sort of corner measure to compare if two corners are same.
	
	if( typeid(rhsShape) != typeid(Polygon) ) 
			return false;
	
	// Polygon to Polygon comparsion
	const Polygon* const pS2 = static_cast<const Polygon* const>(&rhsShape);
	assert( pS2 );
	if( !pS2 ) {
		return false;
		}
	
	// Trivial case of exact same vertices
	if( mVertices == pS2->mVertices )
		return true;
	else if( mVertices.size() != pS2->mVertices.size() ) // Even their vertex count do not match!
		return false;
		
	// Need to do actual congruency test
	// We will build a corner measure and then compare if these measures
	// are same in some cyclic order.
	
	// Generate corner measure for both polygons
	std::vector<Point> cornerPoly1, cornerPoly2;
	genCornerPoly(*this, cornerPoly1);
	genCornerPoly(*pS2, cornerPoly2);
	
	// Compare corner measures in cyclic fashion
	if( compareCyclic(cornerPoly1, cornerPoly2, 1) )
		return true;
	else if( compareCyclic(cornerPoly1, cornerPoly2, -1) )
		return true;
	
	// No match found in either direction; two polygons are not congruent!
	return false;
}

struct Space {
	double x,y;
	double w,h;
};

// Checks if this rect will fit in given space
bool fitsSpace(std::vector<Space>& spaces, const int i, const Rect& rect) {
		return (rect.w <= spaces[i].w && rect.h <= spaces[i].h);	
}


// Checks if this rect will end up creating new space
bool newSpace(std::vector<Space>& spaces, const int i, const Rect& rect) {
		return (rect.w != spaces[i].w && rect.h != spaces[i].h);	
}

// Checks if this rect will end up extending existing boundaries
bool extendsBounds(std::vector<Space>& spaces, const int i, const Rect& rect, const double W, const double H) {
		return 	(spaces[i].x + rect.w > W || spaces[i].y + rect.h > H);
} 

void splitSpace(std::vector<Space>& spaces, const int i, const Rect& rect) {
		// Adjust spaces, split if necessary
		Space space = spaces[i];
		
		if( rect.w == space.w && rect.h == space.h ) { // Fits perfectly, remove it
			space = spaces.back();
    		spaces.pop_back();				
		}
		else if( rect.h == space.h ) { // Fits horizontally, adjust space
		    // |-------|---------------|
	        // |  rect | updated space |
	        // |_______|_______________|
	        space.x += rect.w;
	        space.w -= rect.w;				
		}
		else if( rect.w == space.w ) { // Fits vertically, adjust space 
	        // |---------------|
	        // |      rect     |
	        // |_______________|
	        // | updated space |
	        // |_______________|
	        space.y += rect.h;
    		space.h -= rect.h;							
		}
		else { // Need to split
	        // |-------|-----------|
	        // |  rect | new space |
	        // |_______|___________|
	        // | updated space     |
	        // |___________________|
			
			// Add new space
	        spaces.push_back({x: space.x + rect.w, y: space.y, w: space.w - rect.w, h: rect.h});
			//space = spaces[i]; // Container element refernce can change after another elmeent is added!

			// Update existing space
	        space.y += rect.h;
	        space.h -= rect.h;				
		}
		
		// Update the original spaces
		spaces[i] = space;	
}

int findSpace(const std::vector<Space>& spaces, const Rect& rect) {
	// Look for space to fit given rectangle in reverse order in spaces (we find smallest first!)
	for(int i=spaces.size()-1; i >= 0; i--) {
		// Do we have fit?
		if( rect.w <= spaces[i].w && rect.h <= spaces[i].h ) { // Found the space				
			return i;
		}			
	}
	
	// If we could not find anything!
	return 0;			
}

int findLeftmostSpace(const std::vector<Space>& spaces, const Rect& rect) {
	// Look for space to fit given rectangle in reverse order in spaces (we find smallest first!)
	
	if( spaces.size() <= 1 ) // trivial, no need to search
		return 0;
		
	auto itFound = spaces.rbegin();
	auto itFirst = --(spaces.rend()); //reverse iterator end points to just past first element!
	auto itLeftFound = itFound;
	
	// Search in reverse
	itFound = std::find_if(itFound, spaces.rend(), 
							[&rect](const Space& s){ return (rect.w <= s.w && rect.h <= s.h); });
	itLeftFound = itFound;
	
	// See if we can find another match which is leftmost	
	while( itFound != itFirst && itFound != spaces.rend() ) {
		// Search in reverse: our next range start is past previous found: ++itFound (reverse iterator)
		itFound = std::find_if((++itFound), spaces.rend(), 
							[&rect](const Space& s){ return (rect.w <= s.w && rect.h <= s.h); });
							
		if( itFound != itFirst && (*itFound).x < (*itLeftFound).x ) { // Found a new lower left
			itLeftFound = itFound;			
		}									
	}
																				
	if( itLeftFound != spaces.rend() ) { // convert it into index from top
		return (itFirst - itLeftFound); // reverse iterator: They increment in reverse direction!
	}
	else {
		std::cout << "Error: Not found any place for rectangle {" << rect.w << ", " << rect.h << "}" << std::endl;
		return 0;
	}	
}

int findBottommostSpace(const std::vector<Space>& spaces, const Rect& rect) {
	// Look for space to fit given rectangle in order in spaces (we find smallest first!)
	
	if( spaces.size() <= 1 ) // trivial, no need to search
		return 0;
		
	auto itFound = spaces.begin();
	auto itFirst = spaces.begin();
	auto itLeftFound = itFound;
	
	// Search 
	itFound = std::find_if(itFound, spaces.end(), 
							[&rect](const Space& s){ return (rect.w <= s.w && rect.h <= s.h); });
	itLeftFound = itFound;
	
	// See if we can find another match which is leftmost	
	while( itFound != itFirst && itFound != spaces.end() ) {
		// Search in reverse: our next range start is past previous found: ++itFound
		itFound = std::find_if((++itFound), spaces.end(), 
							[&rect](const Space& s){ return (rect.w <= s.w && rect.h <= s.h); });
							
		if( itFound != itFirst && (*itFound).x < (*itLeftFound).x ) { // Found a new lower left
			itLeftFound = itFound;			
		}									
	}
																				
	if( itLeftFound != spaces.end() ) { // convert it into index from top
		return (itLeftFound - itFirst);
	} 
	else {
		std::cout << "Error: Not found any place for rectangle {" << rect.w << ", " << rect.h << "}" << std::endl;
		return 0;
	}	
}

void calculateBounds(double& W, double& H, const std::vector<std::tuple<int, Point, bool>>& packedRects, const std::vector<Rect>& rects) {
	// Our height is topmost space vertical position as it would be first space still remaining
	// since we chose large vertical space initially.
	// Our width is rightmost rect position
	// W = spaces[0].w; H = spaces[0].y;
	int i;
	Point loc;
	bool flipped;
	double rightEdge, bottomEdge;	
	W = 0; H = 0;
	// Find rightmost edge among our placed rects
	for(const auto& r: packedRects) {
		std::tie(i, loc, flipped) = r;
		rightEdge = (!flipped ? loc.x + rects[i].w : loc.x + rects[i].h);
		bottomEdge = (!flipped ? loc.y + rects[i].h : loc.y + rects[i].w);

		if( W < rightEdge ) W = rightEdge;
		if( H < rightEdge ) H = bottomEdge;
	}
}


// Pack given rectangles (w,h) into a compact rectangle sheet (input/output: W, H) and their positions
// (index of rect, x, y, flip status).
// Returns packing density (ratio of packed area to overall output sheet area)
double packRect(std::vector<std::tuple<int, Point, bool>>& packedRects, double& W, double& H,const std::vector<Rect>& rects, const bool flipRect) {
	// We use simple algo: partioning output space into bins by splitting appropriately.
	// 1. Sort the input rectangles in decreasing order of their heights
	// 2. Initialize output rect space to accomodate at least maximum width rect
	// 3. For each input rect in sorted rects: 
	//		-Place it in smallest space bin possible
	//      -If found space is perfect match, remove this bin
	//		-If found space matches rect height/width, adjust its free space accoridingly
	//		-Otherwise split the space (new and adjust current one)
	
	// Create a sorted index on given rectangles in decreasing order of their heights
	// for same heights rectangles, prefer larger width rectangle
    std::vector<int> sortedRectIdx(rects.size());
  	std::iota(sortedRectIdx.begin(), sortedRectIdx.end(), 0);
  	std::sort(sortedRectIdx.begin(), sortedRectIdx.end(),
       		  [&rects](const int i1, const int i2) 
			  {return (rects[i1].h == rects[i2].h ? rects[i1].w > rects[i2].w : rects[i1].h > rects[i2].h);});
			  
	// Find maximum width among given rectangles
	double maxW = 0, totalH = 0, area = 0, maxH = 0, sqrSide =0;
	for(const auto& r: rects) {
		if( maxW < r.w ) maxW = r.w;
		if( maxH < r.h ) maxH = r.h;
		totalH += r.h;
		area += r.w * r.h;
	}

    sqrSide = (int)(sqrt(area*1.1)+0.5);
	
	// Initialize space to go for squarish fit: side of square containing total area (bit larger)	
	if( W == 0 ) 
		W = MAX(maxW, sqrSide);
	else 
		W = MAX(maxW, W);

	if( H == 0 )
		H = totalH;
	else
		H = MAX(totalH, H);	
	
	Rect flippedRect, tRect;
	std::vector<Space> spaces;
	spaces.push_back({x:0, y:0, w:W, h:H});
	packedRects.clear();
	
	// Put in the first rect
	int r = sortedRectIdx[0], s, s1;
	bool flipped = false;
	packedRects.push_back(std::make_tuple(r, Point({spaces[0].x, spaces[0].y, 0}), false));
	W = spaces[0].x + rects[r].w; H = spaces[0].y + rects[r].h;	
	splitSpace(spaces, 0, rects[r]);			
	
	// Place rest of the rectangles
	for(int j=1; j < sortedRectIdx.size(); j++) {
		// Looking for place to fit this rect
		// We look in reverse order to first match with small space
		// int i = findSpace(spaces, rects[r]);
		r = sortedRectIdx[j];
		
		flippedRect = rects[r];
		flipped = false;
				
		s = findLeftmostSpace(spaces, flippedRect);	
		
		// We prefer orientation that does not create new space and extend bounds
		//if( flipRect && (newSpace(spaces, s, flippedRect) || flippedRect.w > flippedRect.h) ) {
		if( flipRect && (newSpace(spaces, s, flippedRect) || extendsBounds(spaces, s, flippedRect, W, H)) ) {
			tRect.w = flippedRect.h; tRect.h = flippedRect.w;
			s1 = findLeftmostSpace(spaces, tRect);				
			if( fitsSpace(spaces, s1, tRect) && (!extendsBounds(spaces, s1, tRect, W, H) || s1 != 0)  ) {
				flippedRect = tRect;
				flipped = !flipped;
				s = s1;
			}			
		}
		
		packedRects.push_back(std::make_tuple(r, Point({spaces[s].x, spaces[s].y, 0}), flipped));
		if( W < spaces[s].x + flippedRect.w ) W = spaces[s].x + flippedRect.w;
		if( H < spaces[s].y + flippedRect.h ) H = spaces[s].y + flippedRect.h;
						
		splitSpace(spaces, s, flippedRect);
		
	}						
	
	// Did we pack all rects?
	assert( packedRects.size() == rects.size() );
	
	
	// We are now tracking W,H live now!
	// calculateBounds(W, H, packedRects, rects);
		
	// return packing density: What % of output rectangle sheet (W*H) was packed with input rectangles
	return (area/(W*H));		
}

// Pack given rectangles (w,h) into a compact squarish sheet (output: W, H) and their positions
// (index of rect, x, y, flip status).
// Returns packing density (ratio of packed area to overall output sheet area)
double packRectSqr(std::vector<std::tuple<int, Point, bool>>& packedRects, double& W, double& H,const std::vector<Rect>& rects, const bool flipRect) {
	// We use simple algo: partioning output space into bins by splitting appropriately.
	// 1. Sort the input rectangles in decreasing order of their heights
	// 2. Initialize output rect space to accomodate at least maximum width rect
	// 3. For each input rect in sorted rects: 
	//		-Place it in smallest space bin possible
	//      -If found space is perfect match, remove this bin
	//		-If found space matches rect height/width, adjust its free space accoridingly
	//		-Otherwise split the space (new and adjust current one)
	
	// Create a sorted index on given rectangles in decreasing order of their heights
	// for same heights rectangles, prefer larger width rectangle
    std::vector<int> sortedRectIdx(rects.size());
  	std::iota(sortedRectIdx.begin(), sortedRectIdx.end(), 0);
  	std::sort(sortedRectIdx.begin(), sortedRectIdx.end(),
       		  [&rects](const int i1, const int i2) 
			  {return (rects[i1].h == rects[i2].h ? rects[i1].w > rects[i2].w : rects[i1].h > rects[i2].h);});
			  
	// Find maximum width among given rectangles
	double maxW = 0, totalH = 0, area = 0, totalW = 0;
	for(const auto& r: rects) {
		if( maxW < r.w ) maxW = r.w;
		totalH += r.h;
		area += r.w * r.h;
		totalW += r.w;
	}
	
	// Initialize space to go for squarish fit: side of square containing total area (bit larger)
	//W = MAX(maxW, (int)(sqrt(area*1.1)+0.5)); H = totalH;	
	//W = totalW; H = totalH;
	W = (int)(sqrt(totalW*totalH) + 0.5);
	H = W;	
	
	Rect flippedRect, tRect;
	std::vector<Space> spaces;
	spaces.push_back({x:0, y:0, w:W, h:H});
	packedRects.clear();
	
	// Put in the first rect
	int r = sortedRectIdx[0], s, s1;
	bool flipped = false;
	packedRects.push_back(std::make_tuple(r, Point({spaces[0].x, spaces[0].y, 0}), false));	
	W = rects[r].w; H = rects[r].h;
	splitSpace(spaces, 0, rects[r]);			
	
	// Place rest of the rectangles
	for(int j=1; j < sortedRectIdx.size(); j++) {
		// Looking for place to fit this rect
		// We look in squarish order to first match with small space
		// int i = findSpace(spaces, rects[r]);
		r = sortedRectIdx[j];
		
		flippedRect = rects[r];
		flipped = false;
		
		if( W + flippedRect.w <= H ) {	
			s = findLeftmostSpace(spaces, flippedRect);	
		}
		else {
			s = findBottommostSpace(spaces, flippedRect);				
		}
		
		// We prefer orientation that does not create new space and extend bounds
		if( flipRect && (newSpace(spaces, s, flippedRect) || extendsBounds(spaces, s, flippedRect, W, H)) ) {
			tRect.w = flippedRect.h; tRect.h = flippedRect.w;
			s1 = findLeftmostSpace(spaces, tRect);				
			if( fitsSpace(spaces, s1, tRect) && (!extendsBounds(spaces, s1, tRect, W, H) || s1 != 0)  ) {
				flippedRect = tRect;
				flipped = !flipped;
				s = s1;
			}			
		}
		
		packedRects.push_back(std::make_tuple(r, Point({spaces[s].x, spaces[s].y, 0}), flipped));
		if( W < spaces[s].x + flippedRect.w ) W = spaces[s].x + flippedRect.w;
		if( H < spaces[s].y + flippedRect.h ) H = spaces[s].y + flippedRect.h;
						
		splitSpace(spaces, s, flippedRect);
	}						
	
	// Did we pack all rects?
	assert( packedRects.size() == rects.size() );
	
	// We are now tracking W,H live now!
	// calculateBounds(W, H, packedRects, rects);
	
	// return packing density: What % of output rectangle sheet (W*H) was packed with input rectangles
	return (area/(W*H));		
}

double packRectFixedSize(std::vector<std::tuple<int, Point, bool>>& packedRects, double& W, double& H,const std::vector<Rect>& rects, const bool flipRect, std::vector<int>& sortedRectIdx, double area) {

	Rect flippedRect, tRect;
	std::vector<Space> spaces;
	spaces.push_back({x:0, y:0, w:W, h:H});
	packedRects.clear();
	
	// Make the first rect to have max height even if it needs to be flipped
	bool firstRectFlipped = false;
	if( flipRect ) {
		
		int m=0;
		double maxHeight = rects[sortedRectIdx[0]].h;
		for(int i=0; i < sortedRectIdx.size(); i++) {
			if( rects[sortedRectIdx[i]].w > maxHeight ) {
				m = i;
				maxHeight = rects[sortedRectIdx[i]].w;
				firstRectFlipped = true;
			}
		}
		
		// Swap with first rect with max width rect
		if( firstRectFlipped && m > 0) {
			std::swap(sortedRectIdx[0], sortedRectIdx[m]);
		}
	}
	
	// Put in the first rect
	int r = sortedRectIdx[0], s, s1;
	bool flipped = firstRectFlipped;
	packedRects.push_back(std::make_tuple(r, Point({spaces[0].x, spaces[0].y, 0}), firstRectFlipped));
	flippedRect = rects[r];
	if(firstRectFlipped) {
		flippedRect.w = rects[r].h;
		flippedRect.h = rects[r].w;
	}
	W = spaces[0].x + flippedRect.w; H = spaces[0].y + flippedRect.h;	
	splitSpace(spaces, 0, flippedRect);
	
	// Place rest of the rectangles
	for(int j=1; j < sortedRectIdx.size(); j++) {
		// Looking for place to fit this rect
		// We look in reverse order to first match with small space
		// int i = findSpace(spaces, rects[r]);
		r = sortedRectIdx[j];
		
		flippedRect = rects[r];
		flipped = false;
				
		s = findLeftmostSpace(spaces, flippedRect);	
		
		// We prefer orientation that does not create new space and extend bounds
		//if( flipRect && (newSpace(spaces, s, flippedRect) || flippedRect.w > flippedRect.h) ) {
		if( flipRect && (newSpace(spaces, s, flippedRect) || extendsBounds(spaces, s, flippedRect, W, H)) ) {
			tRect.w = flippedRect.h; tRect.h = flippedRect.w;
			s1 = findLeftmostSpace(spaces, tRect);				
			if( fitsSpace(spaces, s1, tRect) && (!extendsBounds(spaces, s1, tRect, W, H) || s1 != 0)  ) {
				flippedRect = tRect;
				flipped = !flipped;
				s = s1;
			}			
		}
		
		packedRects.push_back(std::make_tuple(r, Point({spaces[s].x, spaces[s].y, 0}), flipped));
		if( W < spaces[s].x + flippedRect.w ) W = spaces[s].x + flippedRect.w;
		if( H < spaces[s].y + flippedRect.h ) H = spaces[s].y + flippedRect.h;
						
		splitSpace(spaces, s, flippedRect);
		
	}						
	
	// Did we pack all rects?
	assert( packedRects.size() == rects.size() );
	
	
	// We are now tracking W,H live now!
	// calculateBounds(W, H, packedRects, rects);
		
	// return packing density: What % of output rectangle sheet (W*H) was packed with input rectangles
	return (area/(W*H));		
}


// Pack given rectangles (w,h) into a compact rectangle sheet (input/output: W, H) and their positions
// (index of rect, x, y, flip status).
// Returns packing density (ratio of packed area to overall output sheet area)
// Iterate to find best fit (number of iteration limited to NITER)
double packRectBest(std::vector<std::tuple<int, Point, bool>>& packedRects, double& W, double& H,const std::vector<Rect>& rects, const bool flipRect, const int NITER) {
	// We use simple algo: partioning output space into bins by splitting appropriately.
	// 1. Sort the input rectangles in decreasing order of their heights
	// 2. Initialize output rect space to accomodate at least maximum width rect
	// 3. For each input rect in sorted rects: 
	//		-Place it in smallest space bin possible
	//      -If found space is perfect match, remove this bin
	//		-If found space matches rect height/width, adjust its free space accoridingly
	//		-Otherwise split the space (new and adjust current one)
	
	// Create a sorted index on given rectangles in decreasing order of their heights
	// for same heights rectangles, prefer larger width rectangle
    std::vector<int> sortedRectIdx(rects.size());
  	std::iota(sortedRectIdx.begin(), sortedRectIdx.end(), 0);
  	std::sort(sortedRectIdx.begin(), sortedRectIdx.end(),
       		  [&rects](const int i1, const int i2) 
			  {return (rects[i1].h == rects[i2].h ? rects[i1].w > rects[i2].w : rects[i1].h > rects[i2].h);});
			  
	// Find maximum width among given rectangles
	double maxW = 0, totalH = 0, area = 0, maxH = 0, sqrSide =0;
	for(const auto& r: rects) {
		if( maxW < r.w ) maxW = r.w;
		if( maxH < r.h ) maxH = r.h;
		totalH += r.h;
		area += r.w * r.h;
	}

    sqrSide = (int)(sqrt(area*1.1)+0.5);
	
	// Initialize space to go for squarish fit: side of square containing total area (bit larger)	
	if( W == 0 ) 
		W = MAX(maxW, sqrSide);
	else 
		W = MAX(maxW, W);

	if( H == 0 )
		H = totalH;
	else
		H = MAX(totalH, H);	
		
	
	std::vector<std::tuple<int, Point, bool>> packedRectsIter;	
	double density = 0, maxDensity = 0;
	double bottomEdge = 0, maxBottomEdge=0, bottomW, iterW, iterH;
	int bottomRect = 0, i, sign = -1;
	bool flipped;
	Point loc;
	
	maxDensity = packRectFixedSize(packedRects, W, H, rects, flipRect, sortedRectIdx, area);
	packedRectsIter = packedRects; iterW = W; iterH = H;
	for(int k=1; k < NITER; k++) {
		// Adjust width and height and try the next fit
		
		// Find bottommost placed rect
		std::tie(i, loc, flipped) = packedRectsIter[0];
		maxBottomEdge = (!flipped ? loc.y + rects[i].h : loc.y + rects[i].w);
		bottomRect = 0;
		for(int j=1; j < packedRectsIter.size(); j++) {
			std::tie(i, loc, flipped) = packedRectsIter[j];
			bottomEdge = (!flipped ? loc.y + rects[i].h : loc.y + rects[i].w);
			if( bottomEdge > maxBottomEdge ) {
				maxBottomEdge = bottomEdge;
				bottomRect = j;
			}
		}
		
		// Adjust the width by bottommost rect width and try the next fit
		std::tie(i, loc, flipped) = packedRectsIter[bottomRect];
		bottomW = (!flipped ? rects[i].w : rects[i].h);
		
		iterW += sign * bottomW;
		if( iterW < maxW ) {
			iterW = maxW;
			sign = 1;
		}
		else {
			sign = -1;
		}
			
		iterH = totalH;
		density = packRectFixedSize(packedRectsIter, iterW, iterH, rects, flipRect, sortedRectIdx, area);
		
		// Keep track of best fit so far
		if(density > maxDensity) {
			maxDensity = density;
			packedRects = packedRectsIter;
			W = iterW; H = iterH;
		}		
				
	}
	
	return maxDensity;	
}


