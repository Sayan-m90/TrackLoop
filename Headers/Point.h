///////////////////////////////////////////////////////////////////////////////
//
// THIS SOFTWARE IS PROVIDED "AS-IS". THERE IS NO WARRANTY OF ANY KIND.
// NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE
// FOR ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY.
//
// Copyright (c) 2010 Jyamiti Research Group.
// CS&E Department of the Ohio State University, Columbus, OH.
// All rights reserved.
//
// Author: Oleksiy Busaryev.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef POINT_H
#define POINT_H

#include <vector>
#include <string>
//template< typename Kernel_ >
using namespace std;
class Point
{

	public:
	
		double get_coord(int index) const;	// gets the coordinate of the point in the indexth dimension.
		int get_dim() const;			// returns the dimension of the point.
		void set_coord(int dim, double coord);	// sets the dimth dimensoin of the point to coord
		Point(int dim);				// constructor
		Point();				// default constructor
		double get_squared_distance_to(Point const p) const;	// Gets the squared distance between this and point p.
		
	public:

		vector<double> coordinates;		// stores the coordinates of the point
		int dimensions;				// dimension of the point
};	

Point::Point(int dim)
{
	dimensions = dim;
	coordinates.resize(dim);
}

Point::Point()
{
}

double Point::get_coord(int index) const
{
	return coordinates[index];
}

int Point::get_dim() const
{
	return dimensions;
}

void Point::set_coord(int dim, double coord)
{
	coordinates[dim] = coord;
}

double Point::get_squared_distance_to(const Point p) const
{
	double dist;
	double total_dist = 0;
	for(int i=0; i<dimensions; i++)
	{
		dist = coordinates[i] - p.get_coord(i);
		total_dist += dist*dist;
	}
	return total_dist;
}
#endif
