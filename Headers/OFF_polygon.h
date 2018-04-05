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

#ifndef HEADERS_OFF_POLYGON_H
#define HEADERS_OFF_POLYGON_H

#include <Colored.h>

#include <vector>

namespace Headers
{
	class OFF_polygon
		: public Colored
	{

	public:

		// number of vertices in the polygon
		unsigned size() const;
		
		unsigned point_index_at( unsigned i_ ) const;
		void reserve( unsigned number_of_points_ );
		void add_point_index( unsigned point_index_ );

	protected:

		// vector of indices of the vertices.  These indices correspond to the vector
		// indices from the OFF_input_file
		std::vector< unsigned > m_vertex_indices;
	};

	inline unsigned
	OFF_polygon::size() const
	{
		return m_vertex_indices.size();
	}
	
	inline unsigned
	OFF_polygon::point_index_at( unsigned i_ ) const
	{
		return m_vertex_indices.at( i_ );
	}

	inline void
	OFF_polygon::reserve( unsigned number_of_points_ )
	{
		m_vertex_indices.reserve( number_of_points_ );
	}

	inline void
	OFF_polygon::add_point_index( unsigned point_index_ )
	{
		m_vertex_indices.push_back( point_index_ );
	}
}

#endif // COMMON_OFF_POLYGON_H
