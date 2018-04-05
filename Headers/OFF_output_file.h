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

#ifndef HEADERS_OFF_OUTPUT_FILE_H
#define HEADERS_OFF_OUTPUT_FILE_H

#include <Text_output_file.h>

#include <OFF_polygon.h>

#include <vector>
#include <iostream>

namespace Headers
{
	template< typename Kernel_ >
	class OFF_output_file
		: public Text_output_file
	{
	public:

		OFF_output_file( char const *path_cstring_ );
		~OFF_output_file();

		void add_point( Point const &point_ );
		void add_polygon( OFF_polygon const &polygon_ );

	protected:

		std::vector< Point > m_points;
		std::vector< OFF_polygon > m_polygons;
	};

	template< typename Kernel_ >
	inline
	OFF_output_file< Kernel_ >::OFF_output_file( char const *path_cstring_ )
		: Text_output_file( path_cstring_ )
	{		
	}

	template< typename Kernel_ >
	inline
	OFF_output_file< Kernel_ >::~OFF_output_file()
	{		
		fprintf( m_p_file, "OFF\n" );
		fprintf( m_p_file, "%d %d 0\n", m_points.size(), m_polygons.size() );

		// prints the coordinates of each point to the output file
		for ( unsigned i( 0 ); i < m_points.size(); ++i )
		{
			Point const &point( m_points.at( i ) );
			
			for(int j = 0; j < point.get_dim(); j++)
				fprintf( m_p_file, "%lf ", point.get_coord(j));
			fprintf(m_p_file, "\n");
		}

		// prints the size of the polygons and the index of them that are in the basis
		for ( unsigned i( 0 ); i < m_polygons.size(); ++i )
		{
			OFF_polygon const &polygon( m_polygons.at( i ) );
			fprintf( m_p_file, "%d", polygon.size() );

			for ( unsigned j( 0 ); j < polygon.size(); ++j )
				fprintf( m_p_file, " %d", polygon.point_index_at( j ) );

			CGAL::Color const &color( polygon.color() );
			fprintf( m_p_file, " %d %d %d", color.red(),
				color.green(), color.blue() );

			fprintf( m_p_file, "\n" );
		}
	}

	template< typename Kernel_ >
	inline void
	OFF_output_file< Kernel_ >::add_point( Point const &point_ )
	{
		m_points.push_back( point_ );
	}

	template< typename Kernel_ >
	inline void
	OFF_output_file< Kernel_ >::add_polygon( OFF_polygon const &polygon_ )
	{
		m_polygons.push_back( polygon_ );
	}
}

#endif // COMMON_OFF_OUTPUT_FILE_H
