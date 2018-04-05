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

#ifndef HEADERS_OFF_INPUT_FILE_H
#define HEADERS_OFF_INPUT_FILE_H

#include <Text_input_file.h>

#include <OFF_polygon.h>
#include <Point.h>

#include <vector>
#include <string>

#include <CGAL/Bbox_3.h>
#include <CGAL/IO/Color.h>

// this is the distance matrix used for matrix input
extern vector< vector<double> > g_distance_matrix;

// this is the adjacency list used for graph input
extern vector< vector< pair<int, double> > > g_adjacency_list;

// this is set to 0 for OFF input and 1 for MATRIX input.
extern int input_type_flag;

//This reads all of the input points.  The points can be in OFF format or matrix format.
namespace Headers
{
	template< typename Kernel_ >
	class OFF_input_file
		: public Text_input_file
	{
	public:

		OFF_input_file( char const *path_cstring_ );

		unsigned get_number_of_points() const;
		unsigned get_number_of_polygons() const;
		unsigned get_number_of_dimensions() const;	// dimension of the complex

		Point const &point_at( unsigned index_ ) const;
		OFF_polygon const &polygon_at( unsigned index_ ) const;

		// returns the upper and lower bound vectors at the specified dimension
		double get_lower_bound(unsigned dim) const;
		double get_upper_bound(unsigned dim) const;

	protected:

		void error( char const *what_ );

	protected:

		std::vector< Point > m_points;		// vector of points
		std::vector< OFF_polygon > m_polygons;	// vector of polygons
		int number_of_points;
		int number_of_polygons;
		

		// vector of points that form a lower and upper bound for the points in the cloud in each dimension.
		// For example, the lowest coordinate in the 5th dimension will be lower_bound[4] (assuming the first
		// dimension corresponds to the zeroth vector index)
		std::vector< double > lower_bound;
		std::vector< double > upper_bound;

		int num_dimensions;

		unsigned m_lines_read;
	};

	template< typename Kernel_ >
	inline
	OFF_input_file< Kernel_ >::OFF_input_file( char const *path_cstring_ )
		: Text_input_file( path_cstring_ ), m_lines_read( 0 )
	{
		using namespace std;
		using namespace CGAL;

		enum State
		{
			TYPE,
			OFF_HEADER,
			OFF_VERTEX,
			OFF_FACE,
			MATRIX_HEADER,
			MATRIX_READ,
			GRAPH_HEADER,
			GRAPH_READ,
			END
		};

		State state( TYPE );

		unsigned not_used;
		unsigned num_points_read = 0;

		static unsigned const BUFFER_SIZE( 1048576 );
		static char buffer[ BUFFER_SIZE ];
		static char token[ BUFFER_SIZE ];
			
		while ( true )
		{
			if ( !fgets( buffer, BUFFER_SIZE, m_p_file ) )
			{
				if ( state == END )
					break;

				error( "premature end of file" );
			}

			++m_lines_read;			

			char *p_line( buffer );

			// skip spaces
			while ( ( *p_line != '\0' ) && ( isspace( *p_line ) ) )
				++p_line;

			// if empty line
			if ( ( *p_line == '\n' ) || ( *p_line == '\0' ) )
				continue;

			// if comment
			if ( *p_line == '#' )
				continue;

			// number of characters read by sscanf
			unsigned read( 0 );

			switch ( state )
			{
			case TYPE:
				{
					sscanf( p_line, "%s%n", token, &read );

					// the first line of the file should either be "OFF" or "MATRIX"
					if ( strcmp( token, "OFF" ) && strcmp( token, "MATRIX") && strcmp( token, "GRAPH" ))
						error( "OFF, MATRIX, or GRAPH header expected" );

					p_line += read;
					
					if( strcmp(token, "OFF") == 0)
					{
						input_type_flag = 0;
						state = OFF_HEADER;
					}
					else if( strcmp(token, "MATRIX") == 0)
					{
						input_type_flag = 1;
						state = MATRIX_HEADER;
					}
					else if( strcmp(token, "GRAPH") == 0)
					{
						input_type_flag = 2;
						state = GRAPH_HEADER;
					}
					else
						error("Header should either be OFF or MATRIX");
						
				}
				break;
			case OFF_HEADER:
				{
					if ( sscanf( p_line, "%d %d %d%n", &number_of_points,
						&number_of_polygons, &not_used, &read ) <= 0 )
					{
						error( "numbers of vertices, faces and edges expected" );
					}

					p_line += read;

					if ( number_of_points > 0 )
					{
						state = OFF_VERTEX;
						m_points.reserve( number_of_points );
					}
					else
						state = END;

				}
				break;
			case OFF_VERTEX:
				{
					// counts the number of dimensions in the file
					if(m_points.size() == 0 )
					{
						num_dimensions = 0;
						double temp_value;
						char *temp_line = p_line;
						while(sscanf(temp_line, " %lf%n", &temp_value, &read) > 0)
						{
							temp_line += read;
							num_dimensions++;
						}
						// allocate space for our boundary vectors
						lower_bound.resize(num_dimensions);
						upper_bound.resize(num_dimensions);
					}

					// vertex coordinates
					Point P(num_dimensions);
					double coord;
					for(int i=0; i<num_dimensions; i++)
					{
						if ( sscanf( p_line, "%lf%n", &coord, &read ) <= 0 )
							error( "vertex coordinates expected" );
						p_line += read;
						P.set_coord(i, coord);
						// If we haven't read any points yet, the upper and lower bounds for all of our
						// points are set to the coordinates of the first point we read
						if ( m_points.size() == 0 )
						{
							lower_bound[i] = coord;
							upper_bound[i] = coord;
						}
						// Check to see if we need to update our upper/lower bounds
						else
						{
							if(lower_bound[i] > coord)
								lower_bound[i] = coord;
							else if(upper_bound[i] < coord)
								upper_bound[i] = coord;
						}
					}					

					m_points.push_back( P );
					// if done with points
					if ( m_points.size() == number_of_points )
					{
						if ( number_of_polygons > 0 )
						{
							state = OFF_FACE;
							m_polygons.reserve( number_of_polygons );
						}
						else
							state = END;
					}
				}
				break;
			case OFF_FACE:
				{
					unsigned size;
					if ( sscanf( p_line, "%d%n", &size, &read ) <= 0 )
						error( "polygon size expected" );

					p_line += read;

					OFF_polygon polygon;
					polygon.reserve( size );

					// vertex indices
					for ( unsigned i( 0 ); i < size; ++i )
					{
						unsigned index;
						if ( sscanf( p_line, "%d%n", &index, &read ) <= 0 )
							error( "vertex index expected" );

						if ( index > number_of_points - 1 )
							error( "invalid vertex index" );

						p_line += read;
						polygon.add_point_index( index );
					}
					
					// color & alpha
					double r, g, b;
					if ( sscanf( p_line, "%lf %lf %lf%n", &r, &g, &b, &read ) > 0 )
					{
						p_line += read;
						if ( ( r <= 1 ) && ( g <= 1 ) && ( b <= 1 ) )
						{
							r *= 255;
							g *= 255;
							b *= 255;
						}
						polygon.set_color( Color(
							static_cast< unsigned char >( r ),
							static_cast< unsigned char >( g ),
							static_cast< unsigned char >( b ) ) );

						double alpha;
						if ( sscanf( p_line, "%lf%n", &alpha, &read ) > 0 )
						{
							p_line += read;
							// process alpha
						}
					}

					m_polygons.push_back( polygon );

					// if done with polygons
					if ( m_polygons.size() == number_of_polygons )
						state = END;
				}
				break;
			case MATRIX_HEADER:
				{
				if ( sscanf( p_line, "%d%n", &number_of_points,
						&read ) <= 0 )
					{
						error( "numbers of vertices expected" );
					}

					p_line += read;

					if ( number_of_points > 0 )
					{
						g_distance_matrix.resize( number_of_points );
						for(int i=0; i<number_of_points; i++)
						{
							g_distance_matrix[i].resize(number_of_points);
						}
						state = MATRIX_READ;
					}
					else
						state = END;

				}
				break;
			case MATRIX_READ:
				{
					if(num_points_read == 0)
					{
						g_distance_matrix[num_points_read][num_points_read] = 0;
						num_points_read++;
					}

					double dist;
					// loops through all points of the matrix
					for(int j=0; j<num_points_read; j++)
					{
						if ( sscanf( p_line, "%lf%n", &dist, &read ) <= 0 )
							error( "distance expected" );
						p_line += read;
						g_distance_matrix[num_points_read][j] = dist;
						g_distance_matrix[j][num_points_read] = dist;
					}

					g_distance_matrix[num_points_read][num_points_read] = 0;
				
					num_points_read++;
					if(num_points_read == number_of_points)
						state = END;
				}
				break;
			case GRAPH_HEADER:
				{
				// reads in the number of input points
				if ( sscanf( p_line, "%d%n", &number_of_points,
						&read ) <= 0 )
					{
						error( "numbers of vertices expected" );
					}

					p_line += read;

					// resizes our adjacency list and goes to the next phase
					if ( number_of_points > 0 )
					{
						g_adjacency_list.resize( number_of_points );
						state = GRAPH_READ;
					}
					else
						state = END;

				}
				break;
			// reads in the graph data from the file
			case GRAPH_READ:
				{
					int current_node, neighbor_node;
					double dist;

					if ( sscanf( p_line, "%d%n", &current_node, &read ) <= 0)
						error( "point index expected");
					p_line += read;

					// gets the edges from the graph
					while(sscanf(p_line, "%d%n", &neighbor_node, &read) > 0)
					{
						p_line += read;
						if(sscanf( p_line, "%lf%n", &dist, &read ) <= 0)
							error("distance expected");
						p_line += read;
						g_adjacency_list[current_node].push_back(make_pair(neighbor_node, dist));
						
					}
				
					num_points_read++;
					if(num_points_read == number_of_points)
						state = END;
				}
				break;
			case END:
			default:
				{
				}
			}

			// skip spaces
			while ( ( *p_line != '\0' ) && ( isspace( *p_line ) ) )
				++p_line;

			if ( ( *p_line != '\n' ) && ( *p_line != '\0' ) )
				error( "extra characters" );
		}		
	}
	
	template< typename Kernel_ >
	inline unsigned
	OFF_input_file< Kernel_ >::get_number_of_points() const
	{
		return number_of_points;
	}

	template< typename Kernel_ >
	inline unsigned
	OFF_input_file< Kernel_ >::get_number_of_dimensions() const
	{
		return num_dimensions;
	}
	
	template< typename Kernel_ >
	inline unsigned
	OFF_input_file< Kernel_ >::get_number_of_polygons() const
	{
		return number_of_polygons;
	}
	
	template< typename Kernel_ >
	inline Point const &
	OFF_input_file< Kernel_ >::point_at( unsigned index_ ) const
	{
		return m_points.at( index_ );
	}

	template< typename Kernel_ >
	inline OFF_polygon const &
	OFF_input_file< Kernel_ >::polygon_at( unsigned index_ ) const
	{
		return m_polygons.at( index_ );
	}

	template< typename Kernel_ >
	inline double 
	OFF_input_file< Kernel_ >::get_upper_bound(unsigned dim) const
	{
		return upper_bound[dim];
	}

	template< typename Kernel_ >
	inline double 
	OFF_input_file< Kernel_ >::get_lower_bound(unsigned dim) const
	{
		return lower_bound[dim];
	}
	
	template< typename Kernel_ >
	inline void 
	OFF_input_file< Kernel_ >::error( char const *what_ )
	{
		static char where[ 16 ];
		sprintf( where, "%d", m_lines_read );

		throw Exception( ( path().string() + ": " + what_
			+ " at line " + where ).c_str() );
	}
}
#endif // COMMON_OFF_INPUT_FILE_H
