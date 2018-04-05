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

#ifndef HEADERS_GRID_H
#define HEADERS_GRID_H

#include <list>
#include <Point.h>

#include <boost/math/special_functions/round.hpp>

namespace Headers
{	
	template< typename Kernel_, typename Data_ >
	class Grid
	{

	public:

		class Cell
		{
			// this means the grid class can access private/protected members of the Cell class
			friend class Grid;

		public:
			typedef std::list< Data_ > Datas;
			typedef typename Datas::const_iterator Data_const_iterator;
			typedef typename Datas::iterator Data_iterator;

		public:		
			int get_coord(int dimension) const;
			
			unsigned number_of_datas() const;
			Data_const_iterator datas_begin() const;
			Data_const_iterator datas_end() const;
			Data_iterator datas_begin();
			Data_iterator datas_end();
			void add( Data_ const &data_ );

			// does it have a neighbor at the index corresponding to the direction set in the
			// enum above?
			bool has_neighbor_at( vector<int> direction ) const;
			Cell const &neighbor_at( vector<int> direction ) const;
			Cell &neighbor_at( vector<int> direction );

			double get_squared_distance_to( Point const &point_ ) const;

		private:

			// These are the coordinates of the lower corner of the cell in each dimension.
			vector<int> cell_coords;

			// This is a vector of all of the vertices in this particular cell
			Datas m_datas;

			// parent grid
			Grid *m_p_grid;
		};	

	public:

		Grid( vector<double> lower_bound, vector<double> upper_bound, vector<int> resolution, int num_dimensions );
		~Grid();

		double const get_lower_bound(int dim) const;
		double const get_upper_bound(int dim) const;

		int get_resolution(int dim) const;

		// gets the cell that the specified coordinates belong to
		Cell const &get_cell_at( vector<int> coordinates ) const;
		Cell &get_cell_at( vector<int> coordinates );

		// Gets the cell that the point is in.
		Cell const &get_cell_by_point( Point const &point_ ) const;
		Cell &get_cell_by_point( Point const &point_ );

		Cell &add_at( Data_ const &data_, Point const &location_ );

		int get_num_dimensions() const;	// gets the number of dimensions of the grid

	private:
		
		// array of cells in the grid
		Cell *m_p_cells;

		vector<double> lower_bound;
		vector<double> upper_bound;

		vector<int> resolution;
		int num_dimensions;
	};

	template< typename Kernel_, typename Data_ >
	inline int
	Grid< Kernel_, Data_ >::Cell::get_coord(int dimension) const
	{
		return cell_coords.at(dimension);
	}

	// returns the size of the datas list
	template< typename Kernel_, typename Data_ >
	inline unsigned
	Grid< Kernel_, Data_ >::Cell::number_of_datas() const
	{
		return m_datas.size();
	}
	
	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell::Data_const_iterator
	Grid< Kernel_, Data_ >::Cell::datas_begin() const
	{
		return m_datas.begin();
	}

	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell::Data_const_iterator
	Grid< Kernel_, Data_ >::Cell::datas_end() const
	{
		return m_datas.end();
	}

	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell::Data_iterator
	Grid< Kernel_, Data_ >::Cell::datas_begin()
	{
		return m_datas.begin();
	}

	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell::Data_iterator
	Grid< Kernel_, Data_ >::Cell::datas_end()
	{
		return m_datas.end();
	}
	
	//adds the data (a vertex) to the data vector of the cell
	template< typename Kernel_, typename Data_ >
	inline void
	Grid< Kernel_, Data_ >::Cell::add( Data_ const &data_ )
	{
		m_datas.push_back( data_ );		
	}
		
	// checks to see if the current cell has a neighbor at a specific location.  The direction vector
	// specifies which direction we are checking for the neighbor.
	template< typename Kernel_, typename Data_ >
	inline bool
	Grid< Kernel_, Data_ >::Cell::has_neighbor_at( vector<int> direction ) const
	{
		int dimensions = m_p_grid->get_num_dimensions();
		int neighbor_coordinate;
		// checks to see if the neighbor coordinate is in the grid in each dimension
		for(int i=0; i<dimensions; i++)
		{
			neighbor_coordinate = cell_coords.at(i) + direction.at(i);
			if(neighbor_coordinate < 0 || neighbor_coordinate >= m_p_grid->get_resolution(i))
				return false;
		}
		return true;
	}

	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell const &
	Grid< Kernel_, Data_ >::Cell::neighbor_at( vector<int> direction ) const
	{
		typename Grid< Kernel_, Data_ >::Cell const *p_self( this );
		return p_self->neighbor_at(direction);

		return *this;
	}

	// gets the neighboring cell specified by the neighbor vector
	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell &
	Grid< Kernel_, Data_ >::Cell::neighbor_at( vector<int> direction )
	{
		int dimensions = m_p_grid->get_num_dimensions();
		vector<int> neighbor_coordinates(dimensions);
		// gets the cell coordinates of the desired neighbor
		for(int i=0; i<dimensions; i++)
			neighbor_coordinates.at(i) = cell_coords.at(i) + direction.at(i);
		return m_p_grid->get_cell_at(neighbor_coordinates);
	}

	// gets the squared Euclidian distance from the cell to a point.
	template< typename Kernel_, typename Data_ >
	inline double
	Grid< Kernel_, Data_ >::Cell::get_squared_distance_to(Point const &point_ ) const
	{
		Cell const &particle_cell( m_p_grid->get_cell_by_point( point_ ) );
		int dimensions = m_p_grid->get_num_dimensions();

		vector<double> distances(num_dimensions);
		double total_distance = 0;
		double width;			
		double current_distance;
		double lower_bound, upper_bound;

		// Calculates the squared distance in each dimension between the cell and the point.
		for(int i=0; i<num_dimensions; i++)
		{
			upper_bound = m_p_grid->get_upper_bound(i);
			lower_bound = m_p_grid->get_lower_bound(i);
			width = (upper_bound - lower_bound)/m_p_grid->get_resolution(i);
			if(cell_coords.at(i) > particle_cell.get_coord(i))
			{
				current_distance = (lower_bound + cell_coords.at(i)*width) - point_.get_coord(i);
				total_distance += current_distance*current_distance;
			}
			else if(cell_coords.at(i) < particle_cell.get_coord(i))
			{
				current_distance = point_.get_coord(i) - (lower_bound + cell_coords.at(i)*width); 
				total_distance += current_distance*current_distance;
			}
		}
	}

	template< typename Kernel_, typename Data_ >
	inline
	Grid< Kernel_, Data_ >::Grid( vector<double> lower, vector<double> upper, vector<int> res, int dims )
		: lower_bound(lower), upper_bound(upper), resolution(res), num_dimensions(dims)
	{
		int total_resolution = 1;
		int i;
		vector<int> coordinates(num_dimensions);
		int dimension_index;
		
		for(i=0; i<num_dimensions; i++)
		{
			total_resolution *= resolution[i];
			coordinates[i] = 0;
		}

		// we store the entire grid in a 1 dimensional vector.  It will have total_resoluion elements.
		m_p_cells = new Cell [total_resolution];
		
		// initializes each cell in the grid.  Sets their coordinates and their parent grid.
		m_p_cells[0].cell_coords = coordinates;
		m_p_cells[0].m_p_grid = this;
		for(i=1; i<total_resolution; i++)
		{
			dimension_index = 0;
			coordinates[dimension_index]++;
			while(coordinates[dimension_index] >= resolution[dimension_index])
			{
				coordinates[dimension_index] = 0;
				dimension_index++;
				coordinates[dimension_index]++;
			}

			m_p_cells[i].cell_coords = coordinates;
			m_p_cells[i].m_p_grid = this;
		}
		
	}	
	
	template< typename Kernel_, typename Data_ >
	inline
	Grid< Kernel_, Data_ >::~Grid()
	{
		delete [] m_p_cells;
	}	

	template< typename Kernel_, typename Data_ >
	inline double const
	Grid< Kernel_, Data_ >::get_upper_bound(int dim) const
	{
		return upper_bound[dim];
	}

	template< typename Kernel_, typename Data_ >
	inline double const
	Grid< Kernel_, Data_ >::get_lower_bound(int dim) const
	{
		return lower_bound[dim];
	}

	template< typename Kernel_, typename Data_ >
	inline int
	Grid< Kernel_, Data_ >::get_resolution(int dim) const
	{
		return resolution[dim];
	}
	
	// returns the cell at the specified grid coordinates
	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell const &
	Grid< Kernel_, Data_ >::get_cell_at( vector<int> coordinates ) const
	{
		int index = 0;
		int factor = 1;
		for(int i=0; i<num_dimensions; i++)
		{
			index += coordinates[i]*factor;
			factor *= resolution[i];
		}
		return m_p_cells[index];
	}

	// returns the cell at the specified grid coordinates
	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell &
	Grid< Kernel_, Data_ >::get_cell_at( vector<int> coordinates )
	{
		int index = 0;
		int factor = 1;
		for(int i=0; i<num_dimensions; i++)
		{
			index += coordinates[i]*factor;
			factor *= resolution[i];
		}
		return m_p_cells[index];
	}

	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell const &
	Grid< Kernel_, Data_ >::get_cell_by_point( Point const &point_ ) const
	{
		Grid< Kernel_, Data_ > const *p_const_this( this );
		return p_const_this->get_cell_by_point( point_ );
	}

	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell &
	Grid< Kernel_, Data_ >::get_cell_by_point( Point const &point_ )
	{	
		using namespace boost;
		using namespace math;
		
		double width, cell_width, relative_coord;
		vector<int> new_cell_coords(num_dimensions);
		for(int i=0; i<num_dimensions; i++)
		{
			if(point_.get_coord(i) == upper_bound[i])
				new_cell_coords[i] = resolution[i] - 1;
			else
			{
				width = (upper_bound[i] - lower_bound[i])/resolution[i];
				relative_coord = point_.get_coord(i) - lower_bound[i];
				new_cell_coords[i] = iround<double>(floor(relative_coord / width));
			}
		}
		return get_cell_at(new_cell_coords);
	}

	template< typename Kernel_, typename Data_ >
	inline typename Grid< Kernel_, Data_ >::Cell &
	Grid< Kernel_, Data_ >::add_at( Data_ const &data_, Point const &location_ )
	{
		Cell &cell( get_cell_by_point( location_ ) );
		cell.add( data_ );		

		return cell;
	}

	template< typename Kernel_, typename Data_ >
	inline int
	Grid< Kernel_, Data_ >::get_num_dimensions() const
	{
		return num_dimensions;
	}
}

#endif // COMMON_GRID_H
