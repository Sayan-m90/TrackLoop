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

#ifndef HOM_BAS_COMPLEX_H
#define HOM_BAS_COMPLEX_H

#include <vector>
#include <algorithm>
#include <limits>

#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/pending/mutable_queue.hpp>
#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <boost/checked_delete.hpp>
#include <boost/math/special_functions/round.hpp>

#include <CGAL/Homogeneous_d.h>
#include <CGAL/predicates_d.h>

#include <Grid.h>

#include <Indexed.h>
#include <Flagged.h>
#include <Normed.h>

// this is the distance matrix used for matrix input
extern vector< vector<double> > g_distance_matrix;

// this is the adjacency list used for graph input
extern vector< vector< pair<int, double> > > g_adjacency_list;

// this is set to 0 for OFF input and 1 for MATRIX input.
extern int input_type_flag;

namespace trackLoop
{
	#undef max
	#undef min

	template< typename Kernel_ >
	class Edge;

	///////////////////////////////////////////////////////////////////////////
	//
	// Represents a vertex in a simplicial complex.
	//
	///////////////////////////////////////////////////////////////////////////

	template< typename Kernel_ >
	class Vertex
		: public Headers::Indexed, public Headers::Flagged
	{

	public:

		enum
		{
			IS_IN_QUEUE = 1,
			IS_IN_TREE = 1 << 1,
			IS_IN_SAMPLE = 1 << 2
		};
				
		typedef std::vector< Edge< Kernel_ > * > Coboundary;

	public:

		Vertex( Point const &location_ );
		Vertex();

		Point const &location() const;

		Coboundary const &coboundary() const;
		Coboundary &coboundary();

		Vertex< Kernel_ > const &coneighbor( Edge< Kernel_ > const &edge_ ) const;
		Vertex< Kernel_ > &coneighbor( Edge< Kernel_ > &edge_ );
		
		double distance_to_root() const;
		void set_distance_to_root( double distance_to_root_ );
		
		boost::dynamic_bitset<> const &contracted_path_to_root() const;
		boost::dynamic_bitset<> &contracted_path_to_root();
				
		bool has_parent() const;
		Vertex const &parent() const;
		Vertex &parent();
		Edge< Kernel_ > const &edge_to_parent() const;
		Edge< Kernel_ > &edge_to_parent();
		void set_edge_to_parent( Edge< Kernel_ > &edge_to_parent_ );
		void clear_edge_to_parent();

		Vertex const &image() const;
		Vertex &image();
		void set_image( Vertex &image_ );

	public:

		Point m_location;
		// ID for vertices
		int unique_id;
		// Edges having this vertex as a face
		Coboundary m_coboundary;

		// Distance to the root in some shortest path tree
		double m_distance_to_root;

		// Path to the root in the contracted complex
		boost::dynamic_bitset<> m_contracted_path_to_root;

		// Edge connecting this vertex to its parent in some shortest path tree
		Edge< Kernel_ > *m_p_edge_to_parent;

		// Image in the contracted complex
		Vertex *m_p_image;
	};

	///////////////////////////////////////////////////////////////////////////
	//
	// Represents an edge in a simplicial complex.
	//
	///////////////////////////////////////////////////////////////////////////

	template< typename Kernel_ >
	class Edge
		: public Headers::Indexed, public Headers::Flagged
	{

	public:

		enum
		{
			IS_IN_TREE = 1,
			IS_IN_LOOP = 1 << 2
		};

	public:

		Edge( Vertex< Kernel_ > &a_, Vertex< Kernel_ > &b_, double length);

		Vertex< Kernel_ > const &a() const;
		Vertex< Kernel_ > &a();

		Vertex< Kernel_ > const &b() const;
		Vertex< Kernel_ > &b();

		double length() const;

		double canonical_loop_length() const;
		void set_canonical_loop_length( double canonical_loop_length_ );

		std::vector< Edge< Kernel_ > * > const &image() const;
		std::vector< Edge< Kernel_ > * > &image();

	public:

		Vertex< Kernel_ > &m_a, &m_b;
		
		double m_length;

		// Length of the canonical loop in some shortest path tree;
		// canonical loop is formed by the edge and two chains connecting
		// the endpoints of the edge to the root of the shortest path tree
		double m_canonical_loop_length;

		// Image in the contracted complex; may contain multiple edges 
		// (e.g. when contracting a triangle abc, ab may be mapped to {ac, bc}
		std::vector< Edge< Kernel_ > * > m_image;
	};

	///////////////////////////////////////////////////////////////////////////
	//
	// Represents a triangle in a simplicial complex.
	//
	///////////////////////////////////////////////////////////////////////////

	template< typename Kernel_ >
	class Triangle
	{

	public:

		Triangle( Edge< Kernel_ > &ab_, Edge< Kernel_ > &ac_,
			Edge< Kernel_ > &bc_ );

		Edge< Kernel_ > const &ab() const;
		Edge< Kernel_ > &ab();

		Edge< Kernel_ > const &ac() const;
		Edge< Kernel_ > &ac();

		Edge< Kernel_ > const &bc() const;
		Edge< Kernel_ > &bc();

	public:

		Edge< Kernel_ > &m_ab, &m_ac, &m_bc;
	};

	///////////////////////////////////////////////////////////////////////////
	//
	// Defines a canonical loop representation in the contracted complex.
	// Canonical loop is formed by an edge and two chains connecting
	// the endpoints of this edge to the root of a shortest path tree rooted
	// at some vertex. In the contracted complex, canonical loop can be
	// represented in terms of basis loops (essential homology cycles),
	// which have 1-to-1 correspondence with edges of the contracted complex.
	//
	///////////////////////////////////////////////////////////////////////////

	template< typename Kernel_ >
	class Canonical_loop
		: public Headers::Normed
	{

	public:

		Canonical_loop( Vertex< Kernel_ > &vertex_, Edge< Kernel_ > &edge_,
			boost::dynamic_bitset<> const &bits_ );

		Vertex< Kernel_ > const &vertex() const;
		Edge< Kernel_ > const &edge() const;

		Vertex< Kernel_ > &vertex();
		Edge< Kernel_ > &edge();

		boost::dynamic_bitset<> const &bits() const;
		boost::dynamic_bitset<> &bits();

	public:

		Vertex< Kernel_ > &m_vertex;
		Edge< Kernel_ > &m_edge;

		// 1's in this bitset correspond to essential homology cycles
		// (or to the edges in the contracted complex)
		boost::dynamic_bitset<> m_bits;
	};

	///////////////////////////////////////////////////////////////////////////
	//
	// Represents a loop in the shortest homology basis
	//
	///////////////////////////////////////////////////////////////////////////

	template< typename Kernel_ >
	class Basis_loop
		: public Headers::Normed
	{
		typedef std::vector< Edge< Kernel_ > * > Container;

	public:

		typedef typename Container::const_iterator Const_iterator;
		typedef typename Container::iterator Iterator;

	public:

		unsigned size() const;

		Const_iterator begin() const;
		Const_iterator end() const;

		Iterator begin();
		Iterator end();

		void push_back( Edge< Kernel_ > *p_edge_ );

	public:

		Container m_edges;
	};

	///////////////////////////////////////////////////////////////////////////
	//
	// Represents a simplicial complex.
	//
	///////////////////////////////////////////////////////////////////////////

	template< typename Kernel_ >
	class Complex
	{

	public:

		Complex( int dimensions, bool verbose_ );
		Complex( bool verbose_ );
		~Complex();
		
		unsigned number_of_vertices() const;
		Vertex< Kernel_ > &vertex_at( unsigned i_ );
		void insert_vertex( Vertex< Kernel_ > *p_vertex_);

		unsigned number_of_edges() const;
		Edge< Kernel_ > &edge_at( unsigned i_ );
		void create_edge( Vertex< Kernel_ > &a_, Vertex< Kernel_ > &b_, double length);

		unsigned number_of_triangles() const;
		Triangle< Kernel_ > &triangle_at( unsigned i_ );
		void create_triangle( Vertex< Kernel_ > &a_,
			Vertex< Kernel_ > &b_, Vertex< Kernel_ > &c_ );

		// Performs the Rips expansion (only edges)
		void expand( double alpha_ );
		void expand_matrix( double alpha_ );
		void expand_graph( double alpha_ );

		// Computes the contracted complex
		void contract();
		void contract_if_needed( Edge< Kernel_ > &edge_ );
		void contract_if_needed( Triangle< Kernel_ > &triangle_ );

		void sample( double coefficient_ );

		// returns the shortest path to each node from the node "src"
		vector<double> shortest_path_graph(int src, double alpha);

		unsigned basis_rank() const;
		Basis_loop< Kernel_ > const &basis_loop_at( unsigned i_ ) const;
		Basis_loop< Kernel_ > &basis_loop_at( unsigned i_ );
		void compute_basis();

	public:		

		// If the edge survives in the contracted complex, e2b returns
		// the index of corresponding homology cycle (from 0 to 2g-1);
		// otherwise, e2b returns -1
		unsigned e2b( unsigned edge_index_ ) const;

		void compute_canonical_loops_for( Vertex< Kernel_ > &vertex_ );
		void compute_canonical_loop_lengths();
		void compute_shortest_path_tree_for( Vertex< Kernel_ > &vertex_ );
		void clear_shortest_path_tree();
		
	public:

		bool m_verbose;

		std::vector< Vertex< Kernel_ > * > m_vertices;
		std::vector< Edge< Kernel_ > * > m_edges;
		std::vector< Triangle< Kernel_ > * > m_triangles;
		
		typedef std::pair< Vertex< Kernel_ > *, Vertex< Kernel_ > * > VV;
		typedef boost::unordered_map< VV, Edge< Kernel_ > * > VV2E;
		VV2E m_vv2e;

		unsigned m_tree_size;
		
		std::vector< unsigned > m_e2b;
		std::vector< Canonical_loop< Kernel_ > * > m_canonical_loops;
		std::vector< Basis_loop< Kernel_ > * > m_basis_loops;

		boost::timer m_timer;
		boost::progress_display *m_p_progress;

		std::vector< double > lower_bound;
		std::vector< double > upper_bound;

		// number of dimensions in the complex
		int num_dimensions;
		
		bool m_expanded;
	};

	template< typename Kernel_ >
	inline
	Vertex< Kernel_ >::Vertex( Point const &location_ )
		: m_location( location_ ), m_distance_to_root( 0 ),
		m_p_edge_to_parent( 0 )
	{
		// Initially, the image of each vertex is the vertex itself
		m_p_image = this;
	}

	template< typename Kernel_ >
	inline
	Vertex< Kernel_ >::Vertex()
	{
		m_distance_to_root = 0;
		m_p_edge_to_parent = 0;
		// Initially, the image of each vertex is the vertex itself
		m_p_image = this;
	}

	template< typename Kernel_ >
	inline Point const &
	Vertex< Kernel_ >::location() const
	{
		return m_location;
	}

	template< typename Kernel_ >
	inline typename Vertex< Kernel_ >::Coboundary const &
	Vertex< Kernel_ >::coboundary() const
	{
		return m_coboundary;
	}

	template< typename Kernel_ >
	inline typename Vertex< Kernel_ >::Coboundary &
	Vertex< Kernel_ >::coboundary()
	{
		return m_coboundary;
	}
	
	template< typename Kernel_ >
	inline Vertex< Kernel_ > const &
	Vertex< Kernel_ >::coneighbor( Edge< Kernel_ > const &edge_ ) const
	{
		return &edge_.a() == this ? edge_.b() : edge_.a();
	}
	
	template< typename Kernel_ >
	inline Vertex< Kernel_ > &
	Vertex< Kernel_ >::coneighbor( Edge< Kernel_ > &edge_ )
	{
		return &edge_.a() == this ? edge_.b() : edge_.a();
	}

	template< typename Kernel_ >
	inline double
	Vertex< Kernel_ >::distance_to_root() const
	{
		return m_distance_to_root;
	}

	template< typename Kernel_ >
	inline void
	Vertex< Kernel_ >::set_distance_to_root( double distance_to_root_ )
	{
		m_distance_to_root = distance_to_root_;
	}

	template< typename Kernel_ >
	inline boost::dynamic_bitset<> const &
	Vertex< Kernel_ >::contracted_path_to_root() const
	{
		return m_contracted_path_to_root;
	}

	template< typename Kernel_ >
	inline boost::dynamic_bitset<> &
	Vertex< Kernel_ >::contracted_path_to_root()
	{
		return m_contracted_path_to_root;
	}

	template< typename Kernel_ >
	inline bool
	Vertex< Kernel_ >::has_parent() const
	{
		return m_p_edge_to_parent != 0;
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > const & 
	Vertex< Kernel_ >::parent() const
	{
		if ( &edge_to_parent().a() == this )
			return edge_to_parent().b();
		else
			return edge_to_parent().a();
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > & 
	Vertex< Kernel_ >::parent()
	{
		if ( &edge_to_parent().a() == this )
			return edge_to_parent().b();
		else
			return edge_to_parent().a();
	}

	template< typename Kernel_ >
	inline Edge< Kernel_ > const &
	Vertex< Kernel_ >::edge_to_parent() const
	{
		return *m_p_edge_to_parent;
	}

	template< typename Kernel_ >
	inline Edge< Kernel_ > &
	Vertex< Kernel_ >::edge_to_parent()
	{
		return *m_p_edge_to_parent;
	}

	template< typename Kernel_ >
	inline void
	Vertex< Kernel_ >::set_edge_to_parent( Edge< Kernel_ > &edge_to_parent_ )
	{
		m_p_edge_to_parent = &edge_to_parent_;
	}

	template< typename Kernel_ >
	inline void
	Vertex< Kernel_ >::clear_edge_to_parent()
	{
		m_p_edge_to_parent = 0;
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > const &
	Vertex< Kernel_ >::image() const
	{
		return *m_p_image;
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > &
	Vertex< Kernel_ >::image()
	{
		return *m_p_image;
	}	

	template< typename Kernel_ >
	inline void
	Vertex< Kernel_ >::set_image( Vertex &image_ )
	{
		m_p_image = &image_;
	}

	template< typename Kernel_ >
	inline
	Edge< Kernel_ >::Edge( Vertex< Kernel_ > &a_, Vertex< Kernel_ > &b_, double length)
		: m_a( a_ ), m_b( b_ ), m_length(length), m_canonical_loop_length( INFINITY )
	{ 
		// Initially, the image of each edge is the edge itself
		m_image.push_back( this );
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > const &
	Edge< Kernel_ >::a() const
	{
		return m_a;
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > &
	Edge< Kernel_ >::a()
	{
		return m_a;
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > const &
	Edge< Kernel_ >::b() const
	{
		return m_b;
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > &
	Edge< Kernel_ >::b()
	{
		return m_b;
	}

	template< typename Kernel_ >
	inline double
	Edge< Kernel_ >::length() const
	{
		return m_length;
	}

	template< typename Kernel_ >
	inline double
	Edge< Kernel_ >::canonical_loop_length() const
	{
		return m_canonical_loop_length;
	}

	template< typename Kernel_ >
	inline void
	Edge< Kernel_ >::set_canonical_loop_length(
		double canonical_loop_length_ )
	{
		m_canonical_loop_length = canonical_loop_length_;
	}

	template< typename Kernel_ >
	inline std::vector< Edge< Kernel_ > * > const &
	Edge< Kernel_ >::image() const
	{
		return m_image;
	}

	template< typename Kernel_ >
	inline std::vector< Edge< Kernel_ > * > &
	Edge< Kernel_ >::image()
	{
		return m_image;
	}

	template< typename Kernel_ >
	inline
	Triangle< Kernel_ >::Triangle( Edge< Kernel_ > &ab_,
		Edge< Kernel_ > &ac_, Edge< Kernel_ > &bc_ )
		: m_ab( ab_ ), m_ac( ac_ ), m_bc( bc_ )
	{
	}

	template< typename Kernel_ >
	inline Edge< Kernel_ > const &
	Triangle< Kernel_ >::ab() const
	{
		return m_ab;
	}

	template< typename Kernel_ >
	inline Edge< Kernel_ > &
	Triangle< Kernel_ >::ab()
	{
		return m_ab;
	}

	template< typename Kernel_ >
	inline Edge< Kernel_ > const &
	Triangle< Kernel_ >::ac() const
	{
		return m_ac;
	}

	template< typename Kernel_ >
	inline Edge< Kernel_ > &
	Triangle< Kernel_ >::ac()
	{
		return m_ac;
	}

	template< typename Kernel_ >
	inline Edge< Kernel_ > const &
	Triangle< Kernel_ >::bc() const
	{
		return m_bc;
	}

	template< typename Kernel_ >
	inline Edge< Kernel_ > &
	Triangle< Kernel_ >::bc()
	{
		return m_bc;
	}

	template< typename Kernel_ >
	inline
	Canonical_loop< Kernel_ >::Canonical_loop( Vertex< Kernel_ > &vertex_,
		Edge< Kernel_ > &edge_, boost::dynamic_bitset<> const &bits_ )
		: m_vertex( vertex_ ), m_edge( edge_ ), m_bits( bits_ )
	{
		set_norm( edge_.canonical_loop_length() );
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > const &
	Canonical_loop< Kernel_ >::vertex() const
	{
		return m_vertex;
	}
	
	template< typename Kernel_ >
	inline Edge< Kernel_ > const &
	Canonical_loop< Kernel_ >::edge() const
	{
		return m_edge;
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > &
	Canonical_loop< Kernel_ >::vertex()
	{
		return m_vertex;
	}
	
	template< typename Kernel_ >
	inline Edge< Kernel_ > &
	Canonical_loop< Kernel_ >::edge()
	{
		return m_edge;
	}

	template< typename Kernel_ >
	inline boost::dynamic_bitset<> const &
	Canonical_loop< Kernel_ >::bits() const
	{
		return m_bits;
	}
	
	template< typename Kernel_ >
	inline boost::dynamic_bitset<> &
	Canonical_loop< Kernel_ >::bits()
	{
		return m_bits;
	}

	template< typename Kernel_ >
	inline unsigned
	Basis_loop< Kernel_ >::size() const
	{
		return m_edges.size();
	}

	template< typename Kernel_ >
	inline typename Basis_loop< Kernel_ >::Const_iterator
	Basis_loop< Kernel_ >::begin() const
	{
		return m_edges.begin();
	}

	template< typename Kernel_ >
	inline typename Basis_loop< Kernel_ >::Const_iterator
	Basis_loop< Kernel_ >::end() const
	{
		return m_edges.end();
	}

	template< typename Kernel_ >
	inline typename Basis_loop< Kernel_ >::Iterator
	Basis_loop< Kernel_ >::begin()
	{
		return m_edges.begin();
	}

	template< typename Kernel_ >
	inline typename Basis_loop< Kernel_ >::Iterator
	Basis_loop< Kernel_ >::end()
	{
		return m_edges.end();
	}

	template< typename Kernel_ >
	inline void
	Basis_loop< Kernel_ >::push_back( Edge< Kernel_ > *p_edge_ )
	{
		m_edges.push_back( p_edge_ );
		set_norm( norm() + p_edge_->length() );
	}

	template< typename Kernel_ >
	inline
	Complex< Kernel_ >::Complex( int dimensions, bool verbose_ )
		: num_dimensions(dimensions), m_verbose( verbose_ ), m_tree_size( 0 ), m_p_progress( 0 ),
		m_expanded( false )
	{
		// resize our lower and upper bound vectors
		lower_bound.resize(dimensions);
		upper_bound.resize(dimensions);
	}

	// This is the constructor for matrix format
	template< typename Kernel_ >
	inline
	Complex< Kernel_ >::Complex( bool verbose_)
		: m_verbose( verbose_ ), m_tree_size( 0 ), m_p_progress( 0 ),
		m_expanded( false )
	{
	}

	template< typename Kernel_ >
	inline
	Complex< Kernel_ >::~Complex()
	{
		using namespace std;
		using namespace boost;

		if ( m_verbose )
		{
			cout << endl;
			cout << "Destroying simplicial complex..." << flush;
			m_timer.restart();
		}

		for_each( m_vertices.begin(), m_vertices.end(),
			checked_delete< Vertex< Kernel_ > > );

		for_each( m_edges.begin(), m_edges.end(),
			checked_delete< Edge< Kernel_ > > );
		for_each( m_triangles.begin(), m_triangles.end(),
			checked_delete< Triangle< Kernel_ > > );
		for_each( m_canonical_loops.begin(), m_canonical_loops.end(),
			checked_delete< Canonical_loop< Kernel_ > > );
		for_each( m_basis_loops.begin(), m_basis_loops.end(),
			checked_delete< Basis_loop< Kernel_ > > );

		if ( m_verbose )
			cout << "done in " << m_timer.elapsed() << "s" << endl;
	}

	template< typename Kernel_ >
	inline unsigned
	Complex< Kernel_ >::number_of_vertices() const
	{
		return m_vertices.size();
	}

	template< typename Kernel_ >
	inline Vertex< Kernel_ > &
	Complex< Kernel_ >::vertex_at( unsigned i_ )
	{
		return *m_vertices.at( i_ );
	}


	// inserts the vertix into the complex.
	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::insert_vertex( Vertex< Kernel_ > *p_vertex_)
	{
		if(input_type_flag==0)
		{
			int i;
			// updates the upper and lower bounds for the complex in each dimension
			if ( number_of_vertices() == 0 )
			{
				for(i=0; i<num_dimensions; i++)
				{
					lower_bound.at(i) = p_vertex_->location().get_coord(i);
					upper_bound.at(i) = p_vertex_->location().get_coord(i);
				}
			}
			else
			{
				for(i=0; i<num_dimensions; i++)
				{
					if(lower_bound.at(i) > p_vertex_->location().get_coord(i))
						lower_bound.at(i) = p_vertex_->location().get_coord(i);
					else if(upper_bound.at(i) < p_vertex_->location().get_coord(i))
						upper_bound.at(i) = p_vertex_->location().get_coord(i);
				}
			}
		}
		p_vertex_->set_index( m_vertices.size() );
		m_vertices.push_back( p_vertex_ );
	}

	template< typename Kernel_ >
	inline unsigned
	Complex< Kernel_ >::number_of_edges() const
	{
		return m_edges.size();
	}

	template< typename Kernel_ >
	inline Edge< Kernel_ > &
	Complex< Kernel_ >::edge_at( unsigned i_ )
	{
		return *m_edges.at( i_ );
	}

	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::create_edge( Vertex< Kernel_ > &a_,
		Vertex< Kernel_ > &b_, double length)
	{
		Edge< Kernel_ > *p_ab( new Edge< Kernel_ >( a_, b_, length) );
		p_ab->set_index( number_of_edges() );
		
		a_.coboundary().push_back( p_ab );
		b_.coboundary().push_back( p_ab );

		VV ab_key = std::make_pair( std::min( &a_, &b_ ), std::max( &a_, &b_ ) );
		// m_vv2e is an unordered map.  A pair of vertices get mapped to an edge.
		m_vv2e.insert( std::make_pair( ab_key, p_ab ) );
		m_edges.push_back( p_ab );
	}

	template< typename Kernel_ >
	inline unsigned
	Complex< Kernel_ >::number_of_triangles() const
	{
		return m_triangles.size();
	}

	template< typename Kernel_ >
	inline Triangle< Kernel_ > &
	Complex< Kernel_ >::triangle_at( unsigned i_ )
	{
		return *m_triangles.at( i_ );
	}

	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::create_triangle( Vertex< Kernel_ > &a_,
		Vertex< Kernel_ > &b_, Vertex< Kernel_ > &c_ )
	{
		VV ab_key( std::min( &a_, &b_ ), std::max( &a_, &b_ ) );
		VV ac_key( std::min( &a_, &c_ ), std::max( &a_, &c_ ) );
		VV bc_key( std::min( &b_, &c_ ), std::max( &b_, &c_ ) );

		// if edges of the triangle do not exist, create them
		if ( m_vv2e.find( ab_key ) == m_vv2e.end() )
			create_edge( a_, b_, sqrt(a_.location().get_squared_distance_to(b_.location())));
		if ( m_vv2e.find( ac_key ) == m_vv2e.end() )
			create_edge( a_, c_, sqrt(a_.location().get_squared_distance_to(c_.location())));
		if ( m_vv2e.find( bc_key ) == m_vv2e.end() )
			create_edge( b_, c_, sqrt(b_.location().get_squared_distance_to(c_.location())));

		Edge< Kernel_ > &ab( *m_vv2e[ ab_key ] );
		Edge< Kernel_ > &ac( *m_vv2e[ ac_key ] );
		Edge< Kernel_ > &bc( *m_vv2e[ bc_key ] );

		m_triangles.push_back( new Triangle< Kernel_ >( ab, ac, bc ) );
	}

	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::sample( double coefficient_ )
	{
		using namespace std;
		using namespace boost;
		using namespace math;

		if ( m_verbose )
		{
			cout << endl;
			cout << "Sampling simplicial complex..." << flush;
			m_timer.restart();
		}

		if ( coefficient_ == 1 )
		{
			// every vertex should be in the sample
			for ( unsigned i( 0 ); i < number_of_vertices(); ++i )
				vertex_at( i ).set_flag( Vertex< Kernel_ >::IS_IN_SAMPLE );

			return;
		}

		typedef vector< Vertex< Kernel_ > * > Component;
		vector< Component > components( number_of_vertices() );
		for ( unsigned i( 0 ); i != number_of_vertices(); ++i )
		{
			Vertex< Kernel_ > &vertex( vertex_at( i ) );

			Component &component( components.at( vertex.image().index() ) );
			component.push_back( &vertex );
		}
		
		for ( unsigned i( 0 ); i < components.size(); ++i )
		{
			Component &component( components.at( i ) );
			if ( component.empty() )
				continue;
			
			unsigned new_size( iround( component.size() * coefficient_ ) );
			if ( new_size == 0 )
				new_size = 1;
			
			// random sampling
			for ( unsigned j( 0 ); j != new_size; ++j )
			{				
				unsigned offset( component.size() - new_size + j );
				unsigned new_offset( rand() % ( offset + 1 ) );				

				Vertex< Kernel_ > &candidate( *component.at( new_offset ) );
				if ( candidate.has_flag( Vertex< Kernel_ >::IS_IN_SAMPLE ) )
					new_offset = offset;

				Vertex< Kernel_ > &vertex( *component.at( new_offset ) );
				vertex.set_flag( Vertex< Kernel_ >::IS_IN_SAMPLE );
			}

			// max-min sampling
			/*vector< Vertex< Kernel_ > * > sample;
			sample.reserve( new_size );
			
			for ( unsigned j( 0 ); j != new_size; ++j )
			{
				double max_min_squared_distance( 0 );
				Vertex< Kernel_ > *p_new_sample;

				for ( unsigned k( 0 ); k < component.size(); ++k )
				{
					Vertex< Kernel_ > &candidate( *component.at( k ) );
					if ( candidate.has_flag( Vertex< Kernel_ >::IS_IN_SAMPLE ) )
						continue;

					double min_squared_distance( INFINITY );

					for ( unsigned l( 0 ); l < sample.size(); ++l )
					{
						Vertex< Kernel_ > &sample( *sample.at( l ) );
						double squared_distance( CGAL::squared_distance(
							sample.location(), candidate.location() ) );

						if ( squared_distance < min_squared_distance )
							min_squared_distance = squared_distance;
					}

					if ( min_squared_distance > max_min_squared_distance )
					{
						max_min_squared_distance = min_squared_distance;
						p_new_sample = &candidate;
					}
				}

				p_new_sample->set_flag( Vertex< Kernel_ >::IS_IN_SAMPLE );
				sample.push_back( p_new_sample );
			}*/
		}

		if ( m_verbose )
			cout << "done in " << m_timer.elapsed() << "s" << endl;
	}

	template< typename Kernel_ >
	inline unsigned
	Complex< Kernel_ >::basis_rank() const
	{
		return m_basis_loops.size();
	}

	template< typename Kernel_ >
	inline Basis_loop< Kernel_ > const &
	Complex< Kernel_ >::basis_loop_at( unsigned i_ ) const
	{
		return *m_basis_loops.at( i_ );
	}

	template< typename Kernel_ >
	inline Basis_loop< Kernel_ > &
	Complex< Kernel_ >::basis_loop_at( unsigned i_ )
	{
		return *m_basis_loops.at( i_ );
	}

	//Builds the rips complex with parameter alpha
	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::expand( double alpha_ )
	{
		using namespace std;
		using namespace boost;
		using namespace math;
		using namespace Headers;

		if ( m_verbose )		
			m_p_progress = new progress_display(number_of_vertices() );					

		double squared_alpha( alpha_ * alpha_ );
		vector<double> width(num_dimensions);
		for(int i=0; i<num_dimensions; i++)
			width.at(i) = upper_bound.at(i) - lower_bound.at(i);

		vector<int> resolution(num_dimensions);
		unsigned max_resolution( 50 );	// sets a cap on the number of cells in the grid
		for(int i=0; i<num_dimensions; i++)
		{
			// number of cells in the grid.  Smaller alpha will make more cells.
			resolution.at(i) = iround( floor( width.at(i) / alpha_ ) );
			if ( resolution.at(i) == 0 )
				resolution.at(i) = 1;
			else if(resolution.at(i) > max_resolution)
				resolution.at(i) = max_resolution;

		}

		typedef Grid< Kernel_, Vertex< Kernel_ > * > Vertex_grid;
		typedef typename Vertex_grid::Cell Vertex_cell;

		// This makes a grid of equally spaced points that enclose our point cloud.  It will be used to calculate the Rips complex.
		Vertex_grid grid( lower_bound, upper_bound, resolution, num_dimensions );

		// vector of cells.  
		vector< Vertex_cell * > vertex_cells;	// only has the cells with vertices in them.
		vertex_cells.reserve( number_of_vertices() );

		// fill out the grid
		for ( unsigned i( 0 ); i < number_of_vertices(); ++i )
		{
			Vertex< Kernel_ > &vertex( vertex_at( i ) );	// this is just the vertex at index i
			// this just adds the vertex to the cell's vector of vertices (data) and returns
			// the cell that this vertex belongs to.
			Vertex_cell &cell( grid.add_at( &vertex, vertex.location() ) );	
			vertex_cells.push_back( &cell );
		}	

		int num_neighbors = pow(3, num_dimensions);
		
		int dimension_index;
		// this will store the indices of the neighbor in each dimension.
		// -1 corresponds to the "left" neighbor, 0 is the one that shares the
		// same dimension as the current cell, and 1 is the cell to the "right"
		vector<int> current_neighbor(num_dimensions);

		// create edges
		for ( unsigned i( 0 ); i < number_of_vertices(); ++i )
		{
			
			Vertex< Kernel_ > &a( vertex_at( i ) );
			Vertex_cell &a_cell( *vertex_cells.at( i ) );

			for(int j=0; j<num_dimensions; j++)
				current_neighbor.at(j) = -1;
			
			for ( unsigned j( 0 ); j < num_neighbors; ++j )
			{
				if(j != 0)
				{
					dimension_index = 0;
					current_neighbor.at(dimension_index)++;
					while(current_neighbor.at(dimension_index) > 1)
					{
						current_neighbor.at(dimension_index) = -1;
						dimension_index++;
						current_neighbor.at(dimension_index)++;
					}
				}

				if ( !a_cell.has_neighbor_at( current_neighbor ) )
					continue;

				Vertex_cell &b_cell( a_cell.neighbor_at( current_neighbor ) );
				typename Vertex_cell::Data_iterator it_data( b_cell.datas_begin() );
				// iterates through each cell's data
				for ( ; it_data != b_cell.datas_end(); ++it_data )
				{
					Vertex< Kernel_ > &b( **it_data );

					// if b.index <= a.index, we already checked this one
					if ( b.index() <= a.index() )
						continue;

					// if their distance is greater than the parameter specified to build the 
					// rips complex, we don't create an edge between them.
					if ( a.location().get_squared_distance_to(b.location()) > squared_alpha )
						continue;

					// Creates an edge between a and b since we determined it should be in the complex.
					create_edge(a, b, sqrt(a.location().get_squared_distance_to(b.location())));
				}
			}
			if ( m_verbose )
				++( *m_p_progress );
		}		
		
		if ( m_verbose )
			delete m_p_progress;			

		m_expanded = true;
	}

	//Builds the rips complex with parameter alpha for matrix input
	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::expand_matrix( double alpha_ )
	{
		using namespace std;
		using namespace boost;
		using namespace math;
		using namespace Headers;

		if ( m_verbose )		
			m_p_progress = new progress_display(number_of_vertices() );

		// create edges
		for ( unsigned i( 0 ); i < number_of_vertices(); ++i )
		{
			Vertex< Kernel_ > &a( vertex_at( i ) );
			for( unsigned j = i+1; j < number_of_vertices(); j++)
			{
				Vertex< Kernel_ > &b( vertex_at( j ) );
				if(g_distance_matrix[i][j] < alpha_)
					create_edge(a, b, g_distance_matrix[i][j]);
			}
			if ( m_verbose )
				++( *m_p_progress );
		}		
		
		if ( m_verbose )
			delete m_p_progress;

		m_expanded = true;
	}

	//Builds the rips complex with parameter alpha for graph input
	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::expand_graph( double alpha_ )
	{
		using namespace std;
		using namespace boost;
		using namespace math;
		using namespace Headers;

		// We just add all edges in the graph to the complex if alpha was not specified.
		if(alpha_ == INFINITY)
		{
			for(int i=0; i<g_adjacency_list.size(); i++)
			{
				Vertex< Kernel_ > &a( vertex_at( i ) );
				for(int j=0; j<g_adjacency_list[i].size(); j++)
				{
					int v2 = g_adjacency_list[i][j].first;
					Vertex< Kernel_ > &b(vertex_at(v2));
					if(i < v2)
						create_edge(a, b, g_adjacency_list[i][j].second);
				}
			}
		}
		else
		{				
			if ( m_verbose )		
				m_p_progress = new progress_display(number_of_vertices() );

			vector<double> distances(number_of_vertices());

			// builds the rips complex by getting the shortest distances to each point and comparing the distances to the parameter alpha.
			for(int i=0; i<number_of_vertices(); i++)
			{
				// gets the shortest path from node i to all other nodes in the graph.
				distances = shortest_path_graph(i, alpha_);

				Vertex< Kernel_ > &a( vertex_at( i ) );
				for(int j=i+1; j<number_of_vertices(); j++)
				{
					if(distances[j] < alpha_)
					{
						Vertex< Kernel_ > &b( vertex_at( j ) );
						create_edge(a, b, distances[j]);
					}
				}
				if ( m_verbose )
					++( *m_p_progress );
			}
		
			if ( m_verbose )
				delete m_p_progress;
		}

		m_expanded = true;
	}

	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::contract_if_needed( Edge< Kernel_ > &edge_ )
	{
		// contracts one vertex of the edge
		unsigned a_depth( 0 );
		Vertex< Kernel_ > *p_a_image( &edge_.a() );
		while ( &p_a_image->image() != p_a_image )
		{
			p_a_image = &p_a_image->image();
			++a_depth;
		}

		// contracts the other edge
		unsigned b_depth( 0 );
		Vertex< Kernel_ > *p_b_image( &edge_.b() );
		while ( &p_b_image->image() != p_b_image )
		{
			p_b_image = &p_b_image->image();
			++b_depth;
		}

		if ( p_a_image != p_b_image )
		{
			if ( a_depth < b_depth )
				p_a_image->set_image( *p_b_image );
			else
				p_b_image->set_image( *p_a_image );

			edge_.image().clear();
		}
	}

	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::contract_if_needed( Triangle< Kernel_ > &triangle_ )
	{
		using namespace std;
		using namespace Headers;

		vector< Edge< Kernel_ > * > current_boundary;

		deque< Edge< Kernel_ > * > edges;
		edges.push_back( &triangle_.ab() );
		edges.push_back( &triangle_.ac() );
		edges.push_back( &triangle_.bc() );

		sort( edges.begin(), edges.end(), pointee_index_is_less );

		while ( !edges.empty() )
		{
			Edge< Kernel_ > &edge( *edges.front() );
			edges.pop_front();

			vector< Edge< Kernel_ > * > &image( edge.image() );
			if ( image.empty() )
				continue;

			if ( image.size() == 1 && image.front() == &edge )
			{
				vector< Edge< Kernel_ > * > edge_chain( 1, &edge );

				vector< Edge< Kernel_ > * > buffer;
				set_symmetric_difference( current_boundary.begin(),
					current_boundary.end(), edge_chain.begin(),
					edge_chain.end(), back_inserter( buffer ),
					pointee_index_is_less );
				swap( current_boundary, buffer );

				continue;
			}

			copy( image.begin(), image.end(), back_inserter( edges ) );
		}

		if ( !current_boundary.empty() )
		{
			Edge< Kernel_ > &contracted_edge( *current_boundary.back() );

			vector< Edge< Kernel_ > * > buffer;
			set_symmetric_difference( contracted_edge.image().begin(),
				contracted_edge.image().end(), current_boundary.begin(),
				current_boundary.end(), back_inserter( buffer ),
				pointee_index_is_less );
			swap( contracted_edge.image(), buffer );
		}
	}

	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::contract()
	{
		using namespace std;
		using namespace boost;
		using namespace Headers;

		if ( m_verbose )
		{
			cout << endl;
			cout << "Processing simplicial complex..." << flush;
			
			if ( m_expanded )
			{
				m_p_progress = new progress_display( 3 * number_of_edges()
					+ number_of_vertices() );
			}
			else
			{
				m_p_progress = new progress_display( 2 * number_of_edges()
					+ number_of_vertices() + number_of_triangles() );
			}

			m_timer.restart();						
		}

		// We collapse all edges that do not have empty boundary in the
		// contracted complex; we employ union-find algorithm here; eventually,
		// only one vertex per component survive in the contracted complex.
		// This gives us the root of a tree for each component
		for ( unsigned i( 0 ); i != number_of_edges(); ++i )
		{
			contract_if_needed( edge_at( i ) );		

			if ( m_verbose )
				++( *m_p_progress );
		}		

		// Flatten the vertex images
		// set each vetex to its contracted vertex, does recursively till recieve a vertex not contracted
		for ( unsigned i( 0 ); i != number_of_vertices(); ++i )
		{			
			Vertex< Kernel_ > &vertex( vertex_at( i ) );

			Vertex< Kernel_ > *p_image( &vertex.image() );
			while ( &p_image->image() != p_image )
				p_image = &p_image->image();

			vertex.set_image( *p_image );

			if ( m_verbose )
				++( *m_p_progress );
		}		

		// We collapse all triangles that have non-empty boundary in the
		// current contracted complex; whenever a triangle is collapsed,
		// we map the youngest edge in its boundary to the remaining edges
		if ( m_expanded )
		{			
			for ( unsigned i( 0 ); i < number_of_edges(); ++i )
			{
				Edge< Kernel_ > &ab( edge_at( i ) );

				Vertex< Kernel_ > &a( ab.a() );
				Vertex< Kernel_ > &b( ab.b() );				

				for ( unsigned j( 0 ); j < a.coboundary().size(); ++j )
				{
					Edge< Kernel_ > &ac( *a.coboundary().at( j ) );
					if ( ac.index() <= ab.index() )
						continue;

					Vertex< Kernel_ > &c( a.coneighbor( ac ) );					

					for ( unsigned k( 0 ); k < b.coboundary().size(); ++k )
					{
						Edge< Kernel_ > &bc( *b.coboundary().at( k ) );
						if ( bc.index() <= ab.index() )
							continue;

						if ( &b.coneighbor( bc ) != &c )
							continue;

						Triangle< Kernel_ > triangle( ab, bc, ac );
						contract_if_needed( triangle );
					}
				}

				if ( m_verbose )
					++( *m_p_progress );
			}			
		}
		else
		{
			for ( unsigned i( 0 ); i != number_of_triangles(); ++i )
			{
				contract_if_needed( triangle_at( i ) );

				if ( m_verbose )
					++( *m_p_progress );
			}
		}		

		m_e2b.resize( number_of_edges(), -1 );
		unsigned basis_size( 0 );

		// For each edge, we compute its image in the contracted complex
		for ( unsigned i( 0 ); i != number_of_edges(); ++i )
		{
			vector< Edge< Kernel_ > * > new_image;			
			
			deque< Edge< Kernel_ > * > edges;
			edges.push_back( &edge_at( i ) );

			while ( !edges.empty() )
			{
				Edge< Kernel_ > &edge( *edges.front() );
				edges.pop_front();
				
				vector< Edge< Kernel_ > * > &image( edge.image() );
				if ( image.empty() )
					continue;

				if ( image.size() == 1 && image.front() == &edge )
				{
					vector< Edge< Kernel_ > * > edge_chain( 1, &edge );
					vector< Edge< Kernel_ > * > buffer;
					set_symmetric_difference( new_image.begin(),
						new_image.end(), edge_chain.begin(),
						edge_chain.end(), back_inserter( buffer ),
						pointee_index_is_less );
					swap( new_image, buffer );

					continue;
				}

				copy ( image.begin(), image.end(), back_inserter( edges ) );
			}

			edge_at( i ).image() = new_image;

			// basis edge is mapped to itself
			if ( new_image.size() == 1 && new_image.front() == &edge_at( i ) )
				m_e2b.at( i ) = basis_size++;

			if ( m_verbose )
				++( *m_p_progress );
		}

		m_basis_loops.resize( basis_size, 0 );

		if ( m_verbose )
		{
			delete m_p_progress;

			cout << "done in " << m_timer.elapsed() << "s" << endl;			
		}

		unsigned number_of_components( 0 );
		for ( unsigned i( 0 ); i != number_of_vertices(); ++i )
		{
			Vertex< Kernel_ > const &vertex( vertex_at( i ) );
			if ( &vertex.image() == &vertex )
				++number_of_components;
		}

		if ( m_verbose )
		{
			cout << endl;
			cout << "Number of components: " << number_of_components << endl;
		}
	}

	template< typename Kernel_ >
	inline unsigned
	Complex< Kernel_ >::e2b( unsigned edge_index_ ) const
	{
		return m_e2b.at( edge_index_ );
	}

	// Prints shortest paths from src to all other vertices using Dijkstra's algorithm
	template< typename Kernel_ >
	inline vector<double>
	Complex< Kernel_ >::shortest_path_graph(int src, double alpha)
	{
		// Create a set to store vertices that are being
		// prerocessed
		set< pair<double, int> > setds;

		// Create a vector for distances and initialize all
		// distances as infinite (INFINITY)
		vector<double> dist(m_vertices.size(), INFINITY);

		// Insert source itself in Set and initialize its
		// distance as 0.
		setds.insert(make_pair(0, src));
		dist[src] = 0;

		/* Looping till all shortest distance are finalized
		then setds will become empty */
		while (!setds.empty())
		{
			// The first vertex in Set is the minimum distance
			// vertex, extract it from set.
			pair<double, int> tmp = *(setds.begin());
			setds.erase(setds.begin());

			// vertex label is stored in second of pair (it
			// has to be done this way to keep the vertices
			// sorted distance (distance must be first item
			// in pair)
			int u = tmp.second;

			// 'i' is used to get all adjacent vertices of a vertex
			vector< pair<int, double> >::iterator i;
			for (i = g_adjacency_list[u].begin(); i != g_adjacency_list[u].end(); ++i)
			{
				// Get vertex label and weight of current adjacent
				// of u.
				int v = (*i).first;
				double weight = (*i).second;

				if(dist[u] + weight > alpha)
					continue;

				//  If there is shorter path to v through u.
				if (dist[v] > dist[u] + weight)
				{
					/*  If distance of v is not INF then it must be in
					our set, so removing it and inserting again
					with updated less distance.  
					Note : We extract only those vertices from Set
					for which distance is finalized. So for them, 
					we would never reach here.  */
					if (dist[v] != INFINITY)
						setds.erase(setds.find(make_pair(dist[v], v)));

					// Updating distance of v
					dist[v] = dist[u] + weight;
					setds.insert(make_pair(dist[v], v));
				}
			}
		}
 
		// returns the distance vector
		return dist;
	}
	
	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::compute_basis()
	{
		using namespace std;
		using namespace CGAL;
		using namespace boost;
		using namespace Headers;						

		if ( m_verbose )
		{
			cout << endl;
			cout << "Basis rank: " << basis_rank() << endl;
		}

		if ( basis_rank() == 0 )
			return;

		// COMPUTES THE SAMPLE SIZE OF THE COMPUTED COMPLEX
		unsigned sample_size( 0 );
		for ( unsigned i( 0 ); i != number_of_vertices(); ++i )
		{
			if ( vertex_at( i ).has_flag( Vertex< Kernel_ >::IS_IN_SAMPLE ) )
				++sample_size;
		}

		if ( m_verbose )
		{
			cout << endl;
			cout << "Computing canonical loops..." << flush;

			m_p_progress = new progress_display( sample_size );
			m_timer.restart();
		}		

		// we can have at most Vb_2 canonical loops (less if >1 components)
		cout<<"compute basis: "<<basis_rank()<<"\n";

		m_canonical_loops.reserve( sample_size * basis_rank() );

		for ( unsigned i( 0 ); i != number_of_vertices(); ++i )
				vertex_at( i ).contracted_path_to_root().resize( basis_rank() );

		for ( unsigned i( 0 ); i != number_of_vertices(); ++i )
		{
			Vertex< Kernel_ > &vertex( vertex_at( i ) );	
			if ( !vertex.has_flag( Vertex< Kernel_ >::IS_IN_SAMPLE ) )
				continue;
			
			compute_canonical_loops_for( vertex );
			if ( m_verbose )
				++( *m_p_progress );
		}
		if ( m_verbose )
		{
			cout << m_canonical_loops.size() << " canonical loops computed in "
				<< m_timer.elapsed() << "s" << endl;

			delete m_p_progress;
		}

		if ( m_verbose )
		{
			cout << endl;
			cout << "Computing independent loops..." << flush;

			m_p_progress = new progress_display( m_canonical_loops.size() );
			m_timer.restart();
		}

		sort( m_canonical_loops.begin(), m_canonical_loops.end(),
			pointee_norm_is_less );

		typedef Vector_d< Homogeneous_d< int > > Vector_d;
		vector< Vector_d > basis;
		basis.reserve( basis_rank() );

		unsigned i( 0 );
		
		for ( ; i != m_canonical_loops.size(); ++i )
		{
			Canonical_loop< Kernel_ > &canonical( *m_canonical_loops.at( i ) );

			// Test if independent; if no, proceed

			dynamic_bitset<> &bits( canonical.bits() );

			// static 
			vector< int > coeffs( basis_rank() );
			for ( unsigned j( 0 ); j != basis_rank(); ++j )
				coeffs.at( j ) = ( bits.test( j ) ? 1 : 0 );
			
			Vector_d current( basis_rank(), coeffs.begin(), coeffs.end() );
			basis.push_back( current );

			if ( !linearly_independent( basis.begin(), basis.end() ) )
			{
				basis.pop_back();

				if ( m_verbose )
					++( *m_p_progress );

				continue;
			}

			Basis_loop< Kernel_ > *p_basis_loop( new Basis_loop< Kernel_ > );

			Vertex< Kernel_ > &vertex( canonical.vertex() );
			Edge< Kernel_ > &edge( canonical.edge() );

			compute_shortest_path_tree_for( vertex );

			// Compute actual canonical loop in the initial complex

			Vertex< Kernel_ > *p_a( &edge.a() );
			while ( p_a->has_parent() )
			{
				p_a->edge_to_parent().toggle_flag(
					Edge< Kernel_ >::IS_IN_LOOP );
				p_a = &p_a->parent();
			}

			Vertex< Kernel_ > *p_b( &edge.b() );
			while ( p_b->has_parent() )
			{
				p_b->edge_to_parent().toggle_flag(
					Edge< Kernel_ >::IS_IN_LOOP );
				p_b = &p_b->parent();
			}

			p_a = &edge.a();
			while ( p_a->has_parent() )
			{
				Edge< Kernel_ > &edge_to_parent( p_a->edge_to_parent() );
				if ( edge_to_parent.has_flag( Edge< Kernel_ >::IS_IN_LOOP ) )
				{
					edge_to_parent.clear_flag( Edge< Kernel_ >::IS_IN_LOOP );
					p_basis_loop->push_back( &edge_to_parent );
				}
				p_a = &p_a->parent();
			}

			p_b = &edge.b();
			while ( p_b->has_parent() )
			{
				Edge< Kernel_ > &edge_to_parent( p_b->edge_to_parent() );
				if ( edge_to_parent.has_flag( Edge< Kernel_ >::IS_IN_LOOP ) )
				{
					edge_to_parent.clear_flag( Edge< Kernel_ >::IS_IN_LOOP );
					p_basis_loop->push_back( &edge_to_parent );
				}
				p_b = &p_b->parent();
			}

			p_basis_loop->push_back( &edge );
			m_basis_loops.at( basis.size() - 1 ) = p_basis_loop;

			clear_shortest_path_tree();

			if ( basis.size() == basis_rank() )
				break;

			if ( m_verbose )
				++( *m_p_progress );
		}

		if ( m_verbose )
		{
			for ( ; i != m_canonical_loops.size(); ++i )
				++( *m_p_progress );

			delete m_p_progress;

			cout << m_basis_loops.size() << " independent loops computed in "
				<< m_timer.elapsed() << "s" << endl << endl;
		}
	}

	template< typename Kernel_ >
	inline bool
	pointee_canonical_loop_length_is_less( Edge< Kernel_ > *p_a_,
		Edge< Kernel_ > *p_b_ )
	{
		return p_a_->canonical_loop_length() < p_b_->canonical_loop_length();
	}

	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::compute_canonical_loops_for(
		Vertex< Kernel_ > &vertex_ )
	{
		using namespace std;
		using namespace CGAL;
		using namespace boost;
		compute_shortest_path_tree_for( vertex_ );
		compute_canonical_loop_lengths();
		
		sort( m_edges.begin(), m_edges.end(),
			pointee_canonical_loop_length_is_less< Kernel_ > );
		typedef Vector_d< Homogeneous_d< int > > Vector_d;
		vector< Vector_d > basis;
		basis.reserve( basis_rank() );
		// static dynamic_bitset<> bits;
		for ( unsigned i( m_tree_size ); i != number_of_edges(); ++i )
		{
			Edge< Kernel_ > &edge( edge_at( i ) );
			if ( edge.canonical_loop_length() == INFINITY )
				continue;

			// Represents the canonical loop in the contracted complex
			dynamic_bitset<> bits( basis_rank() );

			// cout<<edge.a().index()<<" basis_rank: "<<basis_rank()<<" bits: "<<bits<<"\n";	
			bits.reset();

			bits ^= edge.a().contracted_path_to_root();

			bits ^= edge.b().contracted_path_to_root();
			
			vector< Edge< Kernel_ > * > &image( edge.image() );
			for ( unsigned i( 0 ); i != image.size(); ++i )
				bits.flip( e2b( image.at( i )->index() ) );

			if ( bits.none() )
				continue;

			// Test if independent; if no, proceed

			vector< int > coeffs( basis_rank() );
			for ( unsigned i( 0 ); i != basis_rank(); ++i )
				coeffs.at( i ) = ( bits.test( i ) ? 1 : 0 );

			Vector_d current( basis_rank(), coeffs.begin(), coeffs.end() );
			basis.push_back( current );

			if ( !linearly_independent( basis.begin(), basis.end() ) )
			{
				basis.pop_back();
				continue;
			}

			// New independent canonical loop found
			m_canonical_loops.push_back( new Canonical_loop< Kernel_ >(
				vertex_, edge, bits ) );

			if ( basis.size() == basis_rank() )
				break;

		}

		clear_shortest_path_tree();
	}

	template< typename Vertex_ >
	struct Pointee_distance_to_root_is_less
	{
		bool operator()( Vertex_ const *p_a_, Vertex_ const *p_b_ ) const
		{
			return p_a_->distance_to_root() < p_b_->distance_to_root();
		}
	};
	
	template< typename Vertex_ >
	struct Property_map	
		: public boost::put_get_helper< size_t, Property_map< Vertex_ > >
	{
		typedef Vertex_ *key_type;
		typedef size_t value_type;
		typedef size_t reference;
		typedef boost::readable_property_map_tag category;

		value_type operator[]( key_type key_ ) const
		{
			return key_->index();
		}
	};

	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::compute_shortest_path_tree_for(
		Vertex< Kernel_ > &vertex_ )
	{
		using namespace std;
		using namespace boost;

		m_tree_size = 0;

		typedef vector< Vertex< Kernel_ > * > Container;
		typedef Pointee_distance_to_root_is_less< Vertex< Kernel_ > > Less;
		typedef Property_map< Vertex< Kernel_ > > Map;

		// we need a priority queue that allows changing the element value
		typedef mutable_queue< Vertex< Kernel_ > *, Container,
			Less, Map > Vertex_queue;
		Vertex_queue vertex_queue( m_vertices.size(), Less(), Map() ) ;

		// Build the shortest path tree using Dijkstra algorithm

		for ( unsigned i( 0 ); i != number_of_vertices(); ++i )
			vertex_at( i ).set_distance_to_root( INFINITY );

		vertex_.set_distance_to_root( 0 );

		vertex_queue.push( &vertex_ );
		vertex_.set_flag( Vertex< Kernel_ >::IS_IN_QUEUE );

		while ( !vertex_queue.empty() )
		{
			Vertex< Kernel_ > &current( *vertex_queue.top() );

			// INFTY means current belongs to another component
			if ( current.distance_to_root() == INFINITY )
				break;

			vertex_queue.pop();
			current.clear_flag( Vertex< Kernel_ >::IS_IN_QUEUE );
			current.set_flag( Vertex< Kernel_ >::IS_IN_TREE );

			if ( current.has_parent() )
			{
				current.contracted_path_to_root() =
					current.parent().contracted_path_to_root();

				Edge< Kernel_ > &edge_to_parent( current.edge_to_parent() );
				edge_to_parent.set_flag( Edge< Kernel_ >::IS_IN_TREE );
				++m_tree_size;

				for ( unsigned i( 0 ); i != edge_to_parent.image().size(); ++i )
				{					
					Edge< Kernel_ > &image( *edge_to_parent.image().at( i ) );
					unsigned index( e2b( image.index() ) );
					assert( index != -1 );
					current.contracted_path_to_root().flip( index );
				}
			}

			for ( unsigned i( 0 ); i != current.coboundary().size(); ++i )
			{
				Edge< Kernel_ > &edge( *current.coboundary().at( i ) );
				Vertex< Kernel_ > &neighbor( &edge.a() == &current ?
					edge.b() : edge.a() );

				if ( neighbor.has_flag( Vertex< Kernel_ >::IS_IN_TREE ) )
					continue;

				double distance_to_root( current.distance_to_root()
					+ edge.length() );
				if ( distance_to_root < neighbor.distance_to_root() )
				{
					neighbor.set_distance_to_root( distance_to_root );
					neighbor.set_edge_to_parent( edge );
				}

				if ( neighbor.has_flag( Vertex< Kernel_ >::IS_IN_QUEUE ) )
					vertex_queue.update( &neighbor );
				else
				{
					neighbor.set_flag( Vertex< Kernel_ >::IS_IN_QUEUE );
					vertex_queue.push( &neighbor );
				}
			}
		}
	}

	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::clear_shortest_path_tree()
	{
		for ( unsigned i( 0 ); i != number_of_vertices(); ++i )
		{
			Vertex< Kernel_ > &vertex( vertex_at( i ) );

			vertex.clear_flag( Vertex< Kernel_ >::IS_IN_TREE );
			vertex.clear_edge_to_parent();

			vertex.contracted_path_to_root().reset();
		}

		for ( unsigned i( 0 ); i != number_of_edges(); ++i )
		{
			Edge< Kernel_ > &edge( edge_at( i ) );
			edge.clear_flag( Edge< Kernel_ >::IS_IN_TREE );
			edge.set_canonical_loop_length( INFINITY );
		}
	}

	template< typename Kernel_ >
	inline void
	Complex< Kernel_ >::compute_canonical_loop_lengths()
	{
		for ( unsigned i( 0 ); i != number_of_edges(); ++i )
		{
			Edge< Kernel_ > &edge( edge_at( i ) );

			if ( edge.has_flag( Edge< Kernel_ >::IS_IN_TREE ) )
				edge.set_canonical_loop_length( 0 );
			else
			{
				Vertex< Kernel_ > &a( edge.a() ), &b( edge.b() );

				if ( a.distance_to_root() == INFINITY )
					edge.set_canonical_loop_length( INFINITY );
				else if ( b.distance_to_root() == INFINITY )
					edge.set_canonical_loop_length( INFINITY );
				else
				{
					edge.set_canonical_loop_length( edge.length() +
						a.distance_to_root() + b.distance_to_root() );
				}
			}
		}
	}
}

#endif // HOM_BAS_COMPLEX_H
