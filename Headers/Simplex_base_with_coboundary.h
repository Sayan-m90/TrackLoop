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

#ifndef HEADERS_SIMPLEX_BASE_WITH_COBOUNDARY_H
#define HEADERS_SIMPLEX_BASE_WITH_COBOUNDARY_H

namespace Headers
{
	template< typename Kernel_, typename Simplex_, typename Base_ = Nothing >
	class Simplex_base_with_coboundary		
		: public Base_
	{

	public:		

		Chain< Simplex_ > const &coboundary() const;
		Chain< Simplex_ > &coboundary();

		Simplex_ const &neighbor( Simplex_ &face_ ) const;
		Simplex_ &neighbor( Simplex_ &face_ );

		Simplex_ const &coneighbor( Simplex_ &coface_ ) const;
		Simplex_ &coneighbor( Simplex_ &coface_ );
		
	private:

		Chain< Simplex_ > m_coboundary;
	};	

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline Chain< Simplex_ > const &
	Simplex_base_with_coboundary< Kernel_, Simplex_, Base_ >::coboundary() const
	{
		return m_coboundary;
	}	

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline Chain< Simplex_ > &
	Simplex_base_with_coboundary< Kernel_, Simplex_, Base_ >::coboundary()
	{
		return m_coboundary;
	}	

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline Simplex_ const &
	Simplex_base_with_coboundary< Kernel_, Simplex_, Base_ >::neighbor(
		Simplex_ &face_ ) const
	{
		Simplex_base_with_coboundary< Kernel_, Simplex_, Base_ > const *const_this( this );
		return const_this->neighbor( face_ );		
	}

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline Simplex_ &
	Simplex_base_with_coboundary< Kernel_, Simplex_, Base_ >::neighbor(
		Simplex_ &face_ )
	{
		Chain< Simplex_ > const &coboundary( face_.coboundary() );
		typename Chain< Simplex_ >::Const_iterator it_neighbor( coboundary.begin() );
		while ( *it_neighbor != this )
			++it_neighbor;
		
		++it_neighbor;
		if ( it_neighbor == coboundary.end() )
			it_neighbor = coboundary.begin();

		return **it_neighbor;
	}

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline Simplex_ const &
	Simplex_base_with_coboundary< Kernel_, Simplex_, Base_ >::coneighbor(
		Simplex_ &coface_ ) const
	{
		Simplex_base_with_coboundary< Kernel_, Simplex_, Base_ > const *const_this( this );
		return const_this->coneighbor( coface_ );		
	}

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline Simplex_ &
	Simplex_base_with_coboundary< Kernel_, Simplex_, Base_ >::coneighbor(
		Simplex_ &coface_ )
	{		
		Chain< Simplex_ > const &boundary( coface_.boundary() );
		typename Chain< Simplex_ >::Const_iterator it_coneighbor( boundary.begin() );
		while ( *it_coneighbor != this )
			++it_coneighbor;
		
		++it_coneighbor;
		if ( it_coneighbor == boundary.end() )
			it_coneighbor = boundary.begin();

		return **it_coneighbor;
	}	
}

#endif // COMMON_SIMPLEX_BASE_WITH_COBOUNDARY_H

