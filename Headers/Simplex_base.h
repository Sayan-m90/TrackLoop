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

#ifndef HEADERS_SIMPLEX_BASE_H
#define HEADERS_SIMPLEX_BASE_H

#include <CGAL/Point_3.h>

#include <Nothing.h>
#include <Polymorphic.h>
#include <Noncopyable.h>
#include <Indexed.h>
#include <Flagged.h>

#include <Chain.h>

namespace Headers
{	
	template< typename Kernel_, typename Simplex_, typename Base_ = Nothing >
	class Simplex_base
		: public Base_, public Noncopyable, public Polymorphic,
		public Indexed, public Flagged
	{	

	public:

		Simplex_base();
		Simplex_base( CGAL::Point_3< Kernel_ > const &location_ );

		std::string name() const;
		unsigned dimension() const;
		CGAL::Point_3< Kernel_ > const &location() const;
		void set_location( CGAL::Point_3< Kernel_ > const &location_ );
		Chain< Simplex_ > const &boundary() const;
		Chain< Simplex_ > &boundary();

		bool is_face_of( Simplex_ const &simplex_ ) const;

	protected:		

		CGAL::Point_3< Kernel_ > m_location;
		Chain< Simplex_ > m_boundary;
	};	

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline
	Simplex_base< Kernel_, Simplex_, Base_ >::Simplex_base()
	{
	}

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline
	Simplex_base< Kernel_, Simplex_, Base_ >::Simplex_base(
		CGAL::Point_3< Kernel_ > const &location_ )
		: m_location( location_ )
	{		
	}	

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline std::string
	Simplex_base< Kernel_, Simplex_, Base_ >::name() const
	{
		using namespace std;

		if ( dimension() == 0 )
			return string( 1, 'a' + index() );

		if ( dimension() == 1 )		
			return boundary().front().name() + boundary().back().name();
		
		string result;

		result.push_back( '{' );
		typename Chain< Simplex_ >::Const_iterator it_face(
			boundary().begin() );
		for ( ; it_face != boundary().end(); ++it_face )
		{
			if ( it_face != boundary().begin() )
				result.push_back( ',' );
			result.append( ( **it_face ).name() );
		}
		result.push_back( '}' );

		return result;
	}	

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline unsigned
	Simplex_base< Kernel_, Simplex_, Base_ >::dimension() const
	{
		return boundary().empty() ? 0 : boundary().size() - 1;
	}

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline CGAL::Point_3< Kernel_ > const &
	Simplex_base< Kernel_, Simplex_, Base_ >::location() const
	{
		return m_location;
	}

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline void
	Simplex_base< Kernel_, Simplex_, Base_ >::set_location(
		CGAL::Point_3< Kernel_ > const &location_ )
	{
		m_location = location_;
	}

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline Chain< Simplex_ > const &
	Simplex_base< Kernel_, Simplex_, Base_ >::boundary() const
	{
		return m_boundary;
	}

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline Chain< Simplex_ > &
	Simplex_base< Kernel_, Simplex_, Base_ >::boundary()
	{
		return m_boundary;
	}	

	template< typename Kernel_, typename Simplex_, typename Base_ >
	inline bool
	Simplex_base< Kernel_, Simplex_, Base_ >::is_face_of(
		Simplex_ const &simplex_ ) const
	{
		typename Chain< Simplex_ >::Const_iterator it_face(
			simplex_.boundary().begin() );
		for ( ; it_face != simplex_.boundary().end(); ++it_face )
		{
			if ( *it_face == this )
				return true;
		}

		return false;
	}
}

#endif // COMMON_SIMPLEX_BASE_H
