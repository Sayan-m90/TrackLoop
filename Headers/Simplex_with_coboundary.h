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

#ifndef HEADERS_SIMPLEX_WITH_COBOUNDARY_H
#define HEADERS_SIMPLEX_WITH_COBOUNDARY_H

#include <Simplex_base.h>
#include <Simplex_base_with_coboundary.h>

namespace Headers
{
	template< typename Kernel_ >
	class Simplex_with_coboundary
		: public Simplex_base< Kernel_, Simplex_with_coboundary<  Kernel_ >,
			Simplex_base_with_coboundary< Kernel_,
				Simplex_with_coboundary<  Kernel_ > > >
	{	
		typedef Simplex_base< Kernel_, Simplex_with_coboundary<  Kernel_ >,
			Simplex_base_with_coboundary< Kernel_,
				Simplex_with_coboundary<  Kernel_ > > > Base;

	public:

		Simplex_with_coboundary();
		Simplex_with_coboundary( CGAL::Point_3< Kernel_ > const &location_ );
	};

	template< typename Kernel_ >
	inline
	Simplex_with_coboundary< Kernel_ >::Simplex_with_coboundary()
		: Base()
	{
	}

	template< typename Kernel_ >
	inline
	Simplex_with_coboundary< Kernel_ >::Simplex_with_coboundary(
		CGAL::Point_3< Kernel_ > const &location_ )
		: Base( location_ )
	{
	}
}

#endif // COMMON_SIMPLEX_WITH_COBOUNDARY_H
