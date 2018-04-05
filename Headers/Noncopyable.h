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

#ifndef NONCOPYABLE_H
#define NONCOPYABLE_H

#include <cassert>

namespace Headers
{
	class Noncopyable
	{

	public:

		Noncopyable();

	private:

		Noncopyable( Noncopyable const & );

		Noncopyable &operator=( Noncopyable const & );
	};

	inline
	Noncopyable::Noncopyable()
	{		
	}

	inline
	Noncopyable::Noncopyable( Noncopyable const & )
	{
		assert( !"Noncopyable::Noncopyable" );
	}

	inline Noncopyable &
	Noncopyable::operator=( Noncopyable const & )
	{
		assert( !"Noncopyable::operator=" );

		return *this;
	}
}

#endif // NONCOPYABLE_H
