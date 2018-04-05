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

#ifndef HEADERS_NONCOMPARABLE_H
#define HEADERS_NONCOMPARABLE_H

namespace Headers
{
	class Noncomparable
	{

	private:

		bool operator<( Noncomparable const & ) const;
	};

	inline bool
	Noncomparable::operator<( Noncomparable const & ) const
	{
		assert( !"Noncomparable::operator<" );

		return true;
	}
}

#endif // COMMON_NONCOMPARABLE_H
