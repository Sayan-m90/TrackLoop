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

#ifndef HEADERS_FLAGGED_H
#define HEADERS_FLAGGED_H

namespace Headers
{
	class Flagged
	{

	public:

		Flagged();

		bool has_flag( unsigned flag_ ) const;
		void set_flag( unsigned flag_ );
		void clear_flag( unsigned flag_ );
		void toggle_flag( unsigned flag_ );

	private:

		mutable unsigned m_flags;
	};

	inline
	Flagged::Flagged()
		: m_flags( 0 )
	{
	}

	inline bool
	Flagged::has_flag( unsigned flag_ ) const
	{
		return ( m_flags & flag_ ) != 0;
	}

	inline void
	Flagged::set_flag( unsigned flag_ )
	{
		m_flags |= flag_;
	}

	inline void
	Flagged::clear_flag( unsigned flag_ )
	{
		m_flags &= ~flag_;
	}

	inline void
	Flagged::toggle_flag( unsigned flag_ )
	{
		m_flags ^= flag_;
	}
}

#endif // COMMON_FLAGGED_H
