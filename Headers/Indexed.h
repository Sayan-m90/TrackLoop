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

#ifndef HEADERS_INDEXED_H
#define HEADERS_INDEXED_H

namespace Headers
{	
	class Indexed
	{

	public:

		Indexed();

		unsigned index() const;
		void set_index( unsigned index_ );

	private:

		unsigned m_index;
	};

	inline 
	Indexed::Indexed()
		: m_index( 0 )
	{
	}

	inline unsigned
	Indexed::index() const
	{
		return m_index;
	}

	inline void
	Indexed::set_index( unsigned index_ )
	{
		m_index = index_;
	}

	struct Pointee_index_is_less
	{
		bool operator()( Indexed const *p_a_, Indexed const *p_b_ ) const
		{
			return p_a_->index() < p_b_->index();
		}
	};
	
	inline bool
	pointee_index_is_less( Indexed const *p_a_, Indexed const *p_b_ )
	{
		return p_a_->index() < p_b_->index();
	}

	struct Pointee_index_is_greater
	{
		bool operator()( Indexed const *p_a_, Indexed const *p_b_ ) const
		{
			return p_a_->index() > p_b_->index();
		}
	};

	inline bool
	pointee_index_is_greater( Indexed const *p_a_, Indexed const *p_b_ )
	{
		return p_a_->index() > p_b_->index();
	}
}

#endif // COMMON_INDEXED_H
