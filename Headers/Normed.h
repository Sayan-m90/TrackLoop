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

#ifndef HEADERS_NORMED_H
#define HEADERS_NORMED_H

namespace Headers
{	
	class Normed
	{

	public:

		Normed();

		double norm() const;
		void set_norm( double norm_ );

	private:

		double m_norm;
	};

	inline 
	Normed::Normed()
		: m_norm( 0 )
	{
	}

	inline double
	Normed::norm() const
	{
		return m_norm;
	}

	inline void
	Normed::set_norm( double norm_ )
	{
		m_norm = norm_;
	}

	struct Pointee_norm_is_less
	{
		bool operator()( Normed const *p_a_, Normed const *p_b_ ) const
		{
			return p_a_->norm() < p_b_->norm();
		}
	};

	inline bool
	pointee_norm_is_less( Normed const *p_a_, Normed const *p_b_ )
	{
		return p_a_->norm() < p_b_->norm();
	}

	struct Pointee_norm_is_greater
	{
		bool operator()( Normed const *p_a_, Normed const *p_b_ ) const
		{
			return p_a_->norm() > p_b_->norm();
		}
	};

	inline bool
	pointee_norm_is_greater( Normed const *p_a_, Normed const *p_b_ )
	{
		return p_a_->norm() > p_b_->norm();
	}
}

#endif // COMMON_NORMED_H
