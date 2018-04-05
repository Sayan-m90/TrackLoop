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

#ifndef HEADERS_COLORED_H
#define HEADERS_COLORED_H

#include <CGAL/IO/Color.h>

namespace Headers
{
	class Colored
	{

	public:

		Colored();

		CGAL::Color const &color() const;

		void set_color( CGAL::Color const &color_ );

	private:

		CGAL::Color m_color;
	};

	inline
	Colored::Colored()
		: m_color( CGAL::WHITE )
	{
	}

	inline CGAL::Color const &
	Colored::color() const
	{
		return m_color;
	}

	inline void
	Colored::set_color( CGAL::Color const &color_ )
	{
		m_color = color_;
	}
}

#endif // COMMON_COLORED_H
