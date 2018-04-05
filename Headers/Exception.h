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

#ifndef HEADERS_EXCEPTION_H
#define HEADERS_EXCEPTION_H

#include <exception>
#include <string>

namespace Headers
{
	class Exception
		: public std::exception
	{

	public:

		Exception() throw();
		Exception( char const *what_ ) throw();
		Exception( Exception const &e_ ) throw();
		Exception &operator=( Exception const &e_ ) throw();
		/* virtual */ ~Exception() throw();  
		
		/* virtual */ char const *what() const throw();

	public:

		std::string m_what;
	};

	inline
	Exception::Exception() throw()
	{
	}

	inline
	Exception::Exception( char const *what_ ) throw()
		: m_what( what_ )
	{
	}

	inline
	Exception::Exception( Exception const &e_ ) throw()
		: m_what( e_.what() )
	{
	}

	inline
	Exception &
	Exception::operator=( Exception const &e_ ) throw()
	{
		m_what = e_.what();
	}

	inline
	Exception::~Exception() throw()
	{
	}

	inline char const *
	Exception::what() const throw()
	{
		return m_what.c_str();
	}
}

#endif // COMMON_EXCEPTION_H
