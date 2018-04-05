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

#ifndef HEADERS_FILE_H
#define HEADERS_FILE_H

#include <Polymorphic.h>
#include <boost/filesystem/path.hpp>

namespace Headers
{
	class File
		: public Polymorphic
	{

	public:

		File( char const *path_cstring_ );

		boost::filesystem::path const &path() const;

	private:

		boost::filesystem::path m_path;
	};

	inline
	File::File( char const *path_cstring_ )
		: m_path( path_cstring_ )
	{
	}		

	inline
	boost::filesystem::path const &
	File::path() const
	{
		return m_path;
	}
}

#endif // COMMON_FILE_H
