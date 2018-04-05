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

#ifndef HEADERS_TEXT_OUTPUT_FILE_H
#define HEADERS_TEXT_OUTPUT_FILE_H

#include <File.h>
#include <Exception.h>

#include <cstdio>
#include <string>

namespace Headers
{
	class Text_output_file
		: public File
	{

	public:

		Text_output_file( char const *path_cstring_ );
		~Text_output_file();

	protected:

		FILE *m_p_file;
	};

	inline
	Text_output_file::Text_output_file( char const *path_cstring_ )
		: File( path_cstring_ ), m_p_file( fopen( path_cstring_, "w" ) )
	{
		if ( m_p_file == 0 )
		{
			throw Exception( ( std::string( "Cannot open file " )
				+ path_cstring_ ).c_str() );
		}
	}

	inline
	Text_output_file::~Text_output_file()
	{
		fclose( m_p_file );
	}
}

#endif // COMMON_TEXT_OUTPUT_FILE_H
