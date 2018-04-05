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

#ifndef HEADERS_CHAIN_H
#define HEADERS_CHAIN_H

#include <Noncomparable.h>

#include <list>

namespace Headers
{	
	template< typename Value_ >
	class Chain
		: public Noncomparable
	{

	protected:

		enum
		{
			IS_IN_CHAIN = 1 << 30
		};

	public:

		typedef std::list< Value_ * > Values;

		typedef typename Values::const_iterator Const_iterator;
		typedef typename Values::iterator Iterator;
		typedef typename Values::const_reverse_iterator Const_reverse_iterator;
		typedef typename Values::reverse_iterator Reverse_iterator;
	
	public:	

		bool empty() const;
		unsigned size() const;

		Value_ const &front() const;
		Value_ &front();

		Value_ const &back() const;
		Value_ &back();

		Const_iterator begin() const;
		Const_iterator end() const;
		Iterator begin();
		Iterator end();
		Const_reverse_iterator rbegin() const;
		Const_reverse_iterator rend() const;
		Reverse_iterator rbegin();
		Reverse_iterator rend();

		void insert( Value_ &value_ );
		void erase( Iterator const &it_ );
		void clear();

		void operator^=( Chain< Value_ > &chain_ );

		std::string name() const;

	private:

		mutable Values m_values;
	};	

	template< typename Value_ >
	inline bool
	Chain< Value_ >::empty() const
	{
		return m_values.empty();
	}

	template< typename Value_ >
	inline unsigned
	Chain< Value_ >::size() const
	{
		return m_values.size();
	}	

	
	template< typename Value_ >
	inline Value_ const &
	Chain< Value_ >::front() const
	{
		return *m_values.front();
	}
	
	template< typename Value_ >
	inline Value_ &
	Chain< Value_ >::front()
	{
		return *m_values.front();
	}
	
	template< typename Value_ >
	inline Value_ const &
	Chain< Value_ >::back() const
	{
		return *m_values.back();
	}

	template< typename Value_ >
	inline Value_ &
	Chain< Value_ >::back()
	{
		return *m_values.back();
	}
	
	template< typename Value_ >
	inline typename Chain< Value_ >::Const_iterator
	Chain< Value_ >::begin() const
	{
		return m_values.begin();
	}

	template< typename Value_ >
	inline typename Chain< Value_ >::Const_iterator
	Chain< Value_ >::end() const
	{
		return m_values.end();
	}

	template< typename Value_ >
	inline typename Chain< Value_ >::Iterator
	Chain< Value_ >::begin()
	{
		return m_values.begin();
	}

	template< typename Value_ >
	inline typename Chain< Value_ >::Iterator
	Chain< Value_ >::end()
	{
		return m_values.end();
	}	

	template< typename Value_ >
	inline typename Chain< Value_ >::Const_reverse_iterator
	Chain< Value_ >::rbegin() const
	{
		return m_values.rbegin();
	}

	template< typename Value_ >
	inline typename Chain< Value_ >::Const_reverse_iterator
	Chain< Value_ >::rend() const
	{
		return m_values.rend();
	}

	template< typename Value_ >
	inline typename Chain< Value_ >::Reverse_iterator
	Chain< Value_ >::rbegin()
	{
		return m_values.rbegin();
	}

	template< typename Value_ >
	inline typename Chain< Value_ >::Reverse_iterator
	Chain< Value_ >::rend()
	{
		return m_values.rend();
	}	

	template< typename Value_ >
	inline void
	Chain< Value_ >::insert( Value_ &value_ )
	{		
		m_values.push_back( &value_ );
	}

	template< typename Value_ >
	inline void
	Chain< Value_ >::erase( Iterator const &it_ )
	{
		m_values.erase( it_ );
	}

	template< typename Value_ >
	inline void
	Chain< Value_ >::clear()
	{
		m_values.clear();
	}	

	template< typename Value_ >
	inline void
	Chain< Value_ >::operator^=( Chain< Value_ > &chain_ )
	{		
		typename Values::iterator it( chain_.m_values.begin() );
		for ( ; it != chain_.m_values.end(); ++it )
			( *it )->set_flag( IS_IN_CHAIN );

		it = m_values.begin();
		while ( it != m_values.end() )
		{
			typename Values::iterator it_next( it );
			++it_next;

			if ( ( *it )->has_flag( IS_IN_CHAIN ) )
			{
				( *it )->clear_flag( IS_IN_CHAIN );
				erase( it );
			}
			
			it = it_next;
		}

		it = chain_.m_values.begin();
		for ( ; it != chain_.m_values.end(); ++it )
		{
			if ( ( *it )->has_flag( IS_IN_CHAIN ) )
			{
				( *it )->clear_flag( IS_IN_CHAIN );
				insert( **it );
			}
		}
	}

	template< typename Value_ >
	inline std::string
	Chain< Value_ >::name() const
	{
		std::string result;
		for ( Const_iterator it( begin() ); it != end(); ++it )
		{
			if ( it != begin() )
				result.push_back( ' ' );
			result.append( ( **it ).name() );
		}

		return result;
	}
}

#endif // COMMON_CHAIN_H
