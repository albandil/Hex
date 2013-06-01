/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cctype>
#include <iostream>
#include <sstream>

#include "misc.h"
#include "vec3d.h"

std::ostream & operator << (std::ostream & os, vec3d const & v)
{
	os << "(" << v.x << " " << v.y << " " << v.z << ")";
	return os;
}

std::istream & operator >> (std::istream & is, vec3d & v)
{
	// allowed characters in a number
	std::string numc = ".0123456789e+-";
	
	// read buffer
	char c;
	
	// component strings
	std::string xs, ys, zs;
	
	// read components
	while (true)
	{
		// read a character
		is >> c;
		
		// skip left parenthesis and quotation marks
		if (c == '(' or c == '"')
			continue;
		
		if (isspace(c))
		{
			// skip leading white spaces
			if (xs.empty())
				continue;
			
			// break on trailing white spaces
			else
				break;
		}
		
		// terminate on non-numerical characters
		if (numc.find(c) == std::string::npos)
			throw exception ("ERROR: Unexpected character '%c'.", c);
	}
	while (true)
	{
		// read a character
		is >> c;
		
		if (isspace(c))
		{
			// skip leading white spaces
			if (xs.empty())
				continue;
			
			// break on trailing white spaces
			else
				break;
		}
		
		// terminate on non-numerical characters
		if (numc.find(c) == std::string::npos)
			throw exception ("ERROR: Unexpected character '%c'.", c);
	}
	while (true)
	{
		// read a character
		is >> c;
		
		// break on right parenthesis or quotation marks
		if (c == ')' or c == '"')
			break;
		
		if (isspace(c))
		{
			// skip leading white spaces
			if (xs.empty())
				continue;
			
			// break on trailing white spaces
			else
				break;
		}
		
		// terminate on other non-numerical characters
		if (numc.find(c) == std::string::npos)
			throw exception ("ERROR: Unexpected character '%c'.", c);
	}
	
	// convert components
	std::istringstream ixs(xs), iys(ys), izs(zs);
	ixs >> v.x;
	iys >> v.y;
	izs >> v.z;
	
	return is;
}
