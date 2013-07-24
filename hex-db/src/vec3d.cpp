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
	os << "( " << v.x << " " << v.y << " " << v.z << " )";
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
		is >> std::noskipws >> c;
		
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
		
		// add character to string
		xs.push_back(c);
	}
	while (true)
	{
		// read a character
		is >> std::noskipws >> c;
		
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
		
		// add character to string
		ys.push_back(c);
	}
	while (true)
	{
		// read a character
		is >> std::noskipws >> c;
		
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
		
		// add character to string
		zs.push_back(c);
	}
	
	// convert components
	std::istringstream ixs(xs), iys(ys), izs(zs);
	ixs >> v.x;
	iys >> v.y;
	izs >> v.z;
	
	return is;
}

double dot (vec3d const & u, vec3d const & v)
{
	return u.x*v.x + u.y*v.y + u.z*v.z;
}

vec3d cross (vec3d const & u, vec3d const & v)
{
	return vec3d (
		{
			u.y * v.z - u.z * v.y,
			u.z * v.x - u.x * v.z,
			u.x * v.y - u.y * v.x
		}
	);
}

vec3d operator - (vec3d const & u, vec3d const & v)
{
	return vec3d (
		{
			u.x - v.x,
			u.y - v.y,
			u.z - v.z
		}
	);
}

vec3d operator * (vec3d const & u, double a)
{
	return vec3d (
		{
			a * u.x,
			a * u.y,
			a * u.z
		}
	);
}

double norm (vec3d const & v)
{
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

vec3d normalize(vec3d const & v)
{
	return v * (1/norm(v));
}
