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

#ifndef HEX_MISC
#define HEX_MISC

#include <exception>
#include <cstdio>

#include <limits>

/// Infinity
#define Inf (std::numeric_limits<double>::infinity())
/// Not-a-number
#define Nan (std::numeric_limits<double>::quiet_NaN())

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

/// Trim from start.
static inline std::string &ltrim(std::string &s)
{
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

/// Trim from end.
static inline std::string &rtrim(std::string &s)
{
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

/// Trim from both ends.
static inline std::string &trim(std::string &s)
{
	return ltrim(rtrim(s));
}

/**
 * Custom exception class with easy printf-like constructor.
 * 
 * Use something like:
 * \code
 * throw exception("[Error %d] Pointed has the value 0x%x!", id, ptr);
 * \endcode
 */
class exception : public std::exception
{
public:
	
	/// Constructor.
	template <class ...Params> exception (Params ...p)
	{
		std::snprintf(message, 256, p...);
	}
	
	/// Return pointer to the exception text.
	const char* what() const noexcept (true)
	{
		return message;
	}
	
private:
	
	/// Text of the exception.
	char message[256];
};

#endif
