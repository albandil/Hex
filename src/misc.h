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

/**
 * Custom exception class with easy printf-like constructor.
 */
class exception : public std::exception
{
public:
	
	/// Constructor.
	template <class ...Params> exception (Params ...p)
	{
		snprintf(message, 256, p...);
	}
	
	/// Return pointer to the exception text.
	const char* what() const noexcept (true) {
		return message;
	}
	
private:
	
	/// Text of the exception.
	char message[256];
};

#endif
