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

// restricted pointers
#ifndef restrict
#ifdef __GNUC__
    #define restrict __restrict
#else
    #define restrict
    #warning "Don't know how to use restricted pointers with this compiler. The resulting code will be slow."
#endif
#endif

// memory alignment
#ifndef alignof
#ifdef __GNUC__
    #define alignof(x)   (__alignof(x))
    #define aligned(x,y) (__builtin_assume_aligned((x),(y)))
#else
    #define alignof(x)   (sizeof(void*))
    #define aligned(x,y) (x)
    #warning "Don't know how to determine memory alignment. Using non-aligned pointers (may forbid vectorization and result in slower code)."
#endif
#endif

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

/// Signum function.
template <class T> int signum (T x)
{
    if (x < T(0)) return -1;
    else if (x == T(0)) return 0;
    else return 1;
}

/// Many-argument "min".
template <typename T> T min (T x)
{
    return x;
}
template <typename T, class ...Params> T min (T x, Params ...p)
{
    T y = min(p...);
    return std::min(x, y);
}

/**
 * @brief printf-like formatting.
 * 
 * @note Hard limit 1024 characters.
 */
template <class ...Params> char const * format (Params ...p)
{
    static char text[1024];
    snprintf(text, sizeof(text), p...);
    return text;
}

/**
 * @brief Exception class.
 * 
 * Custom exception class with easy printf-like constructor.
 * 
 * Use something like:
 * @code
 * throw exception("[Error %d] Pointed has the value 0x%x!", id, ptr);
 * @endcode
 * 
 * @note Hard limit 256 charasters.
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
    const char* what() const noexcept (true)
    {
        return message;
    }
    
private:
    
    /// Text of the exception.
    char message[256];
};

#endif
