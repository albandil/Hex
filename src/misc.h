/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_MISC
#define HEX_MISC

#include <exception>
#include <chrono>
#include <complex>
#include <cstdio>
#include <iostream>
#include <limits>

/**
 * @brief Exception class.
 * 
 * Custom exception class with easy printf-like constructor.
 * 
 * Use something like:
 * @code
 * throw exception("[Error %d] Pointer has the value 0x%x!", id, ptr);
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

//
// Special numbers.
//

/// Infinity
#define Inf (std::numeric_limits<double>::infinity())
/// Not-a-number
#define Nan (std::numeric_limits<double>::quiet_NaN())

/// Pi
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/// 2/Pi
#ifndef M_2_PI
    #define M_2_PI 0.63661977236758134308
#endif

/// square root of 2
#ifndef M_SQRT2
    #define M_SQRT2 1.41421356237309504880
#endif

// # define M_E           2.7182818284590452354   /* e */
// # define M_LOG2E        1.4426950408889634074   /* log_2 e */
// # define M_LOG10E       0.43429448190325182765  /* log_10 e */
// # define M_LN2          0.69314718055994530942  /* log_e 2 */
// # define M_LN10         2.30258509299404568402  /* log_e 10 */
// # define M_PI           3.14159265358979323846  /* pi */
// # define M_PI_2         1.57079632679489661923  /* pi/2 */
// # define M_PI_4         0.78539816339744830962  /* pi/4 */
// # define M_1_PI         0.31830988618379067154  /* 1/pi */
// # define M_2_PI         0.63661977236758134308  /* 2/pi */
// # define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
// # define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
// # define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <type_traits>

//
// Restricted pointers.
//

#ifndef restrict
#ifdef __GNUC__
    #define restrict __restrict
#else
    #define restrict
    #warning "Don't know how to use restricted pointers with this compiler. The resulting code will be slow."
#endif
#endif

//
// Memory alignment.
// - supports intrinsic compiler optimization, mainly the usage of SIMD instructions (AVX etc.)
// - enabled only for Linux systems (those ought to posess "posix_memalign", which is used)
//

#ifndef __linux__
    #define NO_ALIGN
#endif

//
// List all complex types.
//

template <class T>
struct is_complex { static const bool value = false; };
template<> template <class T>
struct is_complex<std::complex<T>> { static const bool value = true; };

//
// List all scalar types variables.
//

#define declareTypeAsScalar(T)                                                \
                                                                              \
template <> struct is_scalar<T>                                               \
{                                                                             \
    static const bool value = true;                                           \
};

// declare general type as non-scalar (if not complex)
template <class T> struct is_scalar
{
    static const bool value = is_complex<T>::value;
};

// explicitly list all basic scalar types
declareTypeAsScalar(int)
declareTypeAsScalar(long)
declareTypeAsScalar(float)
declareTypeAsScalar(double)

//
// White space trimming.
//

/// Trim from start.
static inline std::string & ltrim (std::string & s)
{
    s.erase (
        s.begin(),
        std::find_if (
            s.begin(),
            s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))
        )
    );
    return s;
}

/// Trim from end.
static inline std::string & rtrim (std::string & s)
{
    s.erase (
        std::find_if (
            s.rbegin(),
            s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))
        ).base(), s.end()
    );
    return s;
}

/// Trim from both ends.
static inline std::string & trim (std::string & s)
{
    return ltrim(rtrim(s));
}

//
// Miscellaneous functions.
//

/// Signum function.
template <class T> int signum (T x)
{
    if (x < T(0))
        return -1;
    else if (x == T(0))
        return 0;
    else
        return 1;
}

/// Many-argument "min" function.
template <typename T> T mmin (T x)
{
    return x;
}
template <typename T, class ...Params> T mmin (T x, Params ...p)
{
    T y = mmin(p...);
    return std::min(x, y);
}

/// Many-argument "max" function.
template <typename T> T mmax (T x)
{
    return x;
}
template <typename T, class ...Params> T mmax (T x, Params ...p)
{
    T y = mmax(p...);
    return std::max(x, y);
}

/// Constant-expression max.
template <class T> constexpr T const & larger_of (T const & a, T const & b)
{
    return (a > b) ? a : b;
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
 * @brief Timing class.
 * 
 * The Timer class can be used for a comfortable computation of
 * elapsed time. The usage would be:
 * @code
 *     Timer timer;
 * 
 *     // .. block ...
 * 
 *     std:cout << "Time = " << timer.elapsed() << "secs.\n";
 * @endcode
 */
class Timer
{
    public:
        
        /// Return object reference (singleton interface).
        Timer()
            : start_(std::chrono::system_clock::now()) {}
        
        /// Start timer.
        void reset ()
        {
            start_ = std::chrono::system_clock::now();
        }
        
        /// Stop timer and return elapsed time in seconds.
        int elapsed ()
        {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::now(); // ? steady_clock
            std::chrono::seconds secs = std::chrono::duration_cast<std::chrono::seconds>(end - start_);
            return secs.count();
        }
        
    private:
        
        /// Start time.
        mutable std::chrono::system_clock::time_point start_;
};

#endif
