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

//
// Special numbers.
//

/// Infinity
#define Inf (std::numeric_limits<double>::infinity())
/// Not-a-number
#define Nan (std::numeric_limits<double>::quiet_NaN())

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <type_traits>

//
// Restricted and aligned pointers.
//

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
    #define aligned_ptr(x,y) (__builtin_assume_aligned((x),(y)))
#else
    #define aligned_ptr(x,y) (x)
    #warning "Don't know how to determine memory alignment. Using non-aligned pointers (may forbid vectorization and result in slower code)."
#endif
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
declareTypeAsScalar(int);
declareTypeAsScalar(long);
declareTypeAsScalar(float);
declareTypeAsScalar(double); 

//
// Trimming.
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
 * The Timer class is a singleton that can be used for a comfortable computation of
 * elapsed time. The usage would be:
 * @code
 *     Timer::timer().start();
 * 
 *     // .. block ...
 * 
 *     std:cout << "Time = " << Timer::timer().stop() << "secs.\n";
 * @endcode
 */
class Timer
{
    public:
        
        /// Return object reference (singleton interface).
        static Timer const & timer()
        {
            static Timer timer_;
            return timer_;
        }
        
        /// Start timer.
        void start () const
        {
            start_ = std::chrono::system_clock::now();
        }
        
        /// Stop timer and return elapsed time in seconds.
        int stop () const
        {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
            std::chrono::seconds secs = std::chrono::duration_cast<std::chrono::seconds>(end - start_);
            return secs.count();
        }
        
    private:
        
        /// Start time.
        mutable std::chrono::system_clock::time_point start_;
};

#endif
