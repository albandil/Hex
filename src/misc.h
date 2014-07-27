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

#include <algorithm> 
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <complex>
#include <exception>
#include <fstream>
#include <functional>
#include <ios>
#include <iomanip>
#include <limits>
#include <locale>
#include <iostream>
#include <type_traits>
#include <string>

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
 * @note The behaviour of the function "snprintf" that is used in the constructor
 *       is expected to follow C99 specification, i.e. that it returns the needed
 *       size of the output buffer when called as snprinf(nullptr,0,format,_args_).
 */
class exception : public std::exception
{
public:
    
    /// Constructor.
    template <class ...Params> exception (Params ...p)
    {
        // get requested size
        int size = snprintf(nullptr, 0, p...);
        
        // allocate buffer
        message = new char [size + 1];
        
        // printf to the allocated storage
        snprintf(message, size + 1, p...);
    }
    
    /// Destructor.
    virtual ~exception() noexcept
    {
        delete [] message;
    }
    
    /// Return pointer to the exception text.
    const char* what() const noexcept (true)
    {
        return message;
    }
    
private:
    
    /// Text of the exception.
    char * message;
};

//
// Restricted pointers.
// - allow some compiler optimizations on arrays
//

#ifndef restrict
#ifdef __GNUC__
    #define restrict __restrict
#else
    #define restrict
    #warning "Don't know how to use restricted pointers with this compiler. The resulting code may be slower."
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
// - used in the list of scalar types below
// - by default, a type is NOT complex
// - only std::complex<T> is an exception from that rule
//

template <class T> struct is_complex { static const bool value = false; };

#define declareTypeAsComplex(T) \
    template <> struct is_complex<T> { static const bool value = true; }
#define declareTemplateTypeAsComplex(TT) \
    template <> template <class T> struct is_complex<TT<T>> { static const bool value = true; }

declareTemplateTypeAsComplex(std::complex);

//
// List all scalar types.
// - necessary for blocking some scalar array functions for non-scalar types (see arrays.h)
// - by default, a type is not scalar if it is not complex
// - explicitly list all other basic types considered "scalar"
//

template <class T> struct is_scalar { static const bool value = is_complex<T>::value; };

#define declareTypeAsScalar(T) \
    template <> struct is_scalar<T> { static const bool value = true; }

declareTypeAsScalar(int);
declareTypeAsScalar(long);
declareTypeAsScalar(float);
declareTypeAsScalar(double);

//
// White space trimming.
//

/**
 * @brief Trim from start.
 * 
 * Removes all leading white-space characters. The classification is
 * done by the function std::isspace.
 */
static inline std::string & ltrim (std::string & s)
{
    s.erase
    (
        s.begin(),
        std::find_if
        (
            s.begin(),
            s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))
        )
    );
    return s;
}

/**
 * @brief Trim from end.
 * 
 * Removes all trailing white-space characters. The classification is
 * done by the function std::isspace.
 */
static inline std::string & rtrim (std::string & s)
{
    s.erase
    (
        std::find_if
        (
            s.rbegin(),
            s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))
        ).base(), s.end()
    );
    return s;
}

/**
 * @brief Trim both ends.
 * 
 * Removes white-space characters from the beginning and end.
 * The routine calls @ref ltrim and @ref rtrim.
 */
static inline std::string & trim (std::string & s)
{
    return ltrim(rtrim(s));
}

//
// Miscellaneous functions.
//

/**
 * @brief Signum function.
 * 
 * This function returns the sign of the argument or zero if the argument
 * is zero.
 */
template <class T> int signum (T x)
{
    if (x < T(0))
        return -1;
    else if (x == T(0))
        return 0;
    else
        return 1;
}

/**
 * @brief Many-argument "min" function.
 * 
 * This function returns the smallest of the arguments. There can be
 * arbitrary number of arguments. The comparison is done by std::min
 * function.
 */
//@{
template <typename T> T mmin (T x)
{
    return x;
}
template <typename T, class ...Params> T mmin (T x, Params ...p)
{
    T y = mmin(p...);
    return std::min(x, y);
}
//@}

/**
 * @brief Many-argument "max" function.
 * 
 * This function returns the largest of the arguments. There can be
 * arbitrary number of arguments. The comparison is done by std::max
 * function.
 */
//@{
template <typename T> T mmax (T x)
{
    return x;
}
template <typename T, class ...Params> T mmax (T x, Params ...p)
{
    T y = mmax(p...);
    return std::max(x, y);
}
//@}

/**
 * @brief Constant-expression max.
 * 
 * This is a plain comparison
   @code
   (a > b) ? a : b
   @endcode
 * However, it is marked as "constexpr" to allow its evaluation at compile
 * time. This is used somewhere in the code.
 */
template <class T> constexpr T const & larger_of (T const & a, T const & b)
{
    return (a > b) ? a : b;
}

/**
 * @brief printf-like formatting.
 * 
 * This function takes an arbitrary number of parameters. It is expected that
 * the first one is the formatting string (printf-like syntax). All the
 * arguments are sent to snprintf without change. This functions returns
 * a pointer to a static character string.
 * 
 * @note The maximal size of the string is hard-coded to 1024 characters.
 */
template <class ...Params> char const * format (Params ...p)
{
    static char text[1024];
    snprintf(text, sizeof(text), p...);
    return text;
}

/**
 * @brief Conversion of string to a type.
 * 
 * This is a generic template function that is used to convert a text entry to
 * the specified data type. This generic form is used only if there is no specialization
 * for the given data type. The default action is an error message, then.
 * See the specializations of this function for different types.
 */
template <class T> T string_to (std::string str)
{
    throw exception
    (
        "Conversion to \"%s\" not implemented.", typeid(T).name()
    );
}

/**
 * @brief Conversion of a string to integer number.
 * 
 * This function will return an integer value of the text given as argument.
 * The library routine "strtol" is used. If the conversion fails, the
 * function throws an exception.
 */
template <> inline int string_to (std::string str)
{
    // convert to int
    char* tail; long val = std::strtol(str.c_str(), &tail, 10);
    
    // throw or return
    if (*tail != 0x0)
        throw exception ("The string \"%s\" cannot be converted to integer number.", str.c_str());
    else
        return val;
}

/**
 * @brief Conversion of a string to floating-point number.
 * 
 * This function will return a floating-point value of the text given as argument.
 * The library routine "strtod" is used. If the conversion fails, the
 * function throws an exception.
 */
template <> inline double string_to (std::string str)
{
    // convert to float
    char* tail; double val = std::strtod(str.c_str(), &tail);
    
    // throw or return
    if (*tail != 0x0)
        throw exception ("The string \"%s\" cannot be converted to real number.", str.c_str());
    else
        return val;
}

/**
 * @brief Output structure for the function @ref ReadNext.
 * 
 * This structure is returned by the function @ref ReadNext and it contains
 * the read value and / or some special values. Currently, only the special character '*'
 * is handled -- in that case 
 */
template <class T> struct ReadItem
{
    enum Flags
    {
        none     = 0x00,
        asterisk = 0x01
    };
    
    T val;
    int flags;
};

/**
 * @brief Read next entry from input stream.
 * 
 * Given an input stream and a template parameter the function "read_next" will
 * scan the stream for the next entry and try to interpret next input as the correct type.
 * The characters between a hash symbol (#) and a newline are ignored (i.e. '#' introduces
 * comments). If the read word is equal to asterisk, the "true" boolean value is trown.
 * If an error occurs, @ref exception is thrown. Otherwise the converted data is returned.
 * 
   @code
       // example usage
       std::ifstream inputfile("input.txt");
       double x = read_next<double>(inputfile);
   @endcode
 */
template <class T> ReadItem<T> ReadNext (std::ifstream & f, unsigned allowed_special = ReadItem<T>::none)
{
    // text buffer
    std::string s;
    
    // while there is something in the file
    while (not f.eof())
    {
        // read string (it won't start with white character)
        f >> s;
        
        // check length (skip empty reads)
        if (s.size() == 0)
            continue;
        
        // check if it is a beginning of a comment
        if (s.front() == '#')
        {
            // get the rest of the line
            std::getline(f, s);
            continue;
        }
        
        // otherwise exit the loop (a valid entry was found)
        break;
    }
    
    // check for asterisk
    if (s == "*")
    {
        if (allowed_special & ReadItem<T>::asterisk)
            return ReadItem<T>({ T(0), ReadItem<T>::asterisk });
        else
            throw exception ("Asterisk '*' not allowed here.");
    }
    
    // convert entry to type T
    return ReadItem<T>({ string_to<T>(s.c_str()), ReadItem<T>::none });
}

/**
 * @brief Timing class.
 * 
 * The Timer class can be used for a comfortable computation of
 * elapsed time. The usage would be:
   @code
       Timer timer;
   
       // .. block ...
   
       std:cout << "Time = " << timer.seconds() << "secs.\n";
   @endcode
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
        
        /// Return elapsed time in seconds.
        unsigned seconds ()
        {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::now(); // ? steady_clock
            std::chrono::seconds secs = std::chrono::duration_cast<std::chrono::seconds>(end - start_);
            return secs.count();
        }
        
        /// Return elapsed time in milliseconds.
        unsigned milliseconds ()
        {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::now(); // ? steady_clock
            std::chrono::milliseconds misecs = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_);
            return misecs.count();
        }
        
        /// Return elapsed time in microseconds.
        unsigned microseconds ()
        {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::now(); // ? steady_clock
            std::chrono::microseconds musecs = std::chrono::duration_cast<std::chrono::microseconds>(end - start_);
            return musecs.count();
        }
        
    private:
        
        /// Start time.
        mutable std::chrono::system_clock::time_point start_;
};

/**
 * @brief Output table.
 * 
 * This class enables table-formatted output. An example code is
 * @code
 * std::ofstream file ("output.txt");
 * OutputTable table (file);
 * table.setWidth (5, 5);
 * table.setAlign (OutputTable::left, OutputTable::right);
 * table.write ("abc", "xyz");
 * table.write (10, 17.5);
 * table.write (" x", "xx");
 * @endcode
 * The resulting output will be
 * @code
 * abc    xyz
 * 10    17.5
 *  x      xx
 * @endcode
 */
class OutputTable
{
    public:
        
        OutputTable (std::ostream & out = std::cout)
            : width_(1, 15), alignment_(1, OutputTable::left), out_(out.rdbuf()) {}
        
        typedef enum
        {
            none = 0,
            left = 1,
            center = 2, // not implemented yet
            right = 3
        } align;
        
        void setStream (std::ostream & out)
        {
            out_.rdbuf(out.rdbuf());
        }
        
        template <class ...Params> void setWidth (Params ...p)
        {
            // clear old widths
            width_.clear();
            
            // set new widths
            append_(width_, p...);
        }
        
        template <class ...Params> void setAlignment (Params ...p)
        {
            // clear old alignments
            alignment_.clear();
            
            // set new alignments
            append_(alignment_, p...);
        }
        
        template <int i> void write ()
        {
            out_ << std::endl;
        }
        
        template <int i = 0, class T, class ...Params> void write (T item, Params ...p)
        {
            // use field width and alignment
            int width = (i < width_.size() ? width_[i] : width_.back());
            int alignment = (i < alignment_.size() ? alignment_[i] : alignment_.back());
            
            out_ << std::setw(width);
            
            switch (alignment)
            {
                case none:
                    break;
                case left:
                    out_ << std::left;
                    break;
                case center:
                    throw exception ("[OutputTable] Item centering not implemented yet!");
                    break;
                case right:
                    out_ << std::right;
                    break;
            }
            
            // write the item
            out_ << item;
            
            // write the rest of the items
            write<i+1>(p...);
        }
        
    private:
        
        std::vector<int> width_;
        std::vector<align> alignment_;
        
        std::ostream out_;
        
        template <class T> void append_ (std::vector<T> & v)
        {
            // do nothing
        }
        
        template <class T, class ...Params> void append_ (std::vector<T> & v, T w, Params ...p)
        {
            v.push_back(w);
            append_(v, p...);
        }
};

/**
 * @brief Range class.
 * 
 * This class holds a pair of values - the first and last one. It can be initialized
 * by a string in the form "&lt;first&gt;-&lt;last&gt;".
 */
template <class T> class Range
{
    public:
        
        // the first element of the range
        T first;
        
        // the last element of the range
        T last;
        
        // constructor
        Range (std::string Lin)
        {
            // find separator ('-')
            std::string::const_iterator sep = std::find(Lin.begin(),Lin.end(),'-');
            
            // check if the separator was found
            if (sep != Lin.end())
            {
                // split the string at the separator
                std::string strfirst = Lin.substr(0,sep-Lin.end());
                std::string strlast  = Lin.substr(sep-Lin.begin()+1,Lin.end()-sep-1);
                
                // convert portions to the correct type
                std::istringstream infirst(strfirst), inlast(strlast);
                if (!(infirst >> first)) throw exception ("\"%s\" is not a number.", strfirst.c_str());
                if (!(inlast >> last)) throw exception ("\"%s\" is not a number.", strlast.c_str());
            }
            else
            {
                // convert whole string to the correct type
                std::istringstream in(Lin);
                if (!(in >> first)) throw exception ("\"%s\" is not a number.", Lin.c_str());
                
                // the first and last element is the same
                last = first;
            }
        }
};

/**
 * @brief Underline text.
 * 
 * This simple function will count characters in the supplied text string,
 * append a new-line character and a group of '-' characters whose number
 * will be equal to the length of the supplied text. It is UTF-8-aware, so
 * the following code
   @code
       std::cout << "Title containing UTF-8 characters αβγδ" << std::endl;
   @endcode
 * will produce the output
   @verbatim
       Title containing UTF-8 characters αβγδ
       --------------------------------------
   @endverbatim
 * It may be necessary to set the user locale, first, using the function call
   @code
       std::setlocale(LC_ALL, "");
   @endcode
 */
inline std::string underline (std::string text)
{
    std::mbstate_t state = std::mbstate_t();
    const char * mbstr = text.data();
    
    std::size_t wlen = std::mbsrtowcs(NULL, &mbstr, 0, &state);
    
    text.push_back('\n');
    for (std::size_t i = 0; i < wlen; i++)
        text.push_back('-');
    return text;
}

#endif
