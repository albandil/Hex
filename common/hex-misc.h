//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#ifndef HEX_MISC
#define HEX_MISC

#include <algorithm> 
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
#include <vector>

#ifdef __linux__
    #include <unistd.h>
#endif

/**
 * @brief printf-like formatting.
 * 
 * This function takes an arbitrary number of parameters. It is expected that
 * the first one is the formatting string (printf-like syntax). All the
 * arguments are sent to snprintf without change.
 */
template <class ...Params> std::string format (Params ...p)
{
    // calculate the necessary space
    unsigned n = std::snprintf(nullptr, 0, p...);
    
    // allocate the string
    std::string text;
    text.resize(n);
    
    // compose the text
    std::snprintf(&text[0], n + 1, p...);
    
    // return the text
    return text;
}

/// Debug output.
#define Debug std::cout << __func__ << " (" << __FILE__ << ":" << __LINE__ << "): "

/// Fatal error routine.
#define HexException(...) TerminateWithException(__FILE__, __LINE__, __func__, __VA_ARGS__)

/// Prints backtrace. Implemented only for GNU/Linux.
void print_stack_trace ();

/**
 * @brief Fatal error routine.
 * 
 * Fatal error termination routine with easy printf-like interface.
 * Use the macro @ref HexException to call this function with proper line number
 * and function name.
 */
template <class ...Params> [[noreturn]] void TerminateWithException (const char* file, int line, const char* func, Params ...p)
{
    // print error text
    std::cerr << std::endl << std::endl;
    std::cerr << "Program unsuccessfully terminated (in " << file << ":" << line << ", function \"" << func << "\")" << std::endl;
    
    // POSIX terminal allows colours
#if (_POSIX_C_SOURCE >= 200112L)
    if (isatty(fileno(stderr)))
        std::cerr << "\033[1;31m *** " << format(p...) << "\033[0m" << std::endl;
    else
#endif
        std::cerr << " *** " << format(p...) << std::endl;
    
    // print stack trace (if available)
    print_stack_trace();
    
    // exit the program
    std::terminate();
}

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
        : message_(format(p...)) { }
    
    /// Return pointer to the exception text.
    const char* what () const noexcept (true)
    {
        return message_.c_str();
    }
    
private:
    
    /// Text of the exception.
    std::string message_;
};

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
    template <class T> struct is_complex<TT<T>> { static const bool value = true; }

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
declareTypeAsScalar(std::int64_t);
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
 * @brief Conversion of string to a type.
 * 
 * This is a generic template function that is used to convert a text entry to
 * the specified data type. This generic form is used only if there is no specialization
 * for the given data type. The default action is an error message, then.
 * See the specializations of this function for different types.
 */
template <class T> T string_to (std::string str)
{
    HexException("Conversion of string to \"%s\" not implemented.", typeid(T).name());
}

/**
 * @brief Conversion of a string to string.
 */
template <> inline std::string string_to (std::string str)
{
    return str;
}

/**
 * @brief Conversion of a string to character.
 * 
 * This function fails if length of the string is different from one.
 */
template <> inline char string_to (std::string str)
{
    // check length
    if (str.size() != 1)
        HexException("The string \"%s\" cannot be converted to character.", str.c_str());
    else
        return str[0];
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
    char* tail; int val = std::strtol(str.c_str(), &tail, 10);
    
    // throw or return
    if (*tail != 0x0)
        HexException("The string \"%s\" cannot be converted to integer number.", str.c_str());
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
        HexException("The string \"%s\" cannot be converted to real number.", str.c_str());
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
            HexException("Asterisk '*' not allowed here.");
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
        
        /// Constructor (with optional offset in seconds).
        Timer (int offset = 0)
            : offset_(offset), start_(std::chrono::system_clock::now()) {}
        
        /// Start timer.
        void reset ()
        {
            start_ = std::chrono::system_clock::now();
            offset_ = 0;
        }
        
        /// Return elapsed time in seconds.
        unsigned seconds () const
        {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::now(); // ? steady_clock
            std::chrono::seconds secs = std::chrono::duration_cast<std::chrono::seconds>(end - start_);
            return secs.count() + offset_;
        }
        
        /// Return formatted time.
        std::string nice_time () const
        {
            // get elapsed time
            unsigned secs = seconds();
            
            // return formatted time
            return format("%02d:%02d:%02d", secs / 3600, (secs % 3600) / 60, secs % 60);
        }
        
        /// Return elapsed time in milliseconds.
        unsigned milliseconds () const
        {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::now(); // ? steady_clock
            std::chrono::milliseconds misecs = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_);
            return misecs.count() + offset_ * 1000;
        }
        
        /// Return elapsed time in microseconds.
        unsigned microseconds () const
        {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::now(); // ? steady_clock
            std::chrono::microseconds musecs = std::chrono::duration_cast<std::chrono::microseconds>(end - start_);
            return musecs.count() + offset_ * 1000000;
        }
        
    private:
        
        /// Time offset (number of seconds to be added to the time).
        unsigned offset_;
        
        /// Start time.
        mutable std::chrono::system_clock::time_point start_;
};

/**
   @brief Output table.
   
   This class enables table-formatted output. An example code is
   @code
   std::ofstream file ("output.txt");
   OutputTable table (file);
   table.setWidth (5, 5);
   table.setAlign (OutputTable::left, OutputTable::right);
   table.write ("abc", "xyz");
   table.write (10, 17.5);
   table.write (" x", "xx");
   @endcode
   The resulting output will be
   @code
   abc    xyz
   10    17.5
    x      xx
   @endcode
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
                    HexException("[OutputTable] Item centering not implemented yet!");
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
                if (!(infirst >> first)) HexException("\"%s\" is not a number.", strfirst.c_str());
                if (!(inlast >> last)) HexException("\"%s\" is not a number.", strlast.c_str());
            }
            else
            {
                // convert whole string to the correct type
                std::istringstream in(Lin);
                if (!(in >> first)) HexException("\"%s\" is not a number.", Lin.c_str());
                
                // the first and last element is the same
                last = first;
            }
        }
};

/**
   @brief Underline text.
   
   This simple function will count characters in the supplied text string,
   append a new-line character and a group of '-' characters whose number
   will be equal to the length of the supplied text. It supports UTF-8, so
   the following code
   @code
       std::cout << underline("Title containing UTF-8 characters αβγδ") << std::endl;
   @endcode
   will produce the output
   @verbatim
       Title containing UTF-8 characters αβγδ
       --------------------------------------
   @endverbatim
   It may be necessary to set the user locale, first, using the function call
   @code
       std::setlocale(LC_ALL, "en_GB.utf8");
   @endcode
 */
inline std::string underline (std::string text)
{
    std::mbstate_t state = std::mbstate_t();
    const char * mbstr = text.data();
    
    std::size_t wlen = std::mbsrtowcs(nullptr, &mbstr, 0, &state);
    
    text.push_back('\n');
    for (std::size_t i = 0; i < wlen; i++)
        text.push_back('-');
    return text;
}

#endif
