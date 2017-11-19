//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2017, Jakub Benda, Charles University in Prague                    //
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

#include <iostream>
#include <x86intrin.h>

// --------------------------------------------------------------------------------- //

#ifdef __linux__
    #include <sys/stat.h>
    #include <unistd.h>
    
    #ifdef __GNUC__
        #include <execinfo.h>
        #include <cxxabi.h>
    #endif
#endif

// --------------------------------------------------------------------------------- //

#include "hex-memory.h"
#include "hex-misc.h"

// --------------------------------------------------------------------------------- //

#ifdef _WIN32
    #include <windows.h>
#endif

// --------------------------------------------------------------------------------- //

std::string get_executable_path ()
{
    // allocate output string
    std::string path;
    path.resize(256);
    
#if defined(__linux__)
    std::size_t n = path.size();
    while (n == path.size())
    {
        path.resize(2 * n);
        n = readlink("/proc/self/exe", &path[0], path.size());
    }
    path.resize(n);
#elif defined(_WIN32)
    GetModuleFileNameA(GetModuleHandle(nullptr), &path[0], sizeof(path));
#endif
    
    return path;
}

void print_stack_trace ()
{
#if (defined(__linux__) && defined(__GNUC__))
    
    std::cerr << "Stack trace:" << std::endl;
    
    // retrieve the backtrace
    void* addrlist[65];
    int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));
    if (addrlen == 0)
    {
        std::cerr << "  <empty, possibly corrupt>" << std::endl;
#ifdef NDEBUG
        std::exit(EXIT_FAILURE);
#else
        std::abort();
#endif
    }
    
    // translate the backtrace
    char** symbollist = backtrace_symbols(addrlist, addrlen);
    
    // for all lines of the list
    for (int i = 1; i < addrlen; i++)
    {
        // convert line to string
        std::string line (symbollist[i]);
        
        // get module name
        std::string module;
        std::size_t module_begin = 0, module_end = line.find('(');
        if (module_end != std::string::npos)
            module = line.substr(module_begin, module_end - module_begin);
        
        std::cerr << " #" << i << " " << module << " : ";
        
        // get function name
        std::string function;
        std::size_t function_begin = line.find('(') + 1, function_end = line.find('+');
        if (function_begin != std::string::npos and function_end != std::string::npos and function_end > function_begin)
            function = line.substr(function_begin, function_end - function_begin);
        
        // demangle function name
        int status;
        char* funcname = nullptr;
        std::size_t funcnamesize = 0;
        char* ret = abi::__cxa_demangle(function.c_str(), funcname, &funcnamesize, &status);
        funcname = (status == 0 ? ret : nullptr);
        
        if (not function.empty() and funcname != nullptr)
            std::cerr << "[" << funcname << "] ";
        else if (not function.empty())
            std::cerr << "[" << function << "] ";
        else
            std::cerr << "[unknown] ";
        
        // get address
        std::string address;
        std::size_t address_begin = line.find('[') + 1, address_end = line.find(']');
        if (address_begin != std::string::npos and address_end != std::string::npos and address_end > address_begin)
            address = line.substr(address_begin, address_end - address_begin);
        
        // translated address into line in source file
        std::string srcfile;
        
#if _POSIX_C_SOURCE >= 2
        // compose the command line
        std::string cmd = format
        (
            "addr2line --exe=%s %s 2>&1",
            get_executable_path().c_str(),
            address.c_str()
        );
        
        // launch "addr2line" (if available) and connect its output as a pipe
        FILE * pipe = popen(cmd.c_str(), "r");
        if (pipe != nullptr)
        {
            // read all characters from the pipe
            while (not std::feof(pipe))
                srcfile.push_back(std::getc(pipe));
            
            // check whether the executable has been found
            if (srcfile.find("No such file") != std::string::npos)
            {
                // NO : Erase the error message.
                srcfile.clear();
            }
            else
            {
                // YES : keep only the file name ...
                std::size_t pos = srcfile.rfind('/');
                if (pos != std::string::npos)
                    srcfile = srcfile.substr(pos + 1, std::string::npos);
                
                // ... and erase trailing white/invalid characters
                while (not srcfile.empty() and (std::isspace(srcfile.back()) or srcfile.back() < 0))
                        srcfile.pop_back();
            }
        }
        pclose(pipe);
#endif
        
        // print info
        if (srcfile.empty())
            std::cerr << "@ [" << address << "]" << std::endl;
        else
            std::cerr << "@ " << srcfile << std::endl;
        
        // cleanup
        std::free(funcname);
    }
    std::cerr << std::endl;
    std::free(symbollist);
#endif
}

std::string nice_size (std::size_t bytes)
{
    std::size_t kiB = 1024;
    std::size_t MiB = kiB * kiB;
    std::size_t GiB = kiB * kiB * kiB;
    std::size_t TiB = kiB * kiB * kiB * kiB;
    
    if (bytes < kiB) return format("%d B", bytes);
    if (bytes < MiB) return format("%d kiB", bytes / kiB);
    if (bytes < GiB) return format("%d MiB", bytes / MiB);
    if (bytes < TiB) return format("%d GiB", bytes / GiB);
    
    return format("%d TiB", bytes / TiB);
}

std::string current_time ()
{
    std::time_t result = std::time(nullptr);
    return std::asctime(std::localtime(&result));
}

void create_directory (std::string dir)
{
#ifdef __linux__
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#elif defined(_WIN32)
    CreateDirectoryA(dir.c_str(), NULL);
#endif
}

#if defined(__AVX2__) && defined(__FMA__)
void cmul2xd
(
    std::complex<double>       * const restrict C,
    std::complex<double> const * const restrict A,
    std::complex<double> const * const restrict B
)
{
    // Multiplies two 2-component complex vectors (i.e. each containing 4 double-precision real values).
    //
    // [ a1 ]   [ b1 ]   [ a1 ] * [ b1 ] - [ a2 ] * [ b2 ]
    // [ a2 ]   [ b2 ]   [ a2 ] * [ b1 ] + [ a1 ] * [ b2 ]
    // [ -- ] * [ -- ] = [ -- ]   [ -- ]   [ -- ]   [ -- ]
    // [ a3 ]   [ b3 ]   [ a3 ] * [ b3 ] - [ a4 ] * [ b4 ]
    // [ a4 ]   [ b4 ]   [ a4 ] * [ b3 ] + [ a3 ] * [ b4 ]
    
    __m256d A1234 = _mm256_loadu_pd(reinterpret_cast<const double*>(A));
    __m256d A2143 = _mm256_permute_pd(A1234, 0b01010101);
    
    __m256d B1234 = _mm256_loadu_pd(reinterpret_cast<const double*>(B));
    __m256d B1133 = _mm256_unpacklo_pd(B1234, B1234);
    __m256d B2244 = _mm256_unpackhi_pd(B1234, B1234);
    
    __m256d C1234 = _mm256_loadu_pd(reinterpret_cast<const double*>(C));
    
    C1234 = _mm256_fmaddsub_pd(A2143, B2244, C1234);
    C1234 = _mm256_fmaddsub_pd(A1234, B1133, C1234);
    
    _mm256_storeu_pd(reinterpret_cast<double*>(C), C1234);
}
#endif

#ifdef __AVX512F__
void cmul4xd
(
    std::complex<double>       * const restrict C,
    std::complex<double> const * const restrict A,
    std::complex<double> const * const restrict B
)
{
    // Multiplies two 4-component complex vectors (i.e. each containing 8 double-precision real values).
    //
    // [ a1 ]   [ b1 ]   [ a1 ] * [ b1 ] - [ a2 ] * [ b2 ]
    // [ a2 ]   [ b2 ]   [ a2 ] * [ b1 ] + [ a1 ] * [ b2 ]
    // [ -- ]   [ -- ]   [ -- ]   [ -- ]   [ -- ]   [ -- ]
    // [ a3 ]   [ b3 ]   [ a3 ] * [ b3 ] - [ a4 ] * [ b4 ]
    // [ a4 ]   [ b4 ]   [ a4 ] * [ b3 ] + [ a3 ] * [ b4 ]
    // [ -- ] * [ -- ] = [ -- ]   [ -- ]   [ -- ]   [ -- ]
    // [ a5 ]   [ b5 ]   [ a5 ] * [ b5 ] - [ a6 ] * [ b6 ]
    // [ a6 ]   [ b6 ]   [ a6 ] * [ b5 ] + [ a5 ] * [ b6 ]
    // [ -- ]   [ -- ]   [ -- ]   [ -- ]   [ -- ]   [ -- ]
    // [ a7 ]   [ b7 ]   [ a7 ] * [ b7 ] - [ a8 ] * [ b8 ]
    // [ a8 ]   [ b8 ]   [ a8 ] * [ b7 ] + [ a7 ] * [ b8 ]
    
    __m512d A12345678 = _mm512_loadu_pd(reinterpret_cast<const double*>(A));
    __m512d A21436587 = _mm512_permute_pd(A12345678, 0b01010101);
    
    __m512d B12345678 = _mm512_loadu_pd(reinterpret_cast<const double*>(B));
    __m512d B11335577 = _mm512_unpacklo_pd(B12345678, B12345678);
    __m512d B22446688 = _mm512_unpackhi_pd(B12345678, B12345678);
    
    __m512d C12345678 = _mm512_loadu_pd(reinterpret_cast<const double*>(C));
    
    C12345678 = _mm512_fmaddsub_pd(A21436587, B22446688, C12345678);
    C12345678 = _mm512_fmaddsub_pd(A12345678, B11335577, C12345678);
    
    _mm512_storeu_pd(reinterpret_cast<double*>(C), C12345678);
}
#endif
