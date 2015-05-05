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

#include <iostream>

#if (defined(__linux__) && defined(__GNUC__))
    #include <execinfo.h>
    #include <unistd.h>
    #include <cxxabi.h>
#endif

#include "misc.h"

#ifdef _WIN32
    #include <windows.h>
#endif

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
    GetModuleFileNameA(GetModuleHandle(nullptr), path, sizeof(path));
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
        std::terminate();
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
