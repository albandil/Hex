//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_RTS_H

// --------------------------------------------------------------------------------- //

#include <memory>
#include <string>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-misc.h"

// --------------------------------------------------------------------------------- //

#define baseClassRunTimeSelectionDefinitions(BASE,CTORARGS) \
    virtual BASE * New CTORARGS const = 0; \
    virtual std::string name () const = 0; \
    static std::unique_ptr<std::vector<BASE*>> RTS_Table; \
    static int Add (BASE * ptr) \
    { \
        if (RTS_Table.get() == nullptr) \
            RTS_Table.reset(new std::vector<BASE*>()); \
        RTS_Table->push_back(ptr); \
        return RTS_Table->size(); \
    } \
    template <class ...Params> static BASE * Choose (std::string str, Params & ...p) \
    { \
        if (RTS_Table.get() == nullptr) \
            HexException("No run-time selectables in " #BASE "."); \
        for (BASE *obj : *RTS_Table) \
        { \
            if (obj->name() == str or str == "any") \
                return obj->New(p...); \
        } \
        std::cout << "Error!" << std::endl; \
        std::cout << "  The selectable \"" << str << "\" is not available in " #BASE "." << std::endl; \
        std::cout << "  The program may not be compiled with support for that object." << std::endl; \
        std::cout << "  The available selectables are: "; \
        for (BASE *obj : *RTS_Table) \
            std::cout << "\"" << obj->name() << "\" "; \
        std::cout << std::endl << std::endl; \
        std::exit(1); \
    }

#define defineBaseClassRunTimeSelectionTable(BASE) \
    std::unique_ptr<std::vector<BASE*>> BASE::RTS_Table; \

#define derivedClassRunTimeSelectionDefinitions(BASE,CTORARGS,TYPE,ARGS,NAME) \
    virtual BASE * New CTORARGS const { return new TYPE ARGS ; } \
    virtual std::string name () const { return NAME ; };

#define addClassToParentRunTimeSelectionTable(BASE,TYPE) \
    int add_##TYPE##_to_##BASE = BASE::Add(new TYPE());

// --------------------------------------------------------------------------------- //

#endif // HEX_RTS_H
