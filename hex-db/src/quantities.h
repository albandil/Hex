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

#ifndef HEX_DB_QUANTITIES_H
#define HEX_DB_QUANTITIES_H

// --------------------------------------------------------------------------------- //

#include <iostream>
#include <string>
#include <map>
#include <vector>

// --------------------------------------------------------------------------------- //

#include <sqlitepp/sqlitepp.hpp>

// --------------------------------------------------------------------------------- //

class ScatteringQuantity
{
    public:
        
        // object management
        
            /// object factory
            virtual ScatteringQuantity* New () = 0;
        
            /// destructor
            virtual ~ScatteringQuantity () {}
        
        // access
            
            /// String identification of the variable.
            virtual std::string id () const = 0;
            
            /// Longer description text for use in program help.
            virtual std::string description () const = 0;
            
            /// List of all scattering event parameters (with description) that have to be specified by user.
            virtual std::vector<std::pair<std::string,std::string>> deps () const = 0;
            
            /// List of vectorizable scattering event parameters that have to be specified by user.
            virtual std::vector<std::string> vdeps () const = 0;

        // database interface
        
            /// initialize (e.g.) by defining external routines for SQLite
            virtual bool initialize (sqlitepp::session & db) const = 0;
            
            /// SQL statements that create the required table, or empty vector if not needed.
            virtual bool createTable (sqlitepp::session & db) const = 0;
            
            /// SQL statements that update the table after insetion of new data.
            virtual bool updateTable (sqlitepp::session & db) const = 0;
            
            /// write out requested data
            virtual bool run
            (
                sqlitepp::session & db,
                std::map<std::string,std::string> const & params
            ) const = 0;
};

// --------------------------------------------------------------------------------- //

extern std::shared_ptr<std::vector<ScatteringQuantity*>> quantities;

bool register_new_quantity (ScatteringQuantity* Q);

// --------------------------------------------------------------------------------- //

// define new scattering quantity
#define createNewScatteringQuantity(Q) \
class Q : public ScatteringQuantity \
{ \
    public: \
        std::string id () const; \
        std::string description () const; \
        std::vector<std::pair<std::string,std::string>> deps () const; \
        std::vector<std::string> vdeps () const; \
        bool initialize (sqlitepp::session & db) const; \
        bool createTable (sqlitepp::session & db) const; \
        bool updateTable (sqlitepp::session & db) const; \
        bool run \
        ( \
            sqlitepp::session & db, \
            std::map<std::string,std::string> const & params \
        ) const; \
}; \
bool Q##_##registered = register_new_quantity(new Q())

// --------------------------------------------------------------------------------- //

#endif
