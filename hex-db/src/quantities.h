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
    protected:
        
        /// Database session.
        sqlitepp::session * db;
    
    public:
        
        // properties
            
            /// String identification of the variable.
            virtual std::string name () = 0;
            
            /// Longer description text for use in program help.
            virtual std::string description () = 0;
            
            /// List of all scattering event parameters (with description) that have to be specified by user.
            virtual std::vector<std::pair<std::string,std::string>> params () = 0;
            
            /// List of vectorizable scattering event parameters that have to be specified by user.
            virtual std::vector<std::string> vparams () = 0;
            
            /// List of other scattering quantities that need to be initialized/updated/etc before this class
            /// because this class makes use of their data.
            virtual std::vector<std::string> dependencies ()
            {
                return std::vector<std::string>();
            }
            
        // database interface
            
            sqlitepp::session & session ()
            {
                return *db;
            }
            
            /// initialize (e.g.) by defining external routines for SQLite
            virtual bool initialize (sqlitepp::session & db)
            {
                this->db = &db;
                return true;
            }
            
            /// SQL statements that create the required table, or empty vector if not needed.
            virtual bool createTable ()
            {
                return true;
            }
            
            /// SQL statements that update the table after insetion of new data.
            virtual bool updateTable ()
            {
                return true;
            }
            
            /// write out requested data
            virtual bool run (std::map<std::string,std::string> const & params) = 0;
};

// --------------------------------------------------------------------------------- //

extern std::shared_ptr<std::vector<ScatteringQuantity*>> quantities;

bool register_new_quantity (ScatteringQuantity* Q);

ScatteringQuantity * get_quantity (std::string name);

// --------------------------------------------------------------------------------- //

// define new scattering quantity
#define createNewScatteringQuantity(Q) \
class Q : public ScatteringQuantity \
{ \
    public: \
        std::string name (); \
        std::string description (); \
        std::vector<std::pair<std::string,std::string>> params (); \
        std::vector<std::string> vparams (); \
        std::vector<std::string> dependencies (); \
        bool initialize (sqlitepp::session & db); \
        bool createTable (); \
        bool updateTable (); \
        bool run (std::map<std::string,std::string> const & params); \
}; \
bool Q##_##registered = register_new_quantity(new Q())

// --------------------------------------------------------------------------------- //

#endif
