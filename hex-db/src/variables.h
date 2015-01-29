//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_VARIABLES
#define HEX_VARIABLES

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include <sqlitepp/sqlitepp.hpp>

#include "arrays.h"
#include "hex-db.h"
#include "vec3d.h"

/**
 * @brief Energy units change.
 * 
 * Returns factor that can be used to transform from the unit system A
 * to the unit system B.
 */
double change_units (eUnit A, eUnit B);

/**
 * @brief Lengths units change.
 * 
 * Returns factor that can be used to transform from the unit system A
 * to the unit system B.
 */
double change_units (lUnit A, lUnit B);

/**
 * @brief Angular units change.
 * 
 * Returns factor that can be used to transform from the unit system A
 * to the unit system B.
 */
double change_units (aUnit A, aUnit B);

/**
 * @brief Energy unit name.
 * 
 * Return energy unit name as string.
 */
std::string unit_name (eUnit u);

/**
 * @brief Length unit name.
 * 
 * Return length unit name as string.
 */
std::string unit_name (lUnit u);

/**
 * @brief Length unit name.
 * 
 * Return length unit name as string.
 */
std::string unit_name (aUnit u);

/**
 * Write out std::pair.
 */
inline std::ostream & operator << (std::ostream & os, std::pair<geom::vec3d,geom::vec3d> const & p)
{
    os << p.first << " " << p.second;
    return os;
}

/**
 * Read in std::pair.
 */
inline std::istream & operator >> (std::istream & is, std::pair<geom::vec3d,geom::vec3d> & p)
{
    is >> p.first;
    is >> p.second;
    return is;
}

/**
 * Read data from standard input.
 */
template<typename T> std::vector<T> readStandardInput ()
{
    std::vector<T> data;
    
    T x;
    while (not std::cin.eof())
    {
        std::cin >> std::ws;
        std::cin >> x;
        std::cin >> std::ws;
        data.push_back(x);
    }
    
    return data;
}

/**
 * @brief Convert dictionary entry to a numeric type.
 * 
 * Being given a dictionary (= string-string map) and a keyword,
 * the function finds a correct entry and returns its value converted
 * to the template datatype.
 * @param dict Dictionary to search in.
 * @param keyword Entry to look for.
 * @param name Identification of the calling authority for use in error
 *             message if the entry is not find.
 */
template <typename T> T Conv
(
    std::map<std::string,std::string> const & dict,
    std::string const & keyword,
    std::string const & name
)
{
    // check existence of the keyword
    std::map<std::string,std::string>::const_iterator it = dict.find(keyword);
    if (it == dict.end())
        throw exception ("ERROR: \"%s\" requires specifying the parameter \"--%s\"!\n", name.c_str(), keyword.c_str());
    
    // convert to int
    T x;
    std::istringstream ss(it->second);
    ss >> x;
    return x;
}

/**
 * @brief Base class for scatering variables.
 * 
 * This is the heritage base for all scattering quantities that can be
 * computed by hex-db.
 */
class Variable
{
public:
    
    /// destructor
    virtual ~Variable() {}
    
    // getters
    
    /// String identification of the variable.
    virtual std::string const & id() const = 0;
    
    /// Longer description text for use in program help.
    virtual std::string const & description() const = 0;
    
    /// SQL statements that create the required table, or empty vector if not needed.
    virtual std::vector<std::string> const & SQL_CreateTable() const = 0;
    
    /// SQL statements that update the table after insetion of new data.
    virtual std::vector<std::string> const & SQL_Update() const = 0;
    
    /// List of all scattering event parameters (with description) that have to be specified by user.
    virtual std::vector<std::pair<std::string,std::string>> const & deps() const = 0;
    
    /// List of vectorizable scattering event parameters that have to be specified by user.
    virtual std::vector<std::string> const & vdeps() const = 0;

    // others
    
    /// initialize (e.g.) by defining external routines for SQLite
    virtual bool initialize (sqlitepp::session & db) const = 0;
    
    /// write out requested data
    virtual bool run
    (
        std::map<std::string,std::string> const & params
    ) const = 0;
};

/**
 * @brief List of variables.
 * 
 * All available variables are elements of the list. They are being added in
 * the constructor.
 */
class VariableList
{
public:
    
    VariableList();
    ~VariableList();
    
    Variable const * get (std::string const & id) const;
    
    //
    // STL vector interface
    //
    typedef std::vector<Variable*>::iterator iterator;
    typedef std::vector<Variable*>::const_iterator const_iterator;
    inline iterator begin() { return list.begin(); }
    inline iterator end()   { return list.end();   }
    inline const_iterator begin() const { return list.begin(); }
    inline const_iterator end()   const { return list.end();   }
    
private:
    
    std::vector<Variable*> list;
};

/**
 * @brief Create class for a given variable name.
 * 
 * Creates class declaration so that only definition of the unique 
 * members is necessary. This method allows easy addition of new
 * scattering variables; one just needs to add the line
 * @code
 *     AddNewVariableClass(NewVariableClassName);
 * @endcode
 * at the end of "variables.h", add a line
 * @code
 *     new NewVariableClassName,
 * @endcode
 * among others into the constructor of VariableList in "variables.cpp"
 * and implement the methods in a new CPP file, which is then linked
 * to the program. It is not necessary to write a custom header which needs
 * to be included by every other source file.
 */
#define AddNewVariableClass(ClassName) \
    class ClassName : public Variable \
{ \
    public: \
\
        bool initialize (sqlitepp::session & db) const; \
\
        static const std::string Id; \
        std::string const & id () const { return Id; } \
\
        static const std::string Description; \
        std::string const & description () const { return Description; } \
\
        static const std::vector<std::pair<std::string,std::string>> Dependencies; \
        std::vector<std::pair<std::string,std::string>> const & deps () const { return Dependencies; } \
\
        static const std::vector<std::string> VecDependencies; \
        std::vector<std::string> const & vdeps () const { return VecDependencies; } \
\
        std::vector<std::string> const & SQL_CreateTable () const; \
        std::vector<std::string> const & SQL_Update () const; \
\
        bool run \
        ( \
            std::map<std::string,std::string> const & params \
        ) const; \
};

/// Create class for T-matrix (tmat)
AddNewVariableClass(TMatrix);

/// Create class for scattering amplitude (scatamp)
AddNewVariableClass(ScatteringAmplitude);

/// Create class for scattering amplitude with arbitrary impact direction (scatamp-dir)
AddNewVariableClass(ScatteringAmplitudeDir);

/// Create class for differential cross section (dcs)
AddNewVariableClass(DifferentialCrossSection);

/// Create class for integral cross section (ics)
AddNewVariableClass(IntegralCrossSection);

/// Create class for complete cross section (ccs)
AddNewVariableClass(CompleteCrossSection);

/// Create class for extrapolated cross section (xcs)
AddNewVariableClass(ExtrapolatedCrossSection);

/// Create class for collision strength (colls)
AddNewVariableClass(CollisionStrength);

/// Create class for momentum transfer (momtf)
AddNewVariableClass(MomentumTransfer);

/// Create class for total cross section (tcs)
AddNewVariableClass(TotalCrossSection);

/// Create class for ionization amplitude radial part (ionf)
AddNewVariableClass(IonizationF);

/// Create class for the angle dependent full second Born T-matrix (bornf)
AddNewVariableClass(BornFullTMatrix);

/// Create class for ionization amplitude (ionamp)
AddNewVariableClass(IonizationAmplitude);

/// Triple ionization differential cross section (tdcs).
AddNewVariableClass(TripleDifferentialCrossSection);

/// Create class for stokes parameters (stokes)
AddNewVariableClass(StokesParameters);

/// Create class for spin asymmetry (asy)
AddNewVariableClass(SpinAsymmetry);

/// Create class for Born T-matrix (tmatb)
AddNewVariableClass(TMatrixB);

/// forward declaration of the variable list "vlist"
extern VariableList vlist;

/// forward declaration of the global database handle
extern sqlitepp::session db;

#endif
