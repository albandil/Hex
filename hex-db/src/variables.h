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

#ifndef HEX_VARIABLES
#define HEX_VARIABLES

#include <string>
#include <map>
#include <vector>

#include <sqlitepp/sqlitepp.hpp>

#include "arrays.h"

/**
 * \brief Convert dictionary entry to a numeric type.
 * 
 * Being given a dictionary (= string-string map) and a keyword,
 * the function finds a correct entry and returns its value converted
 * to the template datatype.
 * \param dict Dictionary to search in.
 * \param keyword Entry to look for.
 * \param name Identification of the calling authority for use in error
 *             message if the entry is not find.
 */
template <typename T> T As (
	std::map<std::string,std::string> const & dict,
	std::string const & keyword,
	std::string const & name
)
{
	// check existence of the keyword
	std::map<std::string,std::string>::const_iterator it = dict.find(keyword);
	if (it == dict.end())
	{
		std::cerr << "ERROR: \"" << name << "\" requires specifying the parameter \"--" << keyword << "\"!\n";
		exit(-1);
	}
	
	// convert to int
	T x;
	std::stringstream ss(it->second);
	ss >> x;
	return x;
}

/**
 * \brief Base class for scatering variables.
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
	
	/// SQL statement that creates the required table, or empty string if not needed.
	virtual std::string const & SQL_CreateTable() const = 0;
	
	/// SQL statement that updates the table after insetion of new data.
	virtual std::string const & SQL_Update() const = 0;
	
	/// List of scattering event parameters that have to be specified by user.
	virtual std::vector<std::string> const & dependencies() const = 0;

	// others
	
	/// write out requested data
	virtual bool run (
		sqlitepp::session & db,
		std::map<std::string,std::string> const & data1,
		rArray const & data2
	) const = 0;
	
	/// Returns the program logo for use in output.
	std::string logo() const;
};

/**
 * \brief List of variables.
 * 
 * All available variables are elements of the list. They are being added in
 * the constructor.
 */
class VariableList
{
public:
	
	VariableList();
	~VariableList();
	
	Variable const * const get (std::string const & id) const;
	
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
 * \brief Create class for a given variable name.
 * 
 * Creates class declaration so that only definition of the unique 
 * members is necessary. This method allows easy addition of new
 * scattering variables; one just needs to add the line
 * \code
 * AddNewVariableClass(NewVariableClassName);
 * \endcode
 * at the end of "variables.h", add a line
 * \code
 * new NewVariableClassName,
 * \endcode
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
		static const std::string Id; \
		std::string const & id() const { return Id; } \
\
		static const std::string Description; \
		std::string const & description() const { return Description; } \
\
		static const std::vector<std::string> Dependencies; \
		std::vector<std::string> const & dependencies() const { return Dependencies; } \
\
		std::string const & SQL_CreateTable() const; \
		std::string const & SQL_Update() const; \
\
		bool run ( \
			sqlitepp::session & db, \
			std::map<std::string,std::string> const & data1, \
			rArray const & data2 \
		) const; \
};

/// Create class for T-matrix (tmat)
AddNewVariableClass(TMatrix);

/// Create class for scattering amplitude (scatamp)
AddNewVariableClass(ScatteringAmplitude);

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

/// forward declaration of the variable list "vlist"
extern VariableList vlist;

#endif
