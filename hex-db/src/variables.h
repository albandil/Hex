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

class Variable
{
public:
	
	// destructor
	virtual ~Variable() {}
	
	// getters
	virtual std::string const & id() const = 0;
	virtual std::string const & description() const = 0;
	virtual std::string const & SQL_Update() const = 0;
	virtual std::string const & SQL_CreateTable() const = 0;
	virtual std::vector<std::string> const & dependencies() const = 0;

	// write out requested data
	virtual bool run (
		sqlitepp::session & db,
		std::map<std::string,std::string> const & data1,
		rArray const & data2
	) const = 0;
	
	std::string logo() const;
};

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

AddNewVariableClass(TMatrix);
AddNewVariableClass(ScatteringAmplitude);
AddNewVariableClass(DifferentialCrossSection);
AddNewVariableClass(IntegralCrossSection);
AddNewVariableClass(CompleteCrossSection);
AddNewVariableClass(ExtrapolatedCrossSection);
AddNewVariableClass(CollisionStrength);
AddNewVariableClass(MomentumTransfer);
AddNewVariableClass(TotalCrossSection);
	
extern VariableList vlist;

#endif
