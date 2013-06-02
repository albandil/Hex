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

#include <algorithm>

#include "variables.h"
#include "vec3d.h"

VariableList::VariableList() :
	list ({
		new TMatrix,
		new ScatteringAmplitude,
		new DifferentialCrossSection,
		new IntegralCrossSection,
		new CompleteCrossSection,
		new ExtrapolatedCrossSection,
		new CollisionStrength,
		new MomentumTransfer,
		new TotalCrossSection,
	    new IonizationF,
		new IonizationAmplitude,
	    new TripleDifferentialCrossSection,
		new StokesParameters
	})
{
	
}

VariableList::~VariableList()
{
	// delete list items
	for (Variable* var : list)
		delete var;
}

Variable const * const VariableList::get (std::string const & id) const
{
	// search the list for a (pointer to a) variable with the correct ID
	std::vector<Variable*>::const_iterator it = std::find_if (
		list.begin(),
		list.end(),
		[&](Variable* const & constptr) -> bool {
			return constptr->id() == id;
		}
	);
	
	// return pointer to the variable
	if (it != list.end())
		return *it;
	else
		return nullptr;
}

std::string Variable::logo() const
{
	return std::string (
		
		"#\n"
		"#       / /   / /    __    \\ \\  / /\n"
		"#      / /__ / /   / _ \\    \\ \\/ /\n"
		"#     /  ___  /   | |/_/    / /\\ \\\n"
		"#    / /   / /    \\_\\      / /  \\ \\\n"
		"#\n"
		"#             UK MFF (c) 2013\n"
		"#\n"
		
	);
}

double change_units(eUnit A, eUnit B){
	// no change
	if (A == B)
		return 1.;
	
	double ufactor = 1.;
	
	// transform to Rydbergs
	if (A == eUnit_au)
		ufactor *= 2.;
	if (A == eUnit_eV)
		ufactor *= 1./13.605692;
	
	// tranform from Rydbergs
	if (B == eUnit_au)
		ufactor *= 0.5;
	if (B == eUnit_eV)
		ufactor *= 13.605692;
	
	return ufactor;
}

double change_units(lUnit A, lUnit B)
{
	// no change
	if (A == B)
		return 1.;
	
	double ufactor = 1.;
	
	// transform to a.u.
	if (A == lUnit_cgs)
		ufactor *= 1./5.29177211e-9;
	
	// tranform from a.u.
	if (B == lUnit_cgs)
		ufactor *= 5.29177211e-9;
	
	return ufactor;
}

std::string unit_name(eUnit u)
{
	switch (u)
	{
		case eUnit_au:
			return std::string("a.u.");
		case eUnit_eV:
			return std::string("eV");
		case eUnit_Ry:
			return std::string("Ry");
		default:
			return std::string("");
	}
}

std::string unit_name(lUnit u)
{
	switch (u)
	{
		case lUnit_au:
			return std::string("a.u.");
		case lUnit_cgs:
			return std::string("CGS");
		default:
			return std::string("");
	}
}
