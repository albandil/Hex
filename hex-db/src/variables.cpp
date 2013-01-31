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
		new TotalCrossSection
	})
{}

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
