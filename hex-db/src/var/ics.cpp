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

#include <map>
#include <string>
#include <vector>

#include <fftw3.h>
#include <sqlite3.h>

#include "../interpolate.h"
#include "../variables.h"

// -------------------------------------------------------------------------- //

//
// custom function for evaluation of square root within a SQL statement
//

void db_sqrt(sqlite3_context* pdb, int n, sqlite3_value** val)
{
	sqlite3_result_double(pdb, sqrt(sqlite3_value_double(*val)));
}

//
// custom function for evaluation of Gauss-Chebyshev integration of a
// squared Chebyshev expansion
//

void db_gausschebsqr(sqlite3_context* pdb, int n, sqlite3_value** val)
{
	// get blob data as text; reinterpret_cast is save as we are using
	// the low ASCII only
	std::string blob = reinterpret_cast<const char*>(sqlite3_value_text(*val));
	
	// convert text data to binary array
	cArray coeffs;
	coeffs.fromBlob(blob);
	int N = coeffs.size();
	
	// setup FFTW
	cArray mirror_coeffs(4*N+1), evalf(4*N);
	fftw_plan plan = fftw_plan_dft_1d (
		4*N,
		reinterpret_cast<fftw_complex*>(&mirror_coeffs[0]),
		reinterpret_cast<fftw_complex*>(&evalf[0]),
		FFTW_FORWARD,
		0
	);
	
	// mirror oddly around N (3N), evenly around 2N (0,4N)
	for (int k = 0; k < N; k++)
	{
		mirror_coeffs[4*N-k] = mirror_coeffs[k] = coeffs[k];
		mirror_coeffs[2*N+k] = mirror_coeffs[2*N-k] = -coeffs[k];
	}
	
	// integrate
	//    ₁                       n/2-1                
	//   ⌠              dx     2π ===  |       2j+1     |²
	// 2 ⎮ |f(|x|)|² ——————— = —— >    | f(cos(———— π)) |  
	//   ⌡           √(1-x²)   n  ===  |        2n      |
	//  ⁰                         j=0
	// where
	//                        N-1
	//        2j+1       c₀   ===         j+½
	//  f(cos(———— π)) = —— + >   ck cos( ——— kπ )
	//         2n        2    ===          n
	//                        k=1
	// can be evaluated by DCT-III (inverse DCT-II) if full precision is used,
	// i.e. n = N; the result will be stored in odd elements (multiplied by 4)
	
	// evaluate the function using the FFT
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	
	// sum contributions
	double result = 0;
	for (int j = 0; j < N/2; j++)
		result += sqrabs(evalf[2*j+1]);      // (FFTW magic) odd elements only
	result *= 0.0625 * M_PI / coeffs.size(); // (FFTW magic) 1/4²
	
	// use result of the integration
	sqlite3_result_double(pdb, result);
}

// -------------------------------------------------------------------------- //

const std::string IntegralCrossSection::Id = "ics";
const std::string IntegralCrossSection::Description = "Integral cross section.";
const std::vector<std::string> IntegralCrossSection::Dependencies = {
	"ni", "li", "mi",
	"nf", "lf", "mf",
	"L", "S",
	"Ei"
};
const std::vector<std::string> IntegralCrossSection::VecDependencies = { "Ei" };

bool IntegralCrossSection::initialize(sqlitepp::session & db) const
{
	//
	// define SQRT function
	//
	
	sqlite3_create_function (
		db.impl(),
		"sqrt",
		1,              // pass single argument
		SQLITE_UTF8,
		nullptr,
		&db_sqrt,
		nullptr,
		nullptr
	);
	
	//
	// define Gauss-Chebyshev integration of squared Chebyshev expansion
	//
	
	sqlite3_create_function (
		db.impl(),
		"gausschebsqr",
		1,              // pass single argument
		SQLITE_UTF8,
		nullptr,
		&db_gausschebsqr,
		nullptr,
		nullptr
	);
	
	return true;
}

std::vector<std::string> const & IntegralCrossSection::SQL_CreateTable() const
{
	static const std::vector<std::string> cmd = {
		"CREATE TABLE '" + IntegralCrossSection::Id + "' ("
			"ni INTEGER, "
			"li INTEGER, "
			"mi INTEGER, "
			"nf INTEGER, "
			"lf INTEGER, "
			"mf INTEGER, "
			"L  INTEGER, "
			"S  INTEGER, "
			"Ei DOUBLE PRECISION, "
			"sigma DOUBLE PRECISION, "
			"PRIMARY KEY (ni,li,mi,nf,lf,mf,L,S,Ei)"
		")"
	};
	
	return cmd;
}

std::vector<std::string> const & IntegralCrossSection::SQL_Update() const
{
	static const std::vector<std::string> cmd = {
		
		// insert discrete transitions
		
		"INSERT OR REPLACE INTO " + IntegralCrossSection::Id + " "
			"SELECT ni, li, mi, nf, lf, mf, L, S, Ei, "
				"sqrt(Ei-1./(ni*ni)+1./(nf*nf))/sqrt(Ei)*(2*S+1)*SUM(Re_T_ell*Re_T_ell+Im_T_ell*Im_T_ell)/157.91367 " // 16π²
			"FROM " + TMatrix::Id + " "
			"WHERE Ei > 0 AND Ei - 1./(ni*ni) + 1./(nf*nf) > 0 "
			"GROUP BY ni, li, mi, nf, lf, mf, L, S, Ei",
		
		// insert ionization
		
		"INSERT OR REPLACE INTO " + IntegralCrossSection::Id + " "
			"SELECT ni, li, mi, 0, 0, 0, L, S, Ei, "
				"SUM(0.25*(2*S+1)*gausschebsqr(QUOTE(cheb))) "
			"FROM " + IonizationF::Id + " "
			"GROUP BY ni, li, mi, L, S, Ei"
	};
	return cmd;
}

bool IntegralCrossSection::run (
	eUnit Eunits, lUnit Lunits,
	sqlitepp::session & db,
	std::map<std::string,std::string> const & sdata
) const {
	
	// manage units
	double efactor = change_units(Eunits, eUnit_Ry);
	double lfactor = change_units(lUnit_au, Lunits);
	
	// scattering event parameters
	int ni = As<int>(sdata, "ni", Id);
	int li = As<int>(sdata, "li", Id);
	int mi = As<int>(sdata, "mi", Id);
	int nf = As<int>(sdata, "nf", Id);
	int lf = As<int>(sdata, "lf", Id);
	int mf = As<int>(sdata, "mf", Id);
	int  L = As<int>(sdata, "L", Id);
	int  S = As<int>(sdata, "S", Id);
	
	// energies and cross sections
	double E, sigma;
	rArray energies, E_arr, sigma_arr;
	
	// get energy / energies
	try {
		
		// is there a single energy specified using command line ?
		energies.push_back(As<double>(sdata, "Ei", Id));
		
	} catch (std::exception e) {
		
		// are there more energies specified using the STDIN ?
		energies = readStandardInput<double>();
	}
	
	// compose query
	sqlitepp::statement st(db);
	st << "SELECT Ei, sigma FROM " + IntegralCrossSection::Id + " "
			"WHERE ni = :ni "
			"  AND li = :li "
			"  AND mi = :mi "
			"  AND nf = :nf "
			"  AND lf = :lf "
			"  AND mf = :mf "
			"  AND  L = :L  "
			"  AND  S = :S  "
			"ORDER BY Ei ASC",
		sqlitepp::into(E), sqlitepp::into(sigma),
		sqlitepp::use(ni), sqlitepp::use(li), sqlitepp::use(mi),
		sqlitepp::use(nf), sqlitepp::use(lf), sqlitepp::use(mf),
		sqlitepp::use(L), sqlitepp::use(S);
	
	// retrieve data
	while (st.exec())
	{
		E_arr.push_back(E);
		sigma_arr.push_back(sigma);
	}
	
	// write header
	std::cout << this->logo() <<
		"# Integral cross section in " << unit_name(Lunits) << " for\n" <<
		"#     ni = " << ni << ", li = " << li << ", mi = " << mi << ",\n" <<
	    "#     nf = " << nf << ", lf = " << lf << ", mf = " << mf << ",\n" <<
	    "#     L = " << L << ", S = " << S << "\n" <<
	    "# ordered by energy in " << unit_name(Eunits) << "\n" <<
		"#\n" <<
	    "# E\t σ\n";
	
	if (energies[0] < 0.)
	{
		// negative energy indicates full output
		for (size_t i = 0; i < E_arr.size(); i++)
			std::cout << E_arr[i] / efactor << "\t" << sigma_arr[i] * lfactor * lfactor << "\n";
	}
	else
	{
		// threshold for ionization
		double Eion = 1./(ni*ni);
		
		// interpolate
		rArray ics = (efactor * energies.front() < Eion) ? 
			interpolate_real(E_arr, sigma_arr, energies * efactor, o2scl::itp_linear) :
			interpolate_real(E_arr, sigma_arr, energies * efactor, o2scl::itp_cspline);
		
		// output
		for (size_t i = 0; i < energies.size(); i++)
			std::cout << energies[i] << "\t" << ics[i] * lfactor * lfactor << "\n";
	}
	
	return true;
}
