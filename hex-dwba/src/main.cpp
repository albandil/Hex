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

#include <cstdio>
#include <cstdlib>

#include <gsl/gsl_sf.h>

#include "angs.h"
#include "arrays.h"
#include "dwba1.h"
#include "dwba2.h"
#include "hydrogen.h"
#include "potential.h"

int main(int argc, char *argv[])
{
// 	gsl_set_error_handler_off();
	
	if (argc != 8)
	{
		printf("\nusage:\n");
		printf("\tdwba <ni> <li> <nf> <lf> <Ei> <sigmaeps> <order>\n\n");
		return 0;
	}
	
	// extract parameters
	int Ni = strtol(argv[1], 0, 10);
	int Li = strtol(argv[2], 0, 10);
	int Nf = strtol(argv[3], 0, 10);
	int Lf = strtol(argv[4], 0, 10);
	double Ei = strtod(argv[5], 0);
	double ki = sqrt(Ei);
	double kf = sqrt(Ei - 1./(Ni*Ni) + 1./(Nf*Nf));
	double sigmaeps = strtod(argv[6], 0);
	int order = strtol(argv[7], 0, 10);
	
	int MM = (2*Li+1)*(2*Lf+1);		// m⟶m"transition count
	
	cArrays Tdir;
	cArrays Texc;
	cArrays DD;
	cArrays DE;
	cArrays ED;
	cArrays EE;
		
	// accelerators: if a per-lf contribution decays under "accelerator_eps",
	// such term will not be computed anymore for higher lf-contributions
	double accelerator_eps = 1e-8;
	bool compute_DD = true;
	bool compute_DE = true;
	bool compute_ED = true;
	bool compute_EE = true;
	
	// initial and final atomic state
	HydrogenFunction psii(Ni, Li);
	HydrogenFunction psif(Nf, Lf);
	
	// distorting potentials
	DistortingPotential Uf(Nf);
	DistortingPotential Ui(Ni);
	DistortingPotential Ug(1);
	
	for (int lf = 0; ; lf++)
	{
		cArray Tdir_lf(MM), Texc_lf(MM), DD_lf(MM), DE_lf(MM), ED_lf(MM), EE_lf(MM);
		
		printf("lf = %d\n", lf);
		
		DistortedWave chif = Ui.getDistortedWave(kf,lf);
		
		// direct 1e
		if (Ni == Nf and Li == Lf)
		{
			Complex tmat = DWBA1::computeDirect1e(Uf,lf,ki);
			for (int Mi = -Li; Mi <= Li; Mi++)
				for (int Mf = -Lf; Mf <= Lf; Mf++)
					Tdir_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += tmat;
		}
		
		for (int li = 0; ; li++)
		{
			printf("\tli = %d\n", li);
			
			// conserve angular momentum
			if (li < lf - Li - Lf)
				continue;
			
			// conserve angular momentum
			if (li > lf + Li + Lf)
				break;
			
			// conserve parity
			if ((li + Li) % 2 != (lf + Lf) % 2)
				continue;
			
			cArray DD_lf_li(MM), DE_lf_li(MM), ED_lf_li(MM), EE_lf_li(MM);
			DistortedWave chii = Ui.getDistortedWave(ki,li);
			
			// exchange 1e
			if (Li == lf and Lf == li)
			{
				Complex tmat = DWBA1::computeExchange1e(Uf, Ni, Li, ki, Nf, Lf, kf);
				for (int Mi = -Li; Mi <= Li; Mi++)
					for (int Mf = -Lf; Mf <= Lf; Mf++)
						Texc_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += (Mf == 0) ? tmat : 0.;
			}
			
			// direct 2e
			for (int lambda = std::max(abs(Li-Lf),abs(li-lf)); lambda <= std::min(Li+Lf,li+lf); lambda++)
			{
				Complex tmat = DWBA1::computeDirect2e(Uf, lambda, Nf, Lf, kf, lf, Ni, Li, ki, li);
				
				for (int Mi = -Li; Mi <= Li; Mi++)
				{
					for (int Mf = -Lf; Mf <= Lf; Mf++)
					{
						double Gaunts_dir = Gaunt(lambda, Mi - Mf, li, 0, lf, Mi - Mf) * Gaunt(lambda, Mi - Mf, Lf, Mf, Li, Mi);
						Tdir_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += tmat * Gaunts_dir;
					}
				}
			}
			
			// exchange 2e
			for (int lambda = std::max(abs(Li-lf),abs(li-Lf)); lambda <= std::min(Li+lf,li+Lf); lambda++)
			{
				Complex tmat = DWBA1::computeExchange2e(Uf, lambda, Nf, Lf, kf, lf, Ni, Li, ki, li);
				for (int Mi = -Li; Mi <= Li; Mi++)
				{
					for (int Mf = -Lf; Mf <= Lf; Mf++)
					{
						double Gaunts_exc = Gaunt(lambda, -Mf, Li, Mi, lf, Mi - Mf) * Gaunt(lambda, -Mf, Lf, Mf, li, 0);
						Texc_lf[(Mi+Li)*(2*Lf+1)+Mf+Lf] += tmat * Gaunts_exc;
					}
				}
			}
			
			if (order == 1)
				continue;
			
			// DWBA-2 part
			for (int Ln = 0; ; Ln++)
			{
				cArray DD_lf_li_Ln(MM), DE_lf_li_Ln(MM), ED_lf_li_Ln(MM), EE_lf_li_Ln(MM);
				printf("\t\tLn = %d\n", Ln);
				
				DWBA2::DWBA2_Ln (
					Ei, li, lf, ki, kf, Ni, Nf, Li, Lf,
					Ln,
					Ui, Uf, psii, psif, chii, chif,
					DD_lf_li_Ln, DE_lf_li_Ln, ED_lf_li_Ln, EE_lf_li_Ln,
					compute_DD, compute_DE, compute_ED, compute_EE
				);
				
				// update T-matrices
				DD_lf_li = DD_lf_li + DD_lf_li_Ln;
				DE_lf_li = DE_lf_li + DE_lf_li_Ln;
				ED_lf_li = ED_lf_li + ED_lf_li_Ln;
				EE_lf_li = EE_lf_li + EE_lf_li_Ln;
				
				// relative changes of second-order amplitudes
				cArray relchng_singlet = (DD_lf_li_Ln + DE_lf_li_Ln + ED_lf_li_Ln + EE_lf_li_Ln) / (DD_lf_li + DE_lf_li + ED_lf_li + EE_lf_li);
				cArray relchng_triplet = (DD_lf_li_Ln - DE_lf_li_Ln - ED_lf_li_Ln + EE_lf_li_Ln) / (DD_lf_li - DE_lf_li - ED_lf_li + EE_lf_li);
					
				printf("\t\t\tδ %g %g\n", abs(relchng_singlet[0]), abs(relchng_triplet[0]));
					
				if (std::max(max(abs(relchng_singlet)), max(abs(relchng_triplet))) <= sigmaeps)
					break; // OK, converged
				
			} // end for Ln
		} // end for li
		
		// update T-matrices
		Tdir.push_back(Tdir_lf);
		Texc.push_back(Texc_lf);
		DD.push_back(DD_lf);
		DE.push_back(DE_lf);
		ED.push_back(ED_lf);
		EE.push_back(EE_lf);
		
		// convergence check
		rArray sigma_singlet_lf = pow(abs(Tdir_lf + DD_lf + DE_lf + ED_lf + EE_lf + Texc_lf), 2);
		rArray sigma_triplet_lf = pow(abs(Tdir_lf + DD_lf - DE_lf - ED_lf + EE_lf - Texc_lf), 2);
		rArray sigma_singlet = sums(pow(abs(Tdir + DD + DE + ED + EE + Texc), 2));
		rArray sigma_triplet = sums(pow(abs(Tdir + DD - DE - ED + EE - Texc), 2));
		rArray relcng_singlet = sigma_singlet_lf/sigma_singlet;
		rArray relcng_triplet = sigma_triplet_lf/sigma_triplet;
		if (std::max(max(relcng_singlet), max(relcng_triplet)) < sigmaeps)
			break;
		
		// update regulators ("do not unnecessarily refine epsilons")
		double DD_contrib = std::max(max(abs(DD_lf)/sigma_singlet), max(abs(DD_lf)/sigma_triplet));
		double DE_contrib = std::max(max(abs(DE_lf)/sigma_singlet), max(abs(DE_lf)/sigma_triplet));
		double ED_contrib = std::max(max(abs(ED_lf)/sigma_singlet), max(abs(ED_lf)/sigma_triplet));
		double EE_contrib = std::max(max(abs(EE_lf)/sigma_singlet), max(abs(EE_lf)/sigma_triplet));
		if (compute_DD and DD_contrib < accelerator_eps)
		{
			compute_DD = false;
			printf("\tabandoning DD part of DWBA-2\n");
		}
		if (compute_DE and DE_contrib < accelerator_eps)
		{
			compute_DE = false;
			printf("\tabandoning DE part of DWBA-2\n");
		}
		if (compute_ED and ED_contrib < accelerator_eps)
		{
			compute_ED = false;
			printf("\tabandoning ED part of DWBA-2\n");
		}
		if (compute_EE and EE_contrib < accelerator_eps)
		{
			compute_EE = false;
			printf("\tabandoning EE part of DWBA-2\n");
		}
	}
	
	// extract integral cross section for all transitions
	double sumsumsigma = 0;
	printf("\n%g\t", Ei);
	for (int Mi = -Li; Mi <= Li; Mi++)
	{
		double sumsigma = 0;
		for (int Mf = -Lf; Mf <= Lf; Mf++)
		{
			double sigma = 0;
			
			for (int ell = 0; ell < (int)Tdir.size(); ell++)
			{
				Complex tdir = Tdir[ell][(Mi+Li)*(2*Lf+1)+Mf+Lf];
				Complex texc = Texc[ell][(Mi+Li)*(2*Lf+1)+Mf+Lf];
				double tsinglet = abs(tdir + texc);
				double ttriplet = abs(tdir - texc);
				sigma += kf * (0.25 * tsinglet * tsinglet + 0.75 * ttriplet * ttriplet) / ki;
			}
			
			sigma /= 4 * M_PI * M_PI;
			
			printf("%g ", sigma);
			sumsigma += sigma;
		}
		printf("%g ", sumsigma);
		sumsumsigma += sumsigma;
	}
	printf("%g %ld\n", sumsumsigma, Tdir.size());
	
	return EXIT_SUCCESS;
}
