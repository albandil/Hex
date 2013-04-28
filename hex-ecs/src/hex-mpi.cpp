/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2012                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <complex>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <omp.h>
#include <mpi.h>
#include <H5Cpp.h>

#include "amplitudes.h"
#include "angs.h"
#include "arrays.h"
#include "bspline.h"
#include "input.h"
#include "moments.h"
#include "slater.h"
#include "spmatrix.h"
#include "special.h"

#ifndef SuiteSparse_time
double SuiteSparse_time()
{
	return 0.;
}
#endif

int main(int argc, char* argv[])
{
	MPI::Init();
	int Nproc = MPI::COMM_WORLD.Get_size();
	int iproc = MPI::COMM_WORLD.Get_rank();
	
	// Preparations ------------------------------------------------------- //
	//
	
	gsl_set_error_handler_off();
	
	// variables that can be set by the user from the command line
	FILE* inputfile = 0;					// input file
	
	// get input from command line
	parse_command_line(argc, argv, inputfile);
	
	if (iproc == 0)
	{
		// check input file
		if (inputfile == 0 and (inputfile = fopen("hex.inp", "r")) == 0)
		{
			fprintf(stderr, "Input error: Cannot open the file \"hex.inp\".\n");
			abort();
		}
	}
	
	// create main variables
	int order;							// B-spline order
	double R0, Rmax, ecstheta;			// B-spline knot sequence
	rArray real_knots;					// real knot sequence
	rArray complex_knots;				// complex knot sequence
	int ni, minli, maxli, maxnf, maxlf;	// atomic state parameters
	int L;								// global conserved variables
	int maxell;							// angular momentum limits
	rArray Ei;							// energies in Rydbergs
	double B = 0;						// magnetic field in a.u.
	int real_knots_size, complex_knots_size, Ei_size;
	
	if (iproc == 0)
	{
		// get input from input file
		parse_input_file(
			inputfile, order, R0, ecstheta, Rmax,
			real_knots, complex_knots, ni, maxnf, 
			minli, maxli, maxlf,
			L, maxell, Ei, B
		);
		
		real_knots_size = real_knots.size();
		complex_knots_size = complex_knots.size();
		Ei_size = Ei.size();
	}
	
	// broadcast input values
	MPI::COMM_WORLD.Bcast(&order,    1, MPI::INTEGER, 0);
	MPI::COMM_WORLD.Bcast(&ni,       1, MPI::INTEGER, 0);
	MPI::COMM_WORLD.Bcast(&minli,    1, MPI::INTEGER, 0);
	MPI::COMM_WORLD.Bcast(&maxli,    1, MPI::INTEGER, 0);
	MPI::COMM_WORLD.Bcast(&maxnf,    1, MPI::INTEGER, 0);
	MPI::COMM_WORLD.Bcast(&maxlf,    1, MPI::INTEGER, 0);
	MPI::COMM_WORLD.Bcast(&L,        1, MPI::INTEGER, 0);
	MPI::COMM_WORLD.Bcast(&maxell,   1, MPI::INTEGER, 0);
	MPI::COMM_WORLD.Bcast(&Ei_size,  1, MPI::INTEGER, 0);
	MPI::COMM_WORLD.Bcast(&B,        1, MPI::DOUBLE,  0);
	MPI::COMM_WORLD.Bcast(&R0,       1, MPI::DOUBLE,  0);
	MPI::COMM_WORLD.Bcast(&Rmax,     1, MPI::DOUBLE,  0);
	MPI::COMM_WORLD.Bcast(&ecstheta, 1, MPI::DOUBLE,  0);
	MPI::COMM_WORLD.Bcast(&real_knots_size,    1, MPI::INTEGER, 0);
	MPI::COMM_WORLD.Bcast(&complex_knots_size, 1, MPI::INTEGER, 0);
	if (iproc != 0)
	{
		real_knots.resize(real_knots_size);
		complex_knots.resize(complex_knots_size);
		Ei.resize(Ei_size);
	}
	MPI::COMM_WORLD.Bcast(&Ei[0],            Ei_size,            MPI::DOUBLE,  0);
	MPI::COMM_WORLD.Bcast(&real_knots[0],    real_knots_size,    MPI::DOUBLE,  0);
	MPI::COMM_WORLD.Bcast(&complex_knots[0], complex_knots_size, MPI::DOUBLE,  0);
	
	// set MPI environment
	if (Nproc != maxell + 1)
	{
		printf("! Wrong number of procs.\n");
		printf("!    maxell + 1 = %d\n", maxell + 1);
		printf("!    no. of procs = %d\n", Nproc);
		
		abort();
	}
	
	
	// shorthand for total number of energies
	unsigned Nenergy = Ei.size();
	
	// adjust angular momentum limits
	if (minli < 0) minli = 0;
	if (maxli < 0) maxli = maxell;
	if (maxlf < 0) maxlf = maxell;
	if (maxli >= ni) maxli = ni - 1;
	
	// projectile momenta, initial and final
	std::vector<double> ki(Nenergy), kf(Nenergy);		// in a.u.
	std::transform(
		Ei.begin(), Ei.end(),
		ki.begin(),
		[ = ](double E) -> double { return sqrt(E); }
	);
	// --------------------------------------------------------------------- //
	
	
	// Setup B-spline environment ------------------------------------------ //
	//
	setup_knot_sequence(order, real_knots, R0, ecstheta, complex_knots, Rmax);
	
	//
	// --------------------------------------------------------------------- //
	
	fprintf(stdout, "B-spline total count: %d\n\n", Nspline);
	fprintf(stdout, "Loading/precomputing derivative overlaps... ");
	fflush(stdout);
	
	// Precompute matrix of derivative overlaps ---------------------------- //
	//
	CooMatrix D (Nspline,Nspline);
	//
	if (not D.hdfload("D.hdf"))
		D.symm_populate_band(
			order,	// band halfwidth
			[ = ](unsigned i, unsigned j) -> Complex {
				return computeD(i, j, Nknot - 1);
			}
		).hdfsave("D.hdf");
	// --------------------------------------------------------------------- //
	
	if (iproc == 0)
	{
		fprintf(stdout, "ok\nLoading/precomputing integral moments... ");
		fflush(stdout);
	}
	
	// Precompute useful integral moments ---------------------------------- //
	//
	CooMatrix S(Nspline, Nspline);
	CooMatrix Mm1(Nspline, Nspline), Mm1_tr(Nspline, Nspline);
	CooMatrix Mm2(Nspline, Nspline);
	//
	if (not S.hdfload("S.hdf"))
		S.symm_populate_band(
			order,	// band halfwidth
			[ = ](unsigned m, unsigned n) -> Complex {
				return computeM(0, m, n);
			}
		).hdfsave("S.hdf");
	//
	if (not Mm1.hdfload("Mm1.hdf"))
		Mm1.symm_populate_band(
			order,	// band halfwidth
			[ = ](unsigned m, unsigned n) -> Complex {
				return computeM(-1, m, n);
			}
		).hdfsave("Mm1.hdf");
	//
	if (not Mm1_tr.hdfload("Mm1_tr.hdf"))
		Mm1_tr.symm_populate_band(
			order,	// band halfwidth
			[ = ](unsigned m, unsigned n) -> Complex {
				return computeM(-1, m, n, Nreknot - 1);
			}
		).hdfsave("Mm1_tr.hdf");
	//
	if (not Mm2.hdfload("Mm2.hdf"))
		Mm2.symm_populate_band(
			order,	// band halfwidth
			[ = ](unsigned m, unsigned n) -> Complex {
				return computeM(-2, m, n);
			}
		).hdfsave("Mm2.hdf");
	// --------------------------------------------------------------------- //
	
	fprintf(stdout, "ok\n");
	fflush(stdout);
	
	// Precompute two-electron integrals ----------------------------------- //
	//
	// R-type integral of a given multipole moment may be stored in
	// Nspline × Nspline × Nspline × Nspline four-dimensional field.
	// Only Nspline × (2*order+1) × Nspline × (2*order+1) of these elements
	// are nonzero. This four-dimensional array can be compacted
	// to two-dimensional by grouping pairs of indices to two multi-indices,
	// (i,j,k,l) -> ([ij], [kl]), each of these multi-indices having values
	// from the range [0 .. Nspline × (2*order+1)]. The matrix is then symmetric
	// with respect to exchange of the multi-index [ij] and [kl], with respect
	// to the interchange of indices in multiindex and it is dense.
	
	int maxlambda = 2 * maxell;
	std::vector<CooMatrix> R_tr(maxlambda + 1);
	
	//
	//  a) (ijv)-representation of sparse amtrices
	//
	
	std::vector<long>   Rtr_i, Rtr_j;
	std::vector<Complex> Rtr_v;
	
	//
	//  b) iterate over multipoles
	//
	
	#pragma omp parallel
	{
		#pragma omp master
		fprintf(
			stdout,
			"Precomputing multipole integrals (λ = 0 .. %d) using %d threads.\n",
			maxlambda, omp_get_num_threads()
		); fflush(stdout);
	}
	
	for (int lambda = 0; lambda <= maxlambda; lambda++)
	{
		//
		// c) look for precomputed data on disk
		//
		
		std::ostringstream oss2;
		oss2 << "R_tr[" << lambda << "].hdf";
		if (R_tr[lambda].hdfload(oss2.str().c_str()))
		{
			fprintf(stdout, "\t- integrals for λ = %d successfully loaded\n", lambda);
			fflush(stdout);
			continue; // no need to compute
		}
		
		//
		// d) precompute necessary partial integral moments
		// 
		
		// truncated moments λ and -λ-1
		cArray Mtr_L    = computeMi( lambda,   Nreknot-1);
		cArray Mtr_mLm1 = computeMi(-lambda-1, Nreknot-1);
		
		//
		// e) compute elements of Rtr
		//
		
		Rtr_i.clear();
		Rtr_j.clear();
		Rtr_v.clear();
		
		#pragma omp parallel \
			default (none) \
			shared (stdout, Rtr_i, Rtr_j, Rtr_v) \
			firstprivate (Nspline, order, Mtr_L, Mtr_mLm1, lambda, iproc)
		{
			// reserve threads' private storage
			int vol = Nspline * Nspline * order * order;
			std::vector<long>   th_Rtr_i, th_Rtr_j;
			std::vector<Complex> th_Rtr_v;
			th_Rtr_i.reserve(vol);
			th_Rtr_j.reserve(vol);
			th_Rtr_v.reserve(vol);
			
			// [ij] multi-index loop
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < Nspline; i++)
			{
				
				#pragma omp critical
				fprintf(stdout, "\r\t- multipole λ = %d... %3.0f %%", lambda, i * 100. / Nspline);
				fflush(stdout);
				
				for (int j = i; j < Nspline; j++)
				{
					// [ij] and [ji] multi-indices
					int ij_multi = i * Nspline + j;
					int ji_multi = j * Nspline + i;
					
					// [kl] multi-index loop; both "k" and "l" are restricted
					// to relatively narrow band of width 2*order+1
					for (int k = ((i > order)? i - order : 0); k < Nspline and k <= i + order; k++)
					{
						for (int l = ((j > order)? j - order : 0); l < Nspline and l <= j + order; l++)
						{
							// [kl] and [lk] multi-indices
							int kl_multi = k * Nspline + l;
							int lk_multi = l * Nspline + k;
							
							// compute the integral
							Complex Rijkl_tr = computeR(lambda, i, j, k, l, Mtr_L, Mtr_mLm1);
							
							// coordinates of nonzero-elements
							th_Rtr_i.push_back(ij_multi);
							th_Rtr_j.push_back(kl_multi);
							
							// R₀-truncated integral
							th_Rtr_v.push_back(Rijkl_tr.real());
							
							// if not at ij-diagonal, use symmetry
							//  Rλ(i,j,k,l) = Rλ(j,i,l,k)
							if (i != j)
							{
								// coordinates of nonzero-elements
								th_Rtr_i.push_back(ji_multi);
								th_Rtr_j.push_back(lk_multi);
								
								// R₀-truncated integral
								th_Rtr_v.push_back(Rijkl_tr.real());
							}
						}
					} // end of [jl] multi-index loop
					
				}
			} // end of [ik] multi-index loop
			
			// copy threads' results to shared arrays
			#pragma omp critical
			{
				Rtr_i.insert(Rtr_i.begin(), th_Rtr_i.begin(), th_Rtr_i.end());
				Rtr_j.insert(Rtr_j.begin(), th_Rtr_j.begin(), th_Rtr_j.end());
				Rtr_v.insert(Rtr_v.begin(), th_Rtr_v.begin(), th_Rtr_v.end());
			}
			
		} // end of PARALLEL
		
		//
		// f) create coo-matrices
		//
		
		size_t R_size = Nspline * Nspline;
		R_tr[lambda] = CooMatrix (R_size, R_size, Rtr_i, Rtr_j, Rtr_v).shake();
		
		//
		// g) save them to disk
		//
		
		R_tr[lambda].hdfsave(oss2.str().c_str());
		
		fprintf(stdout, "\r\t- multipole λ = %d... ok            \n", lambda);
		fflush(stdout);
		
	}
	
	fprintf(stdout, "\t- R_tr[λ] has %ld nonzero elements\n", R_tr[0].size());
	fflush(stdout);
	// --------------------------------------------------------------------- //
	
	
	unsigned Hsize = Nspline * Nspline * (maxell + 1) * (maxell + 1);
	fprintf(stdout, "Hamiltonian properties:\n");
	fprintf(stdout, "\t-> hamiltonian size: %d\n", Hsize);
	fprintf(stdout, "Loading/computing B-spline expansions... ");
	fflush(stdout);
	
	
	//
	// expansion weights
	//
	auto weight_edge_damp = [ R0 ](Complex z) -> double {
		// this will suppress function value from R0+1 onwards
		// which is useful for expanding (divergent) Ricatti-Bessel function
		return (z.imag() == 0.) ? (1+tanh(R0 - 5 - z.real()))/2 : 0.;
	};
	auto weight_end_damp = [ Rmax ](Complex z) -> double {
		// whis will suppress function value at Rmax
		// which is useful for expanding everywhere-nonzero hydrogenic function
//		return (z.imag() == 0.) ? tanh(Rmax - z.real()) : 0.;
		return tanh(Rmax - z.real());
	};
	
	// Prepare B-spline overlaps and expansions of Ric-Bess functions ------- //
	//
	//  j-overlaps of shape [Nenergy × Nangmom × Nspline]
	cArray ji_overlaps;
	if (not ji_overlaps.hdfload("ji_overlaps_damp.hdf"))
	{
		ji_overlaps = overlapj(maxell,ki,weight_edge_damp);
		ji_overlaps.hdfsave("ji_overlaps_damp.hdf");
	}
	// 
	//  compute expansions; solve the system
	//      S * B_spline_expansion = B_spline_overlap
	//
	unsigned ji_expansion_count = ji_overlaps.size()/Nspline;
	cArray ji_expansion;
 	if (not ji_expansion.hdfload("ji_expansion.hdf"))
	{
		ji_expansion = S.solve(ji_overlaps, ji_expansion_count);
		ji_expansion.hdfsave("ji_expansion.hdf");
	}
	//
	// construct sparse matrices from the data
	//
	std::vector<CooMatrix> ji_coo(ji_overlaps.size()/Nspline);
	for (unsigned idx_j = 0; idx_j < ji_expansion_count; idx_j++)
		ji_coo[idx_j] = CooMatrix(Nspline, 1, ji_expansion.begin() + idx_j * Nspline);
	// --------------------------------------------------------------------- //
	
	fprintf(stdout, "ok\n");
	fflush(stdout);
	
	// Precompute some accelerators ---------------------------------------- //
	//
	// Kronecker producs
	fprintf(stdout, "Creating Kronecker products... "); fflush(stdout);
	CsrMatrix S_kron_S   = kron(S, S).tocsr();
	CooMatrix S_kron_Mm1 = kron(S, Mm1_tr);
	CsrMatrix S_kron_Mm2 = kron(S, Mm2).tocsr().sparse_like(S_kron_S);
	CooMatrix Mm1_kron_S = kron(Mm1_tr, S);
	CsrMatrix Mm2_kron_S = kron(Mm2, S).tocsr().sparse_like(S_kron_S);
	CsrMatrix half_D_minus_Mm1_tr_kron_S = kron(0.5 * D - Mm1_tr, S).tocsr().sparse_like(S_kron_S);
	CsrMatrix S_kron_half_D_minus_Mm1_tr = kron(S, 0.5 * D - Mm1_tr).tocsr().sparse_like(S_kron_S);
	
	// CSR representation of R_tr matrices
	std::vector<CsrMatrix> R_tr_csr(maxlambda + 1);
	std::transform(
		R_tr.begin(),
		R_tr.end(),
		R_tr_csr.begin(),
		[ & ](const CooMatrix& coo) -> CsrMatrix { return coo.tocsr().sparse_like(S_kron_S); }
	);
	fprintf(stdout, "ok\n");

	// For all right hand sides -------------------------------------------- //
	//
	unsigned iterations_done = 0;
	for (unsigned ie = 0; ie < Nenergy; ie++)
	{
		// print progress information
		fprintf(stdout,
			"\rSolving the systems... %3.0f %% (typically %4d CG iterations per energy)  ",
			ie * 100. / Nenergy,
			ie == 0 ? 0 : iterations_done / ie
		); fflush(stdout);
		
		cArray current_solution, previous_solution;
		
		// we may have already computed all solutions for this energy... is it so?
		bool all_done = true;
		for (int Spin = 0; Spin <= 1; Spin++)
		for (int li = minli; li <= maxli; li++)
		for (int mi = -li; mi <= li; mi++)
		{
			// compose filename of the output file for this solution
			std::ostringstream oss;
			oss << "psi-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
		
			// check if there is some precomputed solution on the disk
			if ( not current_solution.hdfload(oss.str().c_str()) )
				all_done = false;
		}
		if (all_done) continue;
		
		// get total energy of the system
		double E = 0.5 * (Ei[ie] - 1./(ni*ni));
		
		// setup preconditioner - the diagonal block LU-factorizations
		std::vector<CsrMatrix> blocks((maxell + 1)*(maxell + 1));
		std::vector<CsrMatrix::LUft> lufts((maxell + 1)*(maxell + 1));
		
		// for each diagonal block assigned to this node
		int l1 = iproc;
//  		# pragma omp parallel for schedule (dynamic,1)
		for (int l2 = 0; l2 <= maxell; l2++)
		{
			int iblock = l1 * (maxell + 1) + l2;
			
			fprintf(stdout, "[%d] Factorization %d start\n", iproc, iblock); fflush(stdout);
			
			// one-electron parts
			CsrMatrix Hdiag;
			Hdiag = half_D_minus_Mm1_tr_kron_S;
			Hdiag &= ((0.5*l1*(l1+1)) * Mm2_kron_S);
			Hdiag &= S_kron_half_D_minus_Mm1_tr;
			Hdiag &= ((0.5*l2*(l2+1)) * S_kron_Mm2);
			
			// two-electron part
			for (int lambda = 0; lambda <= maxlambda; lambda++)
			{
				Complex _f = computef(lambda,l1,l2,l1,l2,L);
				if (_f != 0.)
					Hdiag &= _f * R_tr_csr[lambda];
			}
			
			// finalize the matrix
			blocks[iblock] = (E*S_kron_S) ^ Hdiag;
		
			// compute the LU factorization
			lufts[iblock] = blocks[iblock].factorize();
			
			fprintf(stdout, "[%d] Factorization %d stop\n", iproc, iblock); fflush(stdout);
		}
		
		// For all initial states ------------------------------------------- //
		//
		
		for (int li = minli; li <= maxli; li++)
		{
			// compute P-overlaps and P-expansion
			cArray Pi_overlaps, Pi_expansion;
			Pi_overlaps = overlapP(ni,li,weight_end_damp);
			Pi_overlaps.hdfsave("Pi_overlaps_full.hdf");
			Pi_expansion = S.solve(Pi_overlaps);
			Pi_expansion.hdfsave("Pi_expansion.hdf");
			CooMatrix Pi_coo(Nspline, 1, Pi_expansion.begin());
			
			for (int Spin = 0; Spin <= 1; Spin++)
			for (int mi = -li; mi <= li; mi++)
			{
				// we may have already computed solution for this state and energy... is it so?
				std::ostringstream cur_oss;
				cur_oss << "psi-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
				if ( current_solution.hdfload(cur_oss.str().c_str()) )
					continue;
				
				// create right hand side
				cArray chi ( (maxell+1)*(maxell+1)*Nspline*Nspline );
				
				// CG preconditioner callback
				auto apply_preconditioner = [ & ](const cArray& r, cArray& z) -> void {
					
					// apply a block inversion preconditioner
					
					cArray z_node(Nspline * Nspline * (maxell + 1));
					
					# pragma omp parallel for schedule (dynamic,1)
					for (int l2 = 0; l2 <= maxell; l2++)
					{
						int iblock = l1 * (maxell + 1) + l2;
						
						// create a copy of a RHS segment
						cArray r_block(Nspline * Nspline);
						memcpy(
							&r_block[0],						// dest
							&r[0] + iblock * Nspline * Nspline,	// src
							Nspline * Nspline * sizeof(Complex)	// n
						);
						
						// multiply by an inverted block
						cArray z_block;
						z_block = lufts[iblock].solve(r_block);
						memcpy(
							&z_node[0] + l2 * Nspline * Nspline,	// dest
							&z_block[0],							// src
							Nspline * Nspline * sizeof(Complex)		// n
						);
					}
					
					// assemble other work from other processes
					fprintf(stdout, "[%d] Waiting for PRECOND MPI_Gather... \n", iproc); fflush(stdout);
					MPI::COMM_WORLD.Allgather(
						&z_node[0], Nspline * Nspline * (maxell + 1), MPI::DOUBLE_COMPLEX,
						&z[0], Nspline * Nspline * (maxell + 1), MPI::DOUBLE_COMPLEX
					);
					fprintf(stdout, "[%d] PRECOND MPI_Gather done\n", iproc); fflush(stdout);
				};
				
				// CG matrix multiplication callback
				auto matrix_multiply = [ & ](const cArray& p, cArray& q) -> void {
					
					// multiply by the matrix of the system
					
					// product
					cArray q_node(Nspline * Nspline * (maxell + 1));
					
					# pragma omp parallel for schedule (dynamic,1)
					for (int l2 = 0; l2 <= maxell; l2++)
					{
						int iblock = l1 * (maxell + 1) + l2;
						
						// product
						cArray q_block(Nspline * Nspline);
						
						// multiply block-row of the matrix with "p"
						for (int l1p = 0; l1p <= maxell; l1p++)
						for (int l2p = 0; l2p <= maxell; l2p++)
						{
							int iblockp = l1p * (maxell + 1) + l2p;
							
							// corresponding fragment of "p"
							cArray p_block(Nspline * Nspline);
							memcpy(
								&p_block[0],							// dest
								&p[0] + iblockp * Nspline * Nspline,	// src
								Nspline * Nspline * sizeof(Complex)		// n
							);
							
							// multiply by hamiltonian terms
							if (iblock == iblockp)
							{
								// reuse the diagonal block
								q_block += blocks[iblock].dot(p_block);
							}
							else
							{
								// compute the offdiagonal block
								for (int lambda = 0; lambda <= maxlambda; lambda++)
								{
									Complex _f = computef(lambda, l1, l2, l1p, l2p, L);
									if (_f != 0.)
										q_block -= _f * R_tr_csr[lambda].dot(p_block);
								}
							}
						}
						
						// copy product to output
						memcpy(
							&q_node[0] + l2 * Nspline * Nspline,	// dest
							&q_block[0],							// src
							Nspline * Nspline * sizeof(Complex)		// n
						);
					}
					
					// assemble other work from other processes
// 					fprintf(stdout, "Waiting for MMUL MPI_Gather..."); fflush(stdout);
					MPI::COMM_WORLD.Allgather(
						&q_node[0], Nspline * Nspline * (maxell + 1), MPI::DOUBLE_COMPLEX,
						&q[0], Nspline * Nspline * (maxell + 1), MPI::DOUBLE_COMPLEX
					);
// 					fprintf(stdout, "ok\n"); fflush(stdout);
				};
				
				// for all segments constituting the RHS
				# pragma omp parallel for collapse(2)
				for (int l1 = 0; l1 <= maxell; l1++)
				for (int l2 = 0; l2 <= maxell; l2++)
				{
					int iblock = l1 * (maxell + 1) + l2;
					cArray chi_block(Nspline * Nspline);
					
					// for all allowed angular momenta (by momentum composition) of the projectile
					for (int l = abs(li - L); l <= li + L; l++)
					{
						cArray Pj1, Pj2;
						const CooMatrix& ji_coo_E_l = ji_coo[ie * (maxell + 1) + l];
						Pj1 = kron(Pi_coo, ji_coo_E_l).todense();
						Pj2 = kron(ji_coo_E_l, Pi_coo).todense();
						
						// (anti)symmetrization
						Complex Sign = ((Spin + L + li + l) % 2 == 0) ? 1. : -1.;
						
						// compute energy- and angular momentum-dependent prefactor
						Complex prefactor = pow(Complex(0.,1.),l) * sqrt(2*M_PI*(2*l+1)) / Complex(ki[ie]); 
						prefactor *= ClebschGordan(li,mi,l,0,L,mi);
						
						for (int lambda = 0; lambda <= maxlambda; lambda++)
						{
							Complex _f1 = computef(lambda, l1, l2, li, l, L);
							Complex _f2 = computef(lambda, l1, l2, l, li, L);
							
							chi_block += prefactor * _f1 * R_tr[lambda].dot(Pj1).todense();
							chi_block += prefactor * _f2 * R_tr[lambda].dot(Pj2).todense() * Sign;
						}
						
						if (li == l1 and l == l2)
							chi_block -= prefactor * S_kron_Mm1.dot(Pj1).todense();
						
						if (li == l2 and l == l1)					
							chi_block -= prefactor * Mm1_kron_S.dot(Pj2).todense() * Sign;
					}
					
					// copy the block into the whole array
					memcpy(
						&chi[0] + iblock * Nspline * Nspline,	// dest
						&chi_block[0],							// src
						Nspline * Nspline * sizeof(Complex)		// n
					);
				}
				
				// Solve the Equations ------------------------------------------------- //
				//
				
				// we may have already computed the previous solution - it will serve as initial guess
				current_solution = cArray(chi.size());
				if (ie > 0)
				{
					std::ostringstream prev_oss;
					prev_oss << "psi-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie-1] << ".hdf";
					if (previous_solution.hdfload(prev_oss.str().c_str()))
						current_solution = previous_solution;
				}
				
				// custom conjugate gradients callback-based solver
				fprintf(stdout, "[%d] running callback\n", iproc); fflush(stdout);
				unsigned iterations = cg_callbacks(
					chi,					// right-hand side
					current_solution,		// on input, the initial guess, on return, the solution
					1e-5,					// requested precision, |A·x - b|² < ε·|b|²
					0,						// minimal iteration count
					(maxell+1) * Nspline,	// maximal iteration count
					apply_preconditioner,	// preconditioner callback
					matrix_multiply			// matrix multiply callback
				);
				fprintf(stdout, "[%d] finishing callback\n", iproc); fflush(stdout);
				fprintf(stdout, "."); fflush(stdout);
					
				// update progress
				iterations_done += iterations;
				
				// save solution to disk
				current_solution.hdfsave(cur_oss.str().c_str());
				
			} // end of For mi, Spin
		} // end of For li
		
		// free Numeric objects
		for (int l2 = 0; l2 <= maxell; l2++)
		{
			int iblock = l1 * (maxell + 1) + l2;
			lufts[iblock].free();
		}
		
	} // end of For ie = 0, ..., Nenergy - 1
	
	fprintf(stdout, "\rSolving the systems... ok                                                            \n");
	fprintf(stdout, "\t(typically %4d CG iterations per energy)\n", iterations_done / Nenergy);
	fflush(stdout);
	
	// --------------------------------------------------------------------- //
	
	// Create SQL batch file
	char filename[100];
	sprintf(filename, "%d-%d-proc%d.sql", ni, L, iproc);
	FILE* fsql = fopen(filename, "w");
	fprintf(fsql, "BEGIN TRANSACTION;\n");
	
	// Extract the cross sections ------------------------------------------ //
	//
	
	fprintf(stdout, "\rExtracting T-matrices... %3.0f%%     ", 0.);
	fflush(stdout);
	
	// list of transitions
	std::vector<std::tuple<int,int,int,int,int>> transitions;
	for (int Spin = 0; Spin <= 1; Spin++)
	for (int li = minli; li <= maxli; li++)
	for (int mi = -li; mi <= li; mi++)
	for (int nf = 1; nf <= maxnf; nf++)
	for (int lf = 0; lf <= maxlf and lf < nf; lf++)
		transitions.push_back(std::make_tuple(Spin,li,mi,nf,lf));
	int Ntransitions = transitions.size();
	
	// finished transitions
	int finished = 0;
	
	// max transitions for this process
	int max_finished = Ntransitions % Nproc > iproc ?  Ntransitions / Nproc + 1 : Ntransitions / Nproc;
	
	// for all transitions
	for (int i = 0; i < Ntransitions; i++)
	{
		if (i % Nproc != iproc)
			// some other process will compute this transition
			continue;
		
		int Spin = std::get<0>(transitions[i]);
		int li   = std::get<1>(transitions[i]);
		int mi   = std::get<2>(transitions[i]);
		int nf   = std::get<3>(transitions[i]);
		int lf   = std::get<4>(transitions[i]);
		
		cArray Pf_overlaps = overlapP(nf,lf,weight_end_damp);
		
		// compute radial integrals
		std::vector<cArray> Lambda(2 * lf + 1);
		for (int mf = -lf; mf <= lf; mf++)
		{
			// compute expansions of final wave function for this final state (nf,lf)
			std::transform(
				Ei.begin(), Ei.end(),
				kf.begin(),
				[ = ](double E) -> double { return sqrt(E - 1./(ni*ni) + 1./(nf*nf) + (mf-mi)*B); }
			);
			cArray jf_overlaps = overlapj(maxell,kf,weight_edge_damp);
			
			// compute Λ for transitions to (nf,lf,mf); it will depend on [ie,ℓ]
			Lambda[mf+lf] = computeLambda(kf, ki, maxell, L, Spin, ni, li, mi, Ei, lf, Pf_overlaps, jf_overlaps);
		}
		
		// save the data
		for (int mf = -lf; mf <= lf; mf++)
		{
			// compute Tℓ
			cArray T_ell(Lambda[mf+lf].size());
			for (unsigned i = 0; i < T_ell.size(); i++)
			{
				int ie  = i / (maxell + 1);
				int ell = i % (maxell + 1);
				
				T_ell[i] = Lambda[mf+lf][i] * 4. * M_PI / kf[ie] * pow(Complex(0.,1.), -ell)
								* ClebschGordan(lf, mf, ell, mi - mf, L, mi) / sqrt(2.);
			}
			
			//
			// print out SQL
			//
			
			for (unsigned i = 0; i < T_ell.size(); i++)
			{
				int ie  = i / (maxell + 1);
				int ell = i % (maxell + 1);
				
				if (finite(T_ell[i].real()) and finite(T_ell[i].imag()))
				if (T_ell[i].real() != 0. or T_ell[i].imag() != 0.)
				{
					fprintf(fsql, "INSERT OR REPLACE INTO \"hex\" VALUES (%d,%d,%d, %d,%d,%d, %d,%d, %e, %d, %e,%e);\n",
						ni, li, mi, nf, lf, mf, L, Spin, Ei[ie], ell, T_ell[i].real(), T_ell[i].imag()
					);
				}
			}
			
			
			//
			// print out the total cross section for quick overview
			//
			
			sprintf(filename, "sigma-%d-%d-%d-%d-%d-%d-%d-%d.dat", ni, li, mi, nf, lf, mf, L, Spin);
			FILE* ftxt = fopen(filename, "w");
			fprintf(ftxt, "# Ei [Ry] sigma [a0^2]\n");
			for (unsigned ie = 0; ie < Nenergy; ie++)
			{
				double sigma = 0.;
				for (int ell = 0; ell <= maxell; ell++)
				{
					double Re_f_ell = -T_ell[ie * (maxell + 1) + ell].real() / (2 * M_PI);
					double Im_f_ell = -T_ell[ie * (maxell + 1) + ell].imag() / (2 * M_PI);
					sigma += 0.25 * (2*Spin + 1) * kf[ie] / ki[ie] * (Re_f_ell * Re_f_ell + Im_f_ell * Im_f_ell);
				}
				if (finite(sigma))
				{
					fprintf(ftxt, "%g\t%g\n", Ei[ie], sigma);
				}
			}
			fclose(ftxt);
			
		}
		
		finished++;
		
		fprintf(stdout, "\rExtracting T-matrices... %3.0f%%     ", finished * 100. / max_finished);
		fflush(stdout);
	}
	fprintf(stdout, "\rExtracting T-matrices... ok       \n"); fflush(stdout);
	// --------------------------------------------------------------------- //
	
	fprintf(fsql, "COMMIT;\n");
	fclose(fsql);
	
	fprintf(stdout, "\nDone.\n\n");
	fflush(stdout);
	
	fclose(stdout);
	
	MPI::Finalize();
	
	return 0;
}
