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
#include <cstdio>
#include <cmath>
#include <chrono>
#include <cstring>
#include <vector>
#include <iostream>
#include <tuple>

#include <gsl/gsl_errno.h>
#include <omp.h>
#include <H5Cpp.h>

#ifndef NO_MPI
	#include <mpi.h>
#endif

#include "amplitudes.h"
#include "angs.h"
#include "arrays.h"
#include "bspline.h"
#include "complex.h"
#include "input.h"
#include "misc.h"
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
	// Preparations ------------------------------------------------------- //
	//
	
	gsl_set_error_handler_off();
    
    // disable buffering of the standard output 
    // (so that all text messages are immediatelly visible)
    std::setvbuf(stdout, nullptr, _IONBF, 0);
	
	// variables that can be set by the user from the command line
	std::ifstream inputfile;			// input file
	std::string zipfile;				// HDF solution expansion file to zip
	int  zipcount = 0;					// zip sample count
	bool parallel = false;				// whether to use OpenMPI
	bool stg1 = false;					// whether to computate radial integrals only
	bool stg12 = false;					// whether to computate radial integrals and solutions only
	
	// get input from command line
	parse_command_line (
		argc, argv,
		inputfile,
		zipfile, zipcount,
		parallel,
		stg1, stg12
	);
	
	// setup MPI
	int Nproc = 1, iproc = 0;
	if (parallel)
	{
#ifndef NO_MPI
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &Nproc);
		MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
#else
		std::cout << "WARNING: --mpi has no effect, as the program was build without the MPI support!\n";
#endif
	}
	bool I_am_master = (iproc == 0);
	
	// check input file
	if (not inputfile.is_open())
	{
		inputfile.open("hex.inp");
		if (inputfile.bad())
			throw exception("Input error: Cannot open the file \"hex.inp\".\n");
	}
	
	// create main variables
	int order;							// B-spline order
	double R0, Rmax, ecstheta;			// B-spline knot sequence
	rArray real_knots;					// real knot sequence
	rArray complex_knots;				// complex knot sequence
	int ni;								// atomic state parameters
	int L;								// global conserved variables
	int maxell;							// angular momentum limits
	rArray Ei;							// energies in Rydbergs
	double B = 0;						// magnetic field in a.u.
	
	// initial and final atomic states
	std::vector<std::tuple<int,int,int>> instates, outstates;
	
	// get input from input file
	parse_input_file (
		inputfile, order, R0, ecstheta, Rmax,
		real_knots, complex_knots, ni, 
		instates, outstates,
		L, maxell, Ei, B
	);
	
	// is there something to compute?
	if (Ei.empty() or instates.empty() or outstates.empty())
	{
		std::cout << "Nothing to compute.\n";
		exit(0);
	}
	
	// shorthand for total number of energies
	unsigned Nenergy = Ei.size();
	
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
	Bspline::ECS().init(order, real_knots, R0, ecstheta, complex_knots, Rmax);
	
	// shortcuts
	int Nspline = Bspline::ECS().Nspline();
	int Nknot = Bspline::ECS().Nknot();
	int Nreknot = Bspline::ECS().Nreknot();
	
	// info
	std::cout << "B-spline total count: " << Bspline().Nspline() << "\n\n";
	//
	// --------------------------------------------------------------------- //
	
	
	if (zipfile.size() != 0 and I_am_master)
	{
		cArray sol;			// stored solution expansion
		cArray ev;			// evaluated solution
		rArray grid;		// real evaluation grid
		
		std::cout << "Zipping B-spline expansion of the solution: \"" << zipfile << "\"" << std::endl;
		
		if (not sol.hdfload(zipfile.c_str()))
			throw exception("Cannot load file %s.", zipfile.c_str());
		
		grid = linspace(0., Rmax, zipcount + 1);
		
		for (int l1 = 0; l1 <= maxell; l1++)
		for (int l2 = 0; l2 <= maxell; l2++)
		{
			// skip zero segments
			if (std::abs(l1 - l2) > L or l1 + l2 < L)
				continue;
			
			std::cout << "\t- partial wave l1 = " << l1 << ", l2 = " << l2 << "\n";
			
			// zip this partial wave
			ev = Bspline::ECS().zip (
				cArrayView (
					sol,
					(l1 * (maxell + 1) + l2) * Nspline * Nspline,
					Nspline * Nspline
				),
				grid,
				grid
			);
			
			// setup output filename
			char outf1[3 + zipfile.size()], outf2[3 + zipfile.size()];
			std::sprintf(outf1, "%s-(%d,%d).re", zipfile.c_str(), l1, l2);
			std::sprintf(outf2, "%s-(%d,%d).im", zipfile.c_str(), l1, l2);
			
			// write real part
			write_2D_data (
				zipcount + 1,
				zipcount + 1,
				outf1,
				[&](size_t i, size_t j) -> double {
					return ev[i * (zipcount + 1) + j].real();
				}
			);
			
			// write imaginary part
			write_2D_data (
				zipcount + 1,
				zipcount + 1,
				outf2,
				[&](size_t i, size_t j) -> double {
					return ev[i * (zipcount + 1) + j].imag();
				}
			);
		}
		
		exit(0);
	}
	
	std::cout << "Loading/precomputing derivative overlaps... ";
	
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
	
	
	std::cout << "ok\nLoading/precomputing integral moments... ";
	
	
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
	
	
	std::cout << "ok\n";
	
	
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
		{
			std::cout << "Precomputing multipole integrals (λ = 0 .. " 
			          << maxlambda 
			          << ") using " 
					  << omp_get_num_threads() 
					  << " threads.\n";
		}
	}
	
	for (int lambda = 0; lambda <= maxlambda; lambda++)
	{
		//
		// c) look for precomputed data on disk
		//
		
		std::ostringstream oss2;
		oss2 << "R_tr[" << lambda << "].hdf";
		
		int R_integ_exists; // whether there are valid precomputed data
		
		// master will load radial integrals...
		//   We could let all processess load the data on their own, but if the
		// task is executed on a cluster with remote storage and only symbolic
		// links to precomputed radial integral HDFs are created, every process
		// would download its copy of R_tr[*].hdf files from the remote storge.
		// That would greatly raise the traffic and also delay the start of
		// computation. If only the master process loads the data and then
		// distributes them using local network, everything is smoother and
		// faster ;-)
#ifndef NO_MPI
		if (I_am_master)
		{
#endif
			R_integ_exists = R_tr[lambda].hdfload(oss2.str().c_str());
#ifndef NO_MPI
		}
#endif
		
#ifndef NO_MPI
		if (parallel)
		{
			// master will broadcast existence information to other processes
			MPI_Bcast(&R_integ_exists, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
			// master will broadcast data to other processes
			if (R_integ_exists)
			{
				// get dimensions
				int size = R_tr[lambda].size();
				int m = R_tr[lambda].rows();
				int n = R_tr[lambda].cols();
				
				// master will broadcast dimensions
				MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&m,    1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&n,    1, MPI_INT, 0, MPI_COMM_WORLD);
				
				// get arrays
				std::vector<long>    i = R_tr[lambda].i();    i.resize(size);
				std::vector<long>    j = R_tr[lambda].j();    j.resize(size);
				std::vector<Complex> v = R_tr[lambda].v();    v.resize(size);
				
				// master will broadcast arrays
				MPI_Bcast(&i[0], i.size(), MPI_LONG, 0, MPI_COMM_WORLD);
				MPI_Bcast(&j[0], j.size(), MPI_LONG, 0, MPI_COMM_WORLD);
				MPI_Bcast(&v[0], v.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
				
				// reconstruct objects
				R_tr[lambda] = CooMatrix(m, n, i, j, v);
			}
		}
#endif
		
		if (R_integ_exists)
		{
			std::cout << "\t- integrals for λ = " << lambda << " successfully loaded\n";
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
			firstprivate (Nspline, order, Mtr_L, Mtr_mLm1, lambda)
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
				{
					std::printf ("\r\t- multipole λ = %d... %3.0f %%", lambda, i * 100. / Nspline);
				}
				
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
		
		std::cout << "\r\t- multipole λ = " << lambda << "... ok            \n";
		
	}
	std::cout << "\t- R_tr[λ] has " << R_tr[0].size() << " nonzero elements\n";
	// --------------------------------------------------------------------- //
	
	
	unsigned Hsize = Nspline * Nspline * (maxell + 1) * (maxell + 1);
	std::cout << "Hamiltonian properties:\n";
	std::cout << "\t-> hamiltonian size: " << Hsize << "\n";
	
	// exit if reqested
	if (stg1)
	{
#ifndef NO_MPI
		if (parallel)
			MPI_Finalize();
#endif
		exit(0);
	}
	
	
	// Expansion weights --------------------------------------------------- //
	//
	auto weight_edge_damp = [ R0 ](Complex z) -> double {
		// this will suppress function value from R0+1 onwards
		// which is useful for expanding (divergent) Ricatti-Bessel function
		return (z.imag() == 0.) ? (1+tanh(R0 - 5 - z.real()))/2 : 0.;
	};
	auto weight_end_damp = [ Rmax ](Complex z) -> double {
		// whis will suppress function value at Rmax
		// which is useful for expanding everywhere-nonzero hydrogenic function
		return tanh(Rmax - z.real());
	};
	// --------------------------------------------------------------------- //
	
	
	// Prepare B-spline overlaps and expansions of Ric-Bess functions ------ //
	//
	std::cout << "Computing B-spline expansions... ";
	
	//  j-overlaps of shape [Nenergy × Nangmom × Nspline]
	cArray ji_overlaps = overlapj(maxell,ki,weight_edge_damp);
	ji_overlaps.hdfsave("ji_overlaps_damp.hdf"); // just for debugging
	
	//  compute expansions; solve the system
	//      S * B_spline_expansion = B_spline_overlap
	unsigned ji_expansion_count = ji_overlaps.size()/Nspline;
	cArray ji_expansion = S.solve(ji_overlaps, ji_expansion_count);
	ji_expansion.hdfsave("ji_expansion.hdf"); // just for debugging
	
	// construct sparse matrices from the data
	std::vector<CooMatrix> ji_coo(ji_overlaps.size()/Nspline);
	for (unsigned idx_j = 0; idx_j < ji_expansion_count; idx_j++)
		ji_coo[idx_j] = CooMatrix(Nspline, 1, ji_expansion.begin() + idx_j * Nspline);
	
	std::cout << "ok\n";
	// --------------------------------------------------------------------- //
	
	
	// Precompute some accelerators ---------------------------------------- //
	//
	// Kronecker producs
	std::cout << "Creating Kronecker products... ";
	CsrMatrix S_kron_S   = kron(S, S).tocsr();
	CooMatrix S_kron_Mm1_tr = kron(S, Mm1_tr);
	CsrMatrix S_kron_Mm2 = kron(S, Mm2).tocsr().sparse_like(S_kron_S);
	CooMatrix Mm1_tr_kron_S = kron(Mm1_tr, S);
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
	std::cout << "ok\n";
	// --------------------------------------------------------------------- //
	
	
	// Distribute LU factorizations among processes ------------------------ //
	//
	std::map<int,int> LUs;
	int worker = 0;
	std::cout << "Balancing " << triangle_count(L,maxell)
	          << " diagonal blocks among " << Nproc 
	          << " worker processes...";
	for (int l1 = 0; l1 <= maxell; l1++)
	for (int l2 = 0; l2 <= maxell; l2++)
	{
		// skip those angular momentum pairs that don't compose L
		if (abs(l1 - l2) > L or l1 + l2 < L)
				continue;
		
		int iblock = l1 * (maxell + 1) + l2;  // get block index
		LUs[iblock] = worker;	              // assign this block to worker
		worker = (worker + 1) % Nproc;        // move to next worker
	}
	std::cout << "ok\n";
	// --------------------------------------------------------------------- //
	
	
	// For all right hand sides -------------------------------------------- //
	//
	int iterations_done = 0, computations_done = 0;
	for (unsigned ie = 0; ie < Nenergy; ie++)
	{
		// print progress information
		std::printf (
			"\nSolving the system for Ei[%d] = %g (%3.0f %% finished, typically %4d CG iterations per energy)\n",
			ie,
		    Ei[ie],
			ie * 100. / Nenergy,
			computations_done == 0 ? 0 : iterations_done / computations_done
		);
		
		cArray current_solution, previous_solution;
		
		// we may have already computed all solutions for this energy... is it so?
		bool all_done = true;
		for (int Spin = 0; Spin <= 1; Spin++)
		for (auto instate : instates)
		{
			int li = std::get<1>(instate);
			int mi = std::get<2>(instate);
			
			// skip angular forbidden states
			bool allowed = false;
			for (int l = abs(li - L); l <= li + L; l++)
				allowed = allowed or ClebschGordan(li,mi,l,0,L,mi);
			if (not allowed)
				continue;
			
			// compose filename of the output file for this solution
			std::ostringstream oss;
			oss << "psi-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
		
			// check if there is some precomputed solution on the disk
			if ( not current_solution.hdfload(oss.str().c_str()) )
				all_done = false;
		}
		if (all_done)
		{
			std::cout << "\tAll solutions for Ei[" << ie << "] = " << Ei[ie] << " loaded.\n";
			continue;
		}
		
		// get total energy of the system
		double E = 0.5 * (Ei[ie] - 1./(ni*ni));
		
		// setup preconditioner - the diagonal block LU-factorizations
		std::vector<CsrMatrix> blocks((maxell+1)*(maxell+1));
		std::vector<CsrMatrix::LUft> lufts((maxell+1)*(maxell+1));
		# pragma omp parallel for collapse (2) schedule (dynamic,1)
		for (int l1 = 0; l1 <= maxell; l1++)
		for (int l2 = 0; l2 <= maxell; l2++)
		{
			// get diagonal block index
			int iblock = l1 * (maxell + 1) + l2;
			
			// skip computation of unwanted blocks for this process
			if (LUs.find(iblock) == LUs.end() or LUs[iblock] != iproc)
				continue;
			
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
		}
		
		// compute the LU factorizations
		for (int l1 = 0; l1 <= maxell; l1++)
		for (int l2 = 0; l2 <= maxell; l2++)
		{
			// get diagonal block index
			int iblock = l1 * (maxell + 1) + l2;
			
			// skip computation of unwanted blocks for this process
			if (LUs.find(iblock) == LUs.end() or LUs[iblock] != iproc)
				continue;
			
			// timer info
			std::chrono::steady_clock::time_point start;
			std::chrono::duration<int> sec;
			
			// log output
			std::cout << "\t[" << iproc << "] LU factorization " 
			          << iblock << " of (" << l1 << "," << l2 << ") block started\n";
			start = std::chrono::steady_clock::now();
			
			// factorize
			lufts[iblock] = blocks[iblock].factorize();
			
			// log output
			sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
			std::cout << "\t[" << iproc << "] LU factorization " 
			          << iblock << " of (" << l1 << "," << l2 << ") block done after " 
					  << sec.count() << " s\n";
		}
		
		// For all initial states ------------------------------------------- //
		//
		
		for (auto instate : instates)
		{
			int li = std::get<1>(instate);
			int mi = std::get<2>(instate);
			
			// skip angular forbidden states
			bool allowed = false;
			for (int l = abs(li - L); l <= li + L; l++)
				allowed = allowed or ClebschGordan(li,mi,l,0,L,mi);
			if (not allowed)
				continue;
			
			// compute P-overlaps and P-expansion
			cArray Pi_overlaps, Pi_expansion;
			Pi_overlaps = overlapP(ni,li,weight_end_damp);
			Pi_expansion = S.solve(Pi_overlaps);
			CooMatrix Pi_coo(Nspline, 1, Pi_expansion.begin());
			
			for (int Spin = 0; Spin <= 1; Spin++)
			{
				// we may have already computed solution for this state and energy... is it so?
				std::ostringstream cur_oss;
				cur_oss << "psi-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
				if ( current_solution.hdfload(cur_oss.str().c_str()) )
					continue;
				
				// create right hand side
				cArray chi ( (maxell+1)*(maxell+1)*Nspline*Nspline );
				
				std::cout << "\tCreate RHS for li = " << li << ", mi = " << mi << ", S = " << Spin << "\n";
				
				// for all segments constituting the RHS
				# pragma omp parallel for collapse(2) schedule (dynamic,1)
				for (int l1 = 0; l1 <= maxell; l1++)
				for (int l2 = 0; l2 <= maxell; l2++)
				{
					// skip those angular momentum pairs that don't compose L
					if (abs(l1 - l2) > L or l1 + l2 < L)
						continue;
					
					// setup storage
					int iblock = l1 * (maxell + 1) + l2;
					cArrayView chi_block(chi, iblock * Nspline * Nspline, Nspline * Nspline);
					chi_block.clear();
					
					// for all allowed angular momenta (by momentum composition) of the projectile
					for (int l = abs(li - L); l <= li + L; l++)
					{
						// (anti)symmetrization
						int Sign = ((Spin + L + li + l) % 2 == 0) ? 1. : -1.;
						
						// compute energy- and angular momentum-dependent prefactor
						Complex prefactor = pow(Complex(0.,1.),l) * sqrt(2*M_PI*(2*l+1)) / Complex(ki[ie]); 
						prefactor *= ClebschGordan(li,mi,l,0,L,mi);
						if (prefactor == 0.)
							continue;
						
						cArray Pj1, Pj2;
						const CooMatrix& ji_coo_E_l = ji_coo[ie * (maxell + 1) + l];
						Pj1 = kron(Pi_coo, ji_coo_E_l).todense();
						Pj2 = kron(ji_coo_E_l, Pi_coo).todense();
						
						// skip angular forbidden right hand sides
						for (int lambda = 0; lambda <= maxlambda; lambda++)
						{
							Complex f1 = computef(lambda, l1, l2, li, l, L);
							Complex f2 = computef(lambda, l1, l2, l, li, L);
							
							if (f1 != 0.)
							{
								chi_block += (prefactor * f1) * R_tr[lambda].dot(Pj1).todense();
							}
							
							if (f2 != 0.)
							{
								if (Sign > 0)
									chi_block += (prefactor * f2) * R_tr[lambda].dot(Pj2).todense();
								else
									chi_block -= (prefactor * f2) * R_tr[lambda].dot(Pj2).todense();
							}
						}
						
						if (li == l1 and l == l2)
						{
							// direct contribution
							chi_block -= prefactor * S_kron_Mm1_tr.dot(Pj1).todense();
						}
						
						if (li == l2 and l == l1)
						{
							// exchange contribution
							if (Sign > 0)
								chi_block -= prefactor * Mm1_tr_kron_S.dot(Pj2).todense();
							else
								chi_block += prefactor * Mm1_tr_kron_S.dot(Pj2).todense();
						}
					}
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
				
				// CG preconditioner callback
				auto apply_preconditioner = [ & ](cArray const & r, cArray & z) -> void
				{
					
					// apply a block inversion preconditioner
					# pragma omp parallel for collapse (2) schedule (dynamic,1)
					for (int l1 = 0; l1 <= maxell; l1++)
					for (int l2 = 0; l2 <= maxell; l2++)
					{
						// get diagonal block index
						int iblock = l1 * (maxell + 1) + l2;
						
						// skip computation of unwanted blocks for this process
						if (LUs.find(iblock) == LUs.end() or LUs[iblock] != iproc)
							continue;
						
						// create copy-to view of "z"
						cArrayView zview(z, iblock * Nspline * Nspline, Nspline * Nspline);
						
						// create copy-from view of "r"
						cArrayView rview(r, iblock * Nspline * Nspline, Nspline * Nspline);
						
						// copy the correcponding slice of "r" multiplied by a correct block inversion
						zview = lufts[iblock].solve(rview);
					}
					
#ifndef NO_MPI
					if (parallel)
					{
						// synchronize across processes
						for (unsigned iblock = 0; iblock < z.size() / (Nspline * Nspline); iblock++)
						{
							// skip inactive segments
							if (LUs.find(iblock) == LUs.end())
								continue;
							
							// relevant process will broadcast this segment's data
							MPI_Bcast (
								&z[0] + iblock * Nspline * Nspline,
								Nspline * Nspline,
								MPI_DOUBLE_COMPLEX,
								LUs[iblock],
								MPI_COMM_WORLD
							);
						}
					}
#endif
				};
				
				// CG matrix multiplication callback
				auto matrix_multiply = [ & ](const cArray& p, cArray& q) -> void
				{
					
					// multiply by the matrix of the system
					# pragma omp parallel for collapse (2) schedule (dynamic,1)
					for (int l1 = 0; l1 <= maxell; l1++)
					for (int l2 = 0; l2 <= maxell; l2++)
					{
						// get diagonal block index
						int iblock = l1 * (maxell + 1) + l2;
						
						// skip computation of unwanted blocks for this process
						if (LUs.find(iblock) == LUs.end() or LUs[iblock] != iproc)
							continue;
						
						// product (copy-to view of "q")
						cArrayView q_block(q, iblock * Nspline * Nspline, Nspline * Nspline);
						q_block.clear(); // initialize with zeros
						
						// multiply block-row of the matrix with "p"
						for (int l1p = 0; l1p <= maxell; l1p++)
						for (int l2p = 0; l2p <= maxell; l2p++)
						{
							// skip those angular momentum pairs that don't compose L
							if (abs(l1p - l2p) > L or l1p + l2p < L)
								continue;
						
							// get diagonal block index
							int blockp = l1p * (maxell + 1) + l2p;
							
							// corresponding (copy-from) fragment of "p"
							cArrayView p_block(p, blockp * Nspline * Nspline, Nspline * Nspline);
							
							// multiply by hamiltonian terms
							if (iblock == blockp)
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
					}
					
#ifndef NO_MPI
					if (parallel)
					{
						// synchronize across processes
						for (unsigned iblock = 0; iblock < q.size() / (Nspline * Nspline); iblock++)
						{
							// skip inactive segments
							if (LUs.find(iblock) == LUs.end())
								continue;
							
							// relevant process will broadcast this segment's data
							MPI_Bcast (
								&q[0] + iblock * Nspline * Nspline,
								Nspline * Nspline,
								MPI_DOUBLE_COMPLEX,
								LUs[iblock],
								MPI_COMM_WORLD
							);
						}
					}
#endif
				};
				
				// custom conjugate gradients callback-based solver
				double tolerance = 1e-10;
				std::cout << "\tStart CG callback with tolerance " << tolerance << "\n";
				unsigned iterations = cg_callbacks (
					chi,					// right-hand side
					current_solution,		// on input, the initial guess, on return, the solution
					tolerance,				// requested precision, |A·x - b|² < ε·|b|²
					0,						// minimal iteration count
					(maxell+1) * Nspline,	// maximal iteration count
					apply_preconditioner,	// preconditioner callback
					matrix_multiply			// matrix multiplication callback
				);
				std::cout << "\tEnd CG callback\n";
				
				// update progress
				computations_done++;
				iterations_done += iterations;
				
				// save solution to disk
				current_solution.hdfsave(cur_oss.str().c_str(), true /* = with compression */);
				
			} // end of For mi, Spin
		} // end of For li
		
		// free Numeric objects
		for (auto luft = lufts.begin(); luft != lufts.end(); luft++)
			luft->free();
		
	} // end of For ie = 0, ..., Nenergy - 1
	
	std::cout << "\rSolving the systems... ok                                                            \n";
	std::cout << "\t(typically " << iterations_done/Nenergy << " CG iterations per energy)\n";
	
	// --------------------------------------------------------------------- //
	
	
	// exit if requested
	if (stg12)
	{
#ifndef NO_MPI
		if (parallel)
			MPI_Finalize();
#endif
		exit(0);
	}
	
	// Create SQL batch file
	char filename[100];
	if (parallel)
		std::sprintf(filename, "%d-%d-(%d).sql", ni, L, iproc);
	else
		std::sprintf(filename, "%d-%d.sql", ni, L);
	FILE* fsql = std::fopen(filename, "w");
	std::fprintf(fsql, "BEGIN TRANSACTION;\n");
	
	
	// Extract the cross sections ------------------------------------------ //
	//
	
	std::cout << "Extracting T-matrices...";
	
	std::vector<std::tuple<int,int,int,int,int>> transitions;
	for (int Spin = 0; Spin <= 1; Spin++)
	for (auto instate  : instates)
	for (auto outstate : outstates)
		transitions.push_back (
			std::make_tuple (
				Spin,
				/*li*/ std::get<1>(instate),
				/*mi*/ std::get<2>(instate),
				/*nf*/ std::get<0>(outstate),
				/*lf*/ std::get<1>(outstate)
			)
		);
	
	int finished = 0;
	
	for (int i = 0; i < (int)transitions.size(); i++)
	{
#ifndef NO_MPI
		// if MPI is active, compute only a subset of the transitions, corresponding to iproc
		if (parallel and i % Nproc != iproc)
			continue;
#endif
		
		int Spin = std::get<0>(transitions[i]);
		int li   = std::get<1>(transitions[i]);
		int mi   = std::get<2>(transitions[i]);
		int nf   = std::get<3>(transitions[i]);
		int lf   = std::get<4>(transitions[i]);
		
		// skip angular forbidden states
		bool allowed = false;
		for (int l = abs(li - L); l <= li + L; l++)
			allowed = allowed or ClebschGordan(li,mi,l,0,L,mi);
		if (not allowed)
			continue;
		
		cArray Pf_overlaps = overlapP(nf,lf,weight_end_damp);
		
		// compute radial integrals
		cArrays Lambda(2 * lf + 1);
		for (int mf = -lf; mf <= lf; mf++)
		{
			// compute expansions of final wave function for this final state (nf,lf)
			std::transform (
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
					std::fprintf(fsql, "INSERT OR REPLACE INTO \"tmat\" VALUES (%d,%d,%d, %d,%d,%d, %d,%d, %e, %d, %e,%e);\n",
						ni, li, mi, nf, lf, mf, L, Spin, Ei[ie], ell, T_ell[i].real(), T_ell[i].imag()
					);
				}
			}
			
			//
			// print out the total cross section for quick overview
			//
			
			std::sprintf(filename, "sigma-%d-%d-%d-%d-%d-%d-%d-%d.dat", ni, li, mi, nf, lf, mf, L, Spin);
			FILE* ftxt = std::fopen(filename, "w");
			std::fprintf(ftxt, "# Ei [Ry] sigma [a0^2]\n");
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
					std::fprintf(ftxt, "%g\t%g\n", Ei[ie], sigma);
				}
			}
			std::fclose(ftxt);
		}
		
		finished++;
		
		std::printf("\rExtracting T-matrices... %3.0f%%     ", finished * 100. / transitions.size());
	}
	std::cout << "\rExtracting T-matrices... ok       \n";
	// --------------------------------------------------------------------- //
	
	std::fprintf(fsql, "COMMIT;\n");
	std::fclose(fsql);
	
#ifndef NO_MPI
	if (parallel)
		MPI_Finalize();
#endif
	
	std::cout << "\nDone.\n\n";
	
	return 0;
}
