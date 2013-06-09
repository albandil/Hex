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
#include <cstdio>
#include <cmath>
#include <chrono>
#include <cstring>
#include <vector>
#include <iostream>
#include <string>
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
	setvbuf(stdout, nullptr, _IONBF, 0);
	
	// variables that can be set by the user from the command line
	std::ifstream inputfile;    // input file
	std::string zipfile;        // HDF solution expansion file to zip
	int  zipcount = 0;          // zip sample count
	bool parallel = false;      // whether to use OpenMPI
	bool stg1 = false;          // whether to computate radial integrals only
	bool stg12 = false;         // whether to computate radial integrals and solutions only
	
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
		parallel = false;
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
	int order;              // B-spline order
	double ecstheta;        // ECS rotation
	rArray real_knots;      // real knot sequence
	rArray complex_knots;   // complex knot sequence
	int ni;                 // atomic state parameters
	int L;                  // global conserved variables
	int maxell;             // angular momentum limits
	rArray Ei;              // energies in Rydbergs
	double B = 0;           // magnetic field in a.u.
	
	// initial and final atomic states
	std::vector<std::tuple<int,int,int>> instates, outstates;
	
	// get input from input file
	parse_input_file (
		inputfile, order, ecstheta,
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
	std::transform (
		Ei.begin(), Ei.end(),
		ki.begin(),
		[ = ](double E) -> double { return sqrt(E); }
	);
	// --------------------------------------------------------------------- //
	
	
	// Setup B-spline environment ------------------------------------------ //
	//
	double R0 = real_knots.back();          // end of real grid
	double Rmax = complex_knots.back();     // end of complex grid
	
	if (R0 >= Rmax)
		throw exception("ERROR: Rmax = %g (end of grid) must be greater than R0 = %g (end of real grid)!", Rmax, R0);

	Bspline::ECS().init(order, real_knots, ecstheta, complex_knots);
	
	// shortcuts
	int Nspline = Bspline::ECS().Nspline();
	int Nknot   = Bspline::ECS().Nknot();
	int Nreknot = Bspline::ECS().Nreknot();
	
	// info
	std::cout << "B-spline total count: " << Nspline << "\n\n";
	//
	// --------------------------------------------------------------------- //
	
	
	if (zipfile.size() != 0 and I_am_master)
	{
		cArray sol;     // stored solution expansion
		cArray ev;      // evaluated solution
		rArray grid;    // real evaluation grid
		
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
			std::ostringstream outf1, outf2;
			outf1 << zipfile << "-(" << l1 << "," << l2 << ").re";
			outf2 << zipfile << "-(" << l1 << "," << l2 << ").im";
			
			// write real part
			write_2D_data (
				zipcount + 1,
				zipcount + 1,
				outf1.str().c_str(),
				[&](size_t i, size_t j) -> double {
					return ev[i * (zipcount + 1) + j].real();
				}
			);
			
			// write imaginary part
			write_2D_data (
				zipcount + 1,
				zipcount + 1,
				outf2.str().c_str(),
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
	SymDiaMatrix D(Nspline);
	//
	D.hdfload("D.hdf") or D.populate (
		order, [=](int i, int j) -> Complex { return computeD(i, j, Nknot - 1); }
	).hdfsave("D.hdf");
	// --------------------------------------------------------------------- //
	
	
	std::cout << "ok\nLoading/precomputing integral moments... ";
	
	
	// Precompute useful integral moments ---------------------------------- //
	//
	SymDiaMatrix S(Nspline);
	SymDiaMatrix Mm1(Nspline), Mm1_tr(Nspline);
	SymDiaMatrix Mm2(Nspline);
	//
	S.hdfload("S.hdf") or S.populate (
		order, [=](int m, int n) -> Complex { return computeM(0, m, n); }
	).hdfsave("S.hdf");
	//
	Mm1.hdfload("Mm1.hdf") or Mm1.populate (
		order, [=](int m, int n) -> Complex { return computeM(-1, m, n); }
	).hdfsave("Mm1.hdf");
	//
	Mm1_tr.hdfload("Mm1_tr.hdf") or Mm1_tr.populate (
		order,	[=](int m, int n) -> Complex { return computeM(-1, m, n, Nreknot - 1);}
	).hdfsave("Mm1_tr.hdf");
	//
	Mm2.hdfload("Mm2.hdf") or Mm2.populate (
		order, [=](int m, int n) -> Complex { return computeM(-2, m, n); }
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
	std::vector<SymDiaMatrix> R_tr_dia(maxlambda + 1);
	
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
	
	// for all multipoles
	for (int lambda = 0; lambda <= maxlambda; lambda++)
	{
		// look for precomputed data on disk
		std::ostringstream oss2;
		oss2 << "R_tr_dia_" << lambda << ".hdf";
		
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
			R_integ_exists = R_tr_dia[lambda].hdfload(oss2.str().c_str());
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
				int diagsize = R_tr_dia[lambda].diag().size();
				int datasize = R_tr_dia[lambda].data().size();
				
				// master will broadcast dimensions
				MPI_Bcast(&diagsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&datasize, 1, MPI_INT, 0, MPI_COMM_WORLD);
				
				// get arrays
				Array<int>     diag = R_tr_dia[lambda].diag();
				Array<Complex> data = R_tr_dia[lambda].data();
				diag.resize(diagsize);
				data.resize(datasize);
				
				// master will broadcast arrays
				MPI_Bcast(&diag[0], diag.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
				MPI_Bcast(&data[0], data.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
				
				// reconstruct objects
				R_tr_dia[lambda] = SymDiaMatrix(Nspline * Nspline, diag, data);
			}
		}
#endif
		
		if (R_integ_exists)
		{
			std::cout << "\t- integrals for λ = " << lambda << " successfully loaded\n";
			continue; // no need to compute
		}
		
		// precompute necessary partial integral moments
		// - truncated moments λ and -λ-1
		cArray Mtr_L    = computeMi( lambda,   Nreknot-1);
		cArray Mtr_mLm1 = computeMi(-lambda-1, Nreknot-1);
		
		// elements of R_tr
		std::vector<long> R_tr_i, R_tr_j, th_R_tr_i, th_R_tr_j;
		std::vector<Complex> R_tr_v, th_R_tr_v;
		
		# pragma omp parallel default(none) \
			private (th_R_tr_i, th_R_tr_j, th_R_tr_v) \
			firstprivate (Nspline, order, lambda, Mtr_L, Mtr_mLm1) \
			shared (R_tr_i, R_tr_j, R_tr_v)
		{
			// for all B-spline pairs
			# pragma omp for collapse(2)
			for (int i = 0; i < Nspline; i++)
			for (int j = 0; j < Nspline; j++)
			{
				// for all nonzero, nonsymmetry R-integrals
				for (int k = i; k < i + order and k < Nspline; k++) // enforce i ≤ k
				for (int l = j; l < j + order and l < Nspline; l++) // enforce j ≤ l
				{
					// skip symmetry ijkl <-> jilk (others are accounted for in the limits
					if (i > j and k > l)
						continue;
					
					// evaluate B-spline integral
					Complex Rijkl_tr = computeR(lambda, i, j, k, l, Mtr_L, Mtr_mLm1);
					
					// store all symmetries
					allSymmetries(i, j, k, l, Rijkl_tr, th_R_tr_i, th_R_tr_j, th_R_tr_v);
				}
			}
			
			# pragma omp critical
			{
				// merge the thread local arrays
				R_tr_i.insert(R_tr_i.end(), th_R_tr_i.begin(), th_R_tr_i.end());
				R_tr_j.insert(R_tr_j.end(), th_R_tr_j.begin(), th_R_tr_j.end());
				R_tr_v.insert(R_tr_v.end(), th_R_tr_v.begin(), th_R_tr_v.end());
			}
		}
		
		// create matrices and save them to disk
		R_tr_dia[lambda] = CooMatrix(Nspline*Nspline, Nspline*Nspline, R_tr_i, R_tr_j, R_tr_v).todia();
		R_tr_dia[lambda].hdfsave(oss2.str().c_str());
		
		std::cout << "\r\t- multipole λ = " << lambda << "... ok            \n";
		
	}
	std::cout << "\t- R_tr[λ] has " << R_tr_dia[0].data().size() << " nonzero elements\n";
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
	cArray ji_expansion = S.tocoo().tocsr().solve(ji_overlaps, ji_expansion_count);
	ji_expansion.hdfsave("ji_expansion.hdf"); // just for debugging
	
	std::cout << "ok\n";
	// --------------------------------------------------------------------- //
	
	
	// Precompute some accelerators ---------------------------------------- //
	//
	// Kronecker producs
	std::cout << "Creating Kronecker products... ";
	SymDiaMatrix S_kron_S   = S.kron(S);
	SymDiaMatrix S_kron_Mm1_tr = S.kron(Mm1_tr);
	SymDiaMatrix S_kron_Mm2 = S.kron(Mm2);
	SymDiaMatrix Mm1_tr_kron_S = Mm1_tr.kron(S);
	SymDiaMatrix Mm2_kron_S = Mm2.kron(S);
	SymDiaMatrix half_D_minus_Mm1_tr = 0.5 * D - Mm1_tr;
	SymDiaMatrix half_D_minus_Mm1_tr_kron_S = half_D_minus_Mm1_tr.kron(S);
	SymDiaMatrix S_kron_half_D_minus_Mm1_tr = S.kron(half_D_minus_Mm1_tr);
	std::cout << "ok\n";
	// --------------------------------------------------------------------- //
	
	
	// Distribute LU factorizations among processes ------------------------ //
	//
	std::map<int,int> LUs;
	std::vector<int> info(Nproc);
	int worker = 0;
	std::cout << "\nBalancing " << triangle_count(L,maxell)
	          << " diagonal blocks among " << Nproc 
	          << " worker processes...\n";
	for (int l1 = 0; l1 <= maxell; l1++)
	for (int l2 = 0; l2 <= maxell; l2++)
	{
		// skip those angular momentum pairs that don't compose L
		if (abs(l1 - l2) > L or l1 + l2 < L)
			continue;
		
		info[worker]++;                       // add work to the process 'worker'
		
		int iblock = l1 * (maxell + 1) + l2;  // get block index
		LUs[iblock] = worker;	              // assign this block to worker
		worker = (worker + 1) % Nproc;        // move to next worker
	}
	// print statistics
	std::cout << "\t-> average " << triangle_count(L,maxell)/double(Nproc) << " blocks/process\n";
	std::cout << "\t-> min " << *std::min_element(info.begin(), info.end()) << " blocks/process\n";
	std::cout << "\t-> max " << *std::max_element(info.begin(), info.end()) << " blocks/process\n";
	// --------------------------------------------------------------------- //
	
	
	// For all right hand sides -------------------------------------------- //
	//
	int iterations_done = 0, computations_done = 0;
	for (unsigned ie = 0; ie < Nenergy; ie++)
	{
		// print progress information
		std::cout << "\nSolving the system for Ei[" << ie << "] = " << Ei[ie] << " ("
		          << int(trunc(ie * 100. / Nenergy + 0.5)) << " % finished, typically "
		          << (computations_done == 0 ? 0 : iterations_done / computations_done)
		          << " CG iterations per energy)\n";
		
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
		
		// diagonal blocks in CSR format
		// - these will be used for UMFPACK factorization
		// - need to exist physically for the LUft object to be valid!
		std::vector<CsrMatrix> csr_blocks((maxell+1)*(maxell+1));
		
		// diagonal blocks in DIA format
		// - these will be used in matrix multiplication
		std::vector<SymDiaMatrix> dia_blocks((maxell+1)*(maxell+1));
		
		// LU factorization of the diagonal blocks
		std::vector<CsrMatrix::LUft> lufts((maxell+1)*(maxell+1));
		
		// setup preconditioner - the diagonal block LU-factorizations
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
			SymDiaMatrix Hdiag =
			    half_D_minus_Mm1_tr_kron_S
			    + ((0.5*l1*(l1+1)) * Mm2_kron_S)
			    + S_kron_half_D_minus_Mm1_tr;
			    + (0.5*l2*(l2+1)) * S_kron_Mm2;
			
			// two-electron part
			for (int lambda = 0; lambda <= maxlambda; lambda++)
			{
				Complex f = computef(lambda,l1,l2,l1,l2,L);
				if (f != 0.)
					Hdiag += f * R_tr_dia[lambda];
			}
			
			// finalize the matrix
			dia_blocks[iblock] = E*S_kron_S - Hdiag;
			csr_blocks[iblock] = dia_blocks[iblock].tocoo().tocsr();
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
			lufts[iblock] = csr_blocks[iblock].factorize();
			
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
			Pi_expansion = S.tocoo().tocsr().solve(Pi_overlaps);
			
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
						
						// pick the correct Bessel function expansion
						cArrayView Ji_expansion (
							ji_expansion,
							Nspline * (ie * (maxell + 1) + l),
							Nspline
						);
						
						// compute outer products of B-spline expansions
						cArray Pj1 = outer_product(Pi_expansion, Ji_expansion);
						cArray Pj2 = outer_product(Ji_expansion, Pi_expansion);
						
						// skip angular forbidden right hand sides
						for (int lambda = 0; lambda <= maxlambda; lambda++)
						{
							Complex f1 = computef(lambda, l1, l2, li, l, L);
							Complex f2 = computef(lambda, l1, l2, l, li, L);
							
							if (f1 != 0.)
							{
								chi_block += (prefactor * f1) * R_tr_dia[lambda].dot(Pj1);
							}
							
							if (f2 != 0.)
							{
								if (Sign > 0)
									chi_block += (prefactor * f2) * R_tr_dia[lambda].dot(Pj2);
								else
									chi_block -= (prefactor * f2) * R_tr_dia[lambda].dot(Pj2);
							}
						}
						
						if (li == l1 and l == l2)
						{
							// direct contribution
							chi_block -= prefactor * S_kron_Mm1_tr.dot(Pj1);
						}
						
						if (li == l2 and l == l1)
						{
							// exchange contribution with the correct sign
							if (Sign > 0)
								chi_block -= prefactor * Mm1_tr_kron_S.dot(Pj2);
							else
								chi_block += prefactor * Mm1_tr_kron_S.dot(Pj2);
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
								q_block += dia_blocks[iblock].dot(p_block);
							}
							else
							{
								// compute the offdiagonal block
								for (int lambda = 0; lambda <= maxlambda; lambda++)
								{
									Complex _f = computef(lambda, l1, l2, l1p, l2p, L);
									if (_f != 0.)
										q_block -= _f * R_tr_dia[lambda].dot(p_block);
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
	
	// compose output filename
	std::ostringstream ossfile;
	if (parallel)
		ossfile << ni << "-" << L << "-(" << iproc << ").sql";
	else
		ossfile << ni << "-" << L << ".sql";
	
	// Create SQL batch file
	std::ofstream fsql(ossfile.str().c_str());
	
	// set exponential format for floating point output
	fsql.setf(std::ios_base::scientific);
	
	// write header
	fsql << "BEGIN TRANSACTION;\n";
	
	
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
		
		if (nf > 0)
		{
			//
			// Discrete transition
			//
			
			// precompute hydrogen function overlaps
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
				
				// compute Λ for transitions to (nf,lf,mf); it will depend on [ie,ℓ]
				Lambda[mf+lf] = std::move(computeLambda(kf, ki, maxell, L, Spin, ni, li, mi, Ei, lf, Pf_overlaps));
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
						fsql << "INSERT OR REPLACE INTO \"tmat\" VALUES ("
						     << ni << "," << li << "," << mi << ","
						     << nf << "," << lf << "," << mf << ","
							 << L  << "," << Spin << ","
							 << Ei[ie] << "," << ell << "," 
							 << T_ell[i].real() << "," << T_ell[i].imag()
							 << ");\n";
					}
				}
				
				//
				// print out the total cross section for quick overview
				//
				
				std::ostringstream sigmaname;
				sigmaname << "sigma-" 
				          << ni << "-" << li << "-" << mi << "-"
						  << nf << "-" << lf << "-" << mf << "-"
						  << L << "-" << Spin << ".dat";

				std::ofstream ftxt(sigmaname.str().c_str());
				
				ftxt << "# Ei [Ry] sigma [a0^2]\n";
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
						ftxt << Ei[ie] << "\t" << sigma << "\n";
					}
				}
				ftxt.close();
			}
		}
		else
		{
			//
			// Ionization
			//
			
			cArrays data = std::move(computeXi(maxell, L, Spin, ni, li, mi, Ei));
			cArrays::const_iterator iter = data.begin();
			
			for (size_t ie = 0; ie < Ei.size(); ie++)
			for (int l1 = 0; l1 <= maxell; l1++)
			for (int l2 = std::abs(l1-L); l2 <= l1; l2++)
			{
				// save data as BLOBs
				fsql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
					 << ni << "," << li << "," << mi << ","
					 << L  << "," << Spin << ","
					 << Ei[ie] << "," << l1 << "," << l2 << ","
					 << iter->toBlob() << ");\n";
				
				// move to next data
				iter++;
			}
		}
		
		finished++;
		
		std::cout << "\rExtracting T-matrices... " 
		          << std::setw(3) << int(trunc(finished * 100. / transitions.size() + 0.5))
				  << " %        ";
	}
	std::cout << "\rExtracting T-matrices... ok       \n";
	// --------------------------------------------------------------------- //
	
	fsql << "COMMIT;\n";
	fsql.close();
	
#ifndef NO_MPI
	if (parallel)
		MPI_Finalize();
#endif
	
	std::cout << "\nDone.\n\n";
	
	return 0;
}
