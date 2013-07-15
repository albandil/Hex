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

#include "version.h"

#ifndef SuiteSparse_time
double SuiteSparse_time()
{
	return 0.;
}
#endif

int main(int argc, char* argv[])
{
	std::cout <<
		"\n"
		"       / /   / /    __    \\ \\  / /\n"
		"      / /__ / /   / _ \\    \\ \\/ /\n"
		"     /  ___  /   | |/_/    / /\\ \\\n"
		"    / /   / /    \\_\\      / /  \\ \\\n"
		"\n"
		"             UK MFF (c) 2013\n"
		"\n"
		"          version hash: " << commit_hash << "\n"
		"\n";
	
	// Preparations ------------------------------------------------------- //
	//
	
	gsl_set_error_handler_off();
	H5::Exception::dontPrint();
    
	// disable buffering of the standard output 
	// (so that all text messages are immediatelly visible)
	setvbuf(stdout, nullptr, _IONBF, 0);
	
	// variables that can be set by the user from the command line
	std::ifstream inputfile;    // input file
	std::string zipfile;        // HDF solution expansion file to zip
	int  zipcount = 0;          // zip sample count
	double zipmax = -1;         // zip bounding box
	bool parallel = false;      // whether to use OpenMPI
	
	// which stages to run (default: all)
	int itinerary = StgNone;
	
	// get input from command line
	parse_command_line (
		argc, argv,
		inputfile,
		zipfile, zipcount, zipmax,
		parallel,
		itinerary
	);
	
	// run the whole sequence if nothing specified
	if (itinerary == StgNone)
		itinerary = StgRadial | StgSolve | StgExtract;
	
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
	int L, Spin, Pi;        // global conserved variables
	int levels;             // angular momentum levels
	rArray Ei;              // energies in Rydbergs
	double B = 0;           // magnetic field in a.u.
	
	// initial and final atomic states
	std::vector<std::tuple<int,int,int>> instates, outstates;
	
	// get input from input file
	parse_input_file (
		inputfile, order, ecstheta,
		real_knots, complex_knots, ni, 
		instates, outstates,
		L, Spin, levels, Ei, B
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
	
	
	// Setup angular data -------------------------------------------------- //
	//
	std::vector<std::pair<int,int>> coupled_states;
	std::vector<int> workers;
	
	std::cout << "Setting up the coupled angular states...\n";
	
	Pi = L % 2;
	
	// for given L, Π and levels list all available (ℓ₁ℓ₂) pairs
	for (int ell = 0; ell <= levels; ell++)
	{
		std::cout << "\t-> [" << ell << "] ";
		
		// get sum of the angular momenta for this angular level
		int sum = 2 * ell + L;
		
		// for all angular momentum pairs that do compose L
		for (int l1 = ell; l1 <= sum - ell; l1++)
		{
			std::cout << "(" << l1 << "," << sum - l1 << ") ";
			coupled_states.push_back(std::make_pair(l1, sum - l1));
		}
		std::cout << "\n";
	}
	
	std::cout << "\t-> The matrix of the set contains " << coupled_states.size() << " diagonal blocks.\n";
	
	// skip if nothing to compute
	if (coupled_states.empty())
		exit(0);
	
	// store maximal angular momentum
	int maxell = levels + L;
	
	std::cout << "\n";
	// --------------------------------------------------------------------- //
	
	
	if (zipfile.size() != 0 and I_am_master)
	{
		if (zipmax < 0)
			zipmax = Rmax;
		
		cArray sol;     // stored solution expansion
		cArray ev;      // evaluated solution
		rArray grid;    // real evaluation grid
		
		std::cout << "Zipping B-spline expansion of the solution: \"" << zipfile << "\"" << std::endl;
		
		if (not sol.hdfload(zipfile.c_str()))
			throw exception("Cannot load file %s.", zipfile.c_str());
		
		grid = linspace(0., zipmax, zipcount);
		
		for (unsigned ill = 0; ill < coupled_states.size(); ill++)
		{
			int l1 = coupled_states[ill].first;
			int l2 = coupled_states[ill].second;
			
			std::cout << "\t- partial wave l1 = " << l1 << ", l2 = " << l2 << "\n";
			
			// zip this partial wave
			ev = Bspline::ECS().zip (
				cArrayView (
					sol,
					ill * Nspline * Nspline,
					Nspline * Nspline
				),
				grid, grid
			);
			
			// write VTK header
			std::ostringstream ofname;
			ofname << zipfile << "_(" << l1 << "," << l2 << ").vtk";
			std::ofstream out(ofname.str().c_str());
			out << "# vtk DataFile Version 3.0\n";
			out << "Hex-ecs wave function partial waves\n";
			out << "ASCII\n";
			out << "DATASET RECTILINEAR_GRID\n";
			out << "DIMENSIONS " << zipcount << " " << zipcount << " 1\n";
			out << "X_COORDINATES " << zipcount << " float\n";
			out << grid.string() << "\n";
			out << "Y_COORDINATES " << zipcount << " float\n";
			out << grid.string() << "\n";
			out << "Z_COORDINATES 1 float\n";
			out << "0\n";
			out << "POINT_DATA " << zipcount * zipcount << "\n";
			out << "FIELD wavefunction 2\n";
			
			// save real part
			out << "pw_" << l1 << "_" << l2 << "_re 1 " << zipcount * zipcount << " float\n";
			for (int i = 0; i < zipcount; i++)
			{
				for (int j = 0; j < zipcount; j++)
					out << ev[i * zipcount + j].real() << " ";
				out << "\n";
			}
			
			// save imaginary part
			out << "pw_" << l1 << "_" << l2 << "_im 1 " << zipcount * zipcount << " float\n";
			for (int i = 0; i < zipcount; i++)
			{
				for (int j = 0; j < zipcount; j++)
					out << ev[i * zipcount + j].imag() << " ";
				out << "\n";
			}
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
	
	
	std::cout << "ok\n\nLoading/precomputing integral moments... ";
	
	
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
	
	
	std::cout << "ok\n\n";
	
	
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
	
	int maxlambda = L + 2 * levels;
	std::vector<SymDiaMatrix> R_tr_dia(maxlambda + 1);
	
Stg1:
{
	
	// skip two-electron integration if told so
	if (not (itinerary & StgRadial))
	{
		std::cout << "Skipped computation of two-electron integrals.\n";
		goto Stg2;
	}
	
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
	
	// for all multipoles : compute / load
	for (int lambda = 0; lambda <= maxlambda; lambda++)
	{
		// this process will only compute a subset of radial integrals
		if (lambda % Nproc != iproc)
			continue;
		
		// look for precomputed data on disk
		std::ostringstream oss2;
		oss2 << "R_tr_dia_" << lambda << ".hdf";
		int R_integ_exists = R_tr_dia[lambda].hdfload(oss2.str().c_str());
		
		if (R_integ_exists)
		{
			std::cout << "\t- integrals for λ = " << lambda << " loaded from \"" << oss2.str().c_str() << "\"\n";
			continue; // no need to compute
		}
		
		// precompute necessary partial integral moments
		// - truncated moments λ and -λ-1
		cArray Mtr_L    = computeMi( lambda,   Nreknot-1);
		cArray Mtr_mLm1 = computeMi(-lambda-1, Nreknot-1);
		
		// elements of R_tr
		NumberArray<long> R_tr_i, R_tr_j, th_R_tr_i, th_R_tr_j;
		NumberArray<Complex> R_tr_v, th_R_tr_v;
		
		# pragma omp parallel default(none) \
			private (th_R_tr_i, th_R_tr_j, th_R_tr_v) \
			firstprivate (Nspline, order, lambda, Mtr_L, Mtr_mLm1) \
			shared (R_tr_i, R_tr_j, R_tr_v)
		{
			// for all B-spline pairs
			# pragma omp for schedule(dynamic,1)
			for (int i = 0; i < Nspline; i++)
			for (int j = 0; j < Nspline; j++)
			{
				// for all nonzero, nonsymmetry R-integrals
				for (int k = i; k <= i + order and k < Nspline; k++) // enforce i ≤ k
				for (int l = j; l <= j + order and l < Nspline; l++) // enforce j ≤ l
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
				R_tr_i.append(th_R_tr_i.begin(), th_R_tr_i.end());
				R_tr_j.append(th_R_tr_j.begin(), th_R_tr_j.end());
				R_tr_v.append(th_R_tr_v.begin(), th_R_tr_v.end());
			}
		}
		
		// create matrices and save them to disk
		R_tr_dia[lambda] = CooMatrix(Nspline*Nspline, Nspline*Nspline, R_tr_i, R_tr_j, R_tr_v).todia();
		R_tr_dia[lambda].hdfsave(oss2.str().c_str(), true, 10);
		
		std::cout << "\t- integrals for λ = " << lambda << " computed\n";
	}
	
	// for all multipoles : synchronize
	#ifndef NO_MPI
	if (parallel)
	{
		for (int lambda = 0; lambda <= maxlambda; lambda++)
		{
			// get owner process of this multipole
			int owner = lambda % Nproc;
			
			// get dimensions
			int diagsize = R_tr_dia[lambda].diag().size();
			int datasize = R_tr_dia[lambda].data().size();
			
			// owner will broadcast dimensions
			MPI_Bcast(&diagsize, 1, MPI_INT, owner, MPI_COMM_WORLD);
			MPI_Bcast(&datasize, 1, MPI_INT, owner, MPI_COMM_WORLD);
			
			// get arrays
			NumberArray<int>     diag = R_tr_dia[lambda].diag();
			NumberArray<Complex> data = R_tr_dia[lambda].data();
			diag.resize(diagsize);
			data.resize(datasize);
			
			// master will broadcast arrays
			MPI_Bcast(&diag[0], diag.size(), MPI_INT, owner, MPI_COMM_WORLD);
			MPI_Bcast(&data[0], data.size(), MPI_DOUBLE_COMPLEX, owner, MPI_COMM_WORLD);
			
			// reconstruct objects
			R_tr_dia[lambda] = SymDiaMatrix(Nspline * Nspline, diag, data);
			
			if (owner != iproc)
			{
				std::cout << "\t- integrals for λ = " << lambda << " retrieved from process " << owner << "\n";
				
				// save to disk
				std::ostringstream oss3;
				oss3 << "R_tr_dia_" << lambda << ".hdf";
				R_tr_dia[lambda].hdfsave(oss3.str().c_str(), true, 10);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	#endif
	
	std::cout << "\t- R_tr[λ] has " << R_tr_dia[0].data().size() << " nonzero elements\n";
	// --------------------------------------------------------------------- //
	
	
	std::cout << "Hamiltonian properties:\n";
	std::cout << "\t-> hamiltonian size: " << Nspline * Nspline * coupled_states.size() << "\n\n";
	
}
Stg2:
{
	// skip stage 2 if told so
	if (not (itinerary & StgSolve))
	{
		std::cout << "Skipped solution of the equation.\n";
		goto Stg3;
	}
	
	
	// Prepare B-spline overlaps and expansions of Ric-Bess functions ------ //
	//
	std::cout << "Computing B-spline expansions... ";
	
	//  j-overlaps of shape [Nenergy × Nangmom × Nspline]
	cArray ji_overlaps = overlapj(maxell,ki,weightEdgeDamp());
	ji_overlaps.hdfsave("ji_overlaps_damp.hdf"); // just for debugging
	
	//  compute expansions; solve the system
	//      S * B_spline_expansion = B_spline_overlap
	unsigned ji_expansion_count = ji_overlaps.size()/Nspline;
	cArray ji_expansion = S.tocoo().tocsr().solve(ji_overlaps, ji_expansion_count);
	ji_expansion.hdfsave("ji_expansion.hdf"); // just for debugging
	
	std::cout << "ok\n\n";
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
	std::cout << "ok\n\n";
	// --------------------------------------------------------------------- //
	
	
	// Distribute LU factorizations among processes ------------------------ //
	//
	std::map<int,int> LUs;
	std::vector<int> info(Nproc);
	int worker = 0;
	std::cout << "Balancing " << coupled_states.size()
	          << " diagonal blocks among " << Nproc 
	          << " worker processes...\n";
	for (unsigned ill = 0; ill < coupled_states.size(); ill++)
	{
		info[worker]++;                       // add work to the process 'worker'
		
		LUs[ill] = worker;	              // assign this block to worker
		worker = (worker + 1) % Nproc;        // move to next worker
	}
	// print statistics
	std::cout << "\t-> average " << coupled_states.size()/double(Nproc) << " blocks/process\n";
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
		std::vector<CsrMatrix> csr_blocks(coupled_states.size());
		
		// diagonal blocks in DIA format
		// - these will be used in matrix multiplication
		std::vector<SymDiaMatrix> dia_blocks(coupled_states.size());
		
		// LU factorization of the diagonal blocks
		std::vector<CsrMatrix::LUft> lufts(coupled_states.size());
		
		// setup preconditioner - the diagonal block LU-factorizations
		std::cout << "\tSetup preconditioner blocks... " << std::flush;
		for (unsigned ill = 0; ill < coupled_states.size(); ill++)
		{
			// skip computation of unwanted blocks for this process
			if (LUs.find(ill) == LUs.end() or LUs[ill] != iproc)
				continue;
			
			int l1 = coupled_states[ill].first;
			int l2 = coupled_states[ill].second;
			
			// one-electron parts
			SymDiaMatrix Hdiag =
			    half_D_minus_Mm1_tr_kron_S
			    + ((0.5*l1*(l1+1)) * Mm2_kron_S)
			    + S_kron_half_D_minus_Mm1_tr
			    + (0.5*l2*(l2+1)) * S_kron_Mm2;
			
			// two-electron part
			for (int lambda = 0; lambda <= maxlambda; lambda++)
			{
				Complex f = computef(lambda,l1,l2,l1,l2,L);
				if (f != 0.)
					Hdiag += f * R_tr_dia[lambda];
			}
			
			// finalize the matrix
			dia_blocks[ill] = E*S_kron_S - Hdiag;
			csr_blocks[ill] = dia_blocks[ill].tocoo().tocsr();
		}
		std::cout << "ok\n";
		
		// compute the LU factorizations
		for (unsigned ill = 0; ill < coupled_states.size(); ill++)
		{
			int l1 = coupled_states[ill].first;
			int l2 = coupled_states[ill].second;
			
			// skip computation of unwanted blocks for this process
			if (LUs.find(ill) == LUs.end() or LUs[ill] != iproc)
				continue;
			
			// timer info
			std::chrono::steady_clock::time_point start;
			std::chrono::duration<int> sec;
			
			// log output
			std::cout << "\t[" << iproc << "] LU factorization " 
			          << ill << " of (" << l1 << "," << l2 << ") block started\n";
			start = std::chrono::steady_clock::now();
			
			// factorize
			lufts[ill] = csr_blocks[ill].factorize();
			
			// log output
			sec = std::chrono::duration_cast<std::chrono::duration<int>>(std::chrono::steady_clock::now()-start);
			std::cout << "\t[" << iproc << "] LU factorization " 
			          << ill << " of (" << l1 << "," << l2 << ") block done after " 
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
			Pi_overlaps = overlapP(ni,li,weightEndDamp());
			Pi_expansion = S.tocoo().tocsr().solve(Pi_overlaps);
			
			// we may have already computed solution for this state and energy... is it so?
			std::ostringstream cur_oss;
			cur_oss << "psi-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << "-" << Ei[ie] << ".hdf";
			if ( current_solution.hdfload(cur_oss.str().c_str()) )
				continue;
			
			// create right hand side
			cArray chi ( coupled_states.size()*Nspline*Nspline );
			
			std::cout << "\tCreate RHS for li = " << li << ", mi = " << mi << "\n";
			
			// for all segments constituting the RHS
			# pragma omp parallel for 
			for (unsigned ill = 0; ill < coupled_states.size(); ill++)
			{
				int l1 = coupled_states[ill].first;
				int l2 = coupled_states[ill].second;
				
				// setup storage
				cArrayView chi_block(chi, ill * Nspline * Nspline, Nspline * Nspline);
				chi_block.clear();
				
				// for all allowed angular momenta (by momentum composition) of the projectile
				for (int l = abs(li - L); l <= li + L; l++)
				{
					// skip wrong parity
					if ((li + l) % 2 != Pi % 2)
						continue;
					
					// (anti)symmetrization
					int Sign = (Spin % 2 == 0) ? 1. : -1.;
					
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
			
			std::cout << "\tSolve the equations.\n";
			
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
				# pragma omp parallel for schedule (dynamic,1)
				for (unsigned ill = 0; ill < coupled_states.size(); ill++)
				{
					// skip computation of unwanted blocks for this process
					if (LUs.find(ill) == LUs.end() or LUs[ill] != iproc)
						continue;
					
					// create copy-to view of "z"
					cArrayView zview(z, ill * Nspline * Nspline, Nspline * Nspline);
					
					// create copy-from view of "r"
					cArrayView rview(r, ill * Nspline * Nspline, Nspline * Nspline);
					
					// copy the correcponding slice of "r" multiplied by a correct block inversion
					zview = lufts[ill].solve(rview);
				}
				
#ifndef NO_MPI
				if (parallel)
				{
					// synchronize across processes
					for (unsigned ill = 0; ill < coupled_states.size(); ill++)
					{
						// relevant process will broadcast this segment's data
						MPI_Bcast (
							&z[0] + ill * Nspline * Nspline,
							Nspline * Nspline,
							MPI_DOUBLE_COMPLEX,
							LUs[ill],
							MPI_COMM_WORLD
						);
					}
				}
#endif
			};
			
			// CG matrix multiplication callback
			auto matrix_multiply = [ & ](cArray const & p, cArray & q) -> void
			{
				// clear all output segments that are going to be referenced by this process
				for (unsigned ill = 0; ill < coupled_states.size(); ill++)
					if (LUs.find(ill) != LUs.end() and LUs[ill] == iproc)
						cArrayView(q, ill * Nspline * Nspline, Nspline * Nspline).clear();
				
				// multiply "q" by the matrix of the system
				# pragma omp parallel for schedule (dynamic,1) collapse(2)
				for (unsigned ill = 0; ill < coupled_states.size(); ill++)
				for (unsigned illp = 0; illp < coupled_states.size(); illp++)
				if (LUs.find(ill) != LUs.end() and LUs[ill] == iproc)
				{
					// row multi-index
					int l1 = coupled_states[ill].first;
					int l2 = coupled_states[ill].second;
					
					// column multi-index
					int l1p = coupled_states[illp].first;
					int l2p = coupled_states[illp].second;
					
					// product segment contribution
					cArray q_contrib (Nspline * Nspline);
					
					// copy-from segment of "p"
					cArrayView p_block(p, illp * Nspline * Nspline, Nspline * Nspline);
					
					// multiply by hamiltonian terms
					if (ill == illp)
					{
						// reuse the diagonal block
						q_contrib += dia_blocks[ill].dot(p_block);
					}
					else
					{
						// compute the offdiagonal block
						for (int lambda = 0; lambda <= maxlambda; lambda++)
						{
							Complex f = computef(lambda, l1, l2, l1p, l2p, L);
							if (f != 0.)
								q_contrib -= f * R_tr_dia[lambda].dot(p_block);
						}
					}
					
					// safely update shared output array "q"
					# pragma omp critical
					cArrayView (q, ill * Nspline * Nspline, Nspline * Nspline) += q_contrib;
				}
				
#ifndef NO_MPI
				if (parallel)
				{
					// synchronize across processes
					for (unsigned ill = 0; ill < coupled_states.size(); ill++)
					{
						// relevant process will broadcast this segment's data
						MPI_Bcast (
							&q[0] + ill * Nspline * Nspline,
							Nspline * Nspline,
							MPI_DOUBLE_COMPLEX,
							LUs[ill],
							MPI_COMM_WORLD
						);
					}
				}
#endif
			};
				
			if (chi.norm() == 0.)
			{
				std::cout << "\t! Right-hand-side is zero (probably due to incompatible angular settings).\n";
				current_solution.clear();
			}
			else
			{
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
				iterations_done += iterations;
			}
			
			computations_done++;
			
			// save solution to disk
			current_solution.hdfsave(cur_oss.str().c_str(), true /* = with compression */);
			
		} // end of For li
		
		// free Numeric objects
		for (auto luft = lufts.begin(); luft != lufts.end(); luft++)
			luft->free();
		
	} // end of For ie = 0, ..., Nenergy - 1
	
	std::cout << "\rSolving the systems... ok                                                            \n";
	std::cout << "\t(typically " << iterations_done/Nenergy << " CG iterations per energy)\n";
	
	// --------------------------------------------------------------------- //
	
}
Stg3:
{
	// skip stage 3 if told so
	if (not (itinerary & StgExtract))
	{
		std::cout << "Skipped extraction of amplitudes.\n";
		goto End;
	}
	
	// compose output filename
	std::ostringstream ossfile;
	if (parallel)
		ossfile << ni << "-" << L << "-" << Spin << "-(" << iproc << ").sql";
	else
		ossfile << ni << "-" << L << "-" << Spin << ".sql";
	
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
			cArray Pf_overlaps = overlapP(nf,lf,weightEndDamp());
			
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
				Lambda[mf+lf] = std::move(computeLambda(kf, ki, maxell, L, Spin, ni, li, mi, Ei, lf, Pf_overlaps, coupled_states));
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
			
			rArray ics;
			cArrays data = std::move(computeXi(maxell, L, Spin, ni, li, mi, Ei, ics, coupled_states));
			
			for (size_t ie = 0; ie < Ei.size(); ie++)
			for (unsigned ill = 0; ill < coupled_states.size(); ill++) //??? or triangular
			{
				// save data as BLOBs
				fsql << "INSERT OR REPLACE INTO \"ionf\" VALUES ("
					 << ni << "," << li << "," << mi << ","
					 << L  << "," << Spin << ","
					 << Ei[ie] << "," << coupled_states[ill].first << ","
					 << coupled_states[ill].second << ","
					 << data[ie * coupled_states.size() + ill].toBlob() << ");\n";
			}
			
			// print ionization cross section
			std::ostringstream fname;
			fname << "isigma-" << L << "-" << Spin << "-" << ni << "-" << li << "-" << mi << ".dat";
			write_array(Ei, ics, fname.str().c_str());
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
}
End:
{	
#ifndef NO_MPI
	if (parallel)
		MPI_Finalize();
#endif
	
	std::cout << "\nDone.\n\n";
}	
	return 0;
}
