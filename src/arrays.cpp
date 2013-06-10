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
#include <vector>

#ifndef NO_HDF
	#include <H5Cpp.h>
#endif

#include "arrays.h"
#include "complex.h"

bool all(std::vector<bool> v)
{
	for (auto it = v.begin(); it != v.end(); it++)
		if (not *it)
			return false;
	return true;
}

bool any(std::vector<bool> v)
{
	for (auto it = v.begin(); it != v.end(); it++)
		if (*it)
			return true;
	return false;
}

#ifndef NO_HDF
bool save_array(rArray const & vec, const char* name, const double * const pdelta)
{
	try
	{
		H5::H5File h5file(name, H5F_ACC_TRUNC);
		
		int rank = 1;
		hsize_t length = vec.size();
		
		// save data
		H5::DataSpace dspc(rank, &length);
		H5::IntType dtype(H5::PredType::NATIVE_DOUBLE);
		H5::DataSet dset = h5file.createDataSet("data", dtype, dspc);
		dset.write(&vec[0], H5::PredType::NATIVE_DOUBLE);
		
		// save delta (if wanted)
		if (pdelta != 0)
		{
			length = 1;
			H5::DataSpace dspc(rank, &length);
			H5::IntType dtype(H5::PredType::NATIVE_DOUBLE);
			H5::DataSet dset = h5file.createDataSet("delta", dtype, dspc);
			dset.write(pdelta, H5::PredType::NATIVE_DOUBLE);
		}
		
		return true;
	}
	catch (...)
	{
		return false;
	}
}

bool save_array(cArray const & vec, const char* name)
{
	try
	{
		H5::H5File h5file(name, H5F_ACC_TRUNC);
		
		int rank = 1;
		hsize_t length = 2 * vec.size();
		
		// save data as an interleaved array
		H5::DataSpace dspc(rank, &length);
		H5::IntType dtype(H5::PredType::NATIVE_DOUBLE);
		H5::DataSet dset = h5file.createDataSet("data", dtype, dspc);
		dset.write(reinterpret_cast<const double*>(&vec[0]), H5::PredType::NATIVE_DOUBLE);
		
		return true;
	}
	catch (...)
	{
		return false;
	}
}

bool load_array(rArray & vec, const char* name, double* pdelta)
{
	try
	{
		H5::Exception::dontPrint();
		H5::H5File h5file(name, H5F_ACC_RDONLY);
		
		// load data 
		H5::DataSet dset = h5file.openDataSet("data");
		H5::DataSpace dspc = dset.getSpace();
		size_t N = dspc.getSimpleExtentNpoints();
		vec.resize(N);
		dset.read(&vec[0], H5::PredType::NATIVE_DOUBLE, dspc, dspc);
		
		// load delta (if wanted)
		if (pdelta != 0)
		{
			dset = h5file.openDataSet("delta");
			dspc = dset.getSpace();
			dset.read(pdelta, H5::PredType::NATIVE_DOUBLE, dspc, dspc);
		}
		
		return true;
	}
	catch (...)
	{
		return false;
	}
}

bool load_array(cArray & vec, const char* name)
{
	try
	{
		H5::Exception::dontPrint();
		H5::H5File h5file(name, H5F_ACC_RDONLY);
		
		// load data 
		H5::DataSet dset = h5file.openDataSet("data");
		H5::DataSpace dspc = dset.getSpace();
		size_t N = dspc.getSimpleExtentNpoints();
		vec.resize(N/2);
		dset.read(reinterpret_cast<Complex*>(&vec[0]), H5::PredType::NATIVE_DOUBLE, dspc, dspc);
		
		return true;
	}
	catch (...)
	{
		return false;
	}
}
#endif

void write_array(const std::map<unsigned long long, Complex>& m, const char* filename)
{
	FILE* f = fopen(filename, "w");
	for (auto it = m.begin(); it != m.end(); it++)
	{
		unsigned short lam = it->first % 65536L;
		unsigned short l2 = (it->first/65536L) % 65536L;
		unsigned short l1 = (it->first/(65536L*65536L)) % 65536L;
		fprintf(f, "%d %d %d %g %g\n", l1, l2, lam, it->second.real(), it->second.imag());
	}
	fclose(f);
}

void write_array(rArray const & array, const char* filename)
{
	FILE* f = fopen(filename, "w");
	for (size_t i = 0; i < array.size(); i++)
		fprintf(f, "%ld %g\n", i, array[i]);
	fclose(f);
}

void write_array(rArray const & grid, rArray const & array, const char* filename)
{
	FILE* f = fopen(filename, "w");
	for (size_t i = 0; i < array.size(); i++)
		fprintf(f, "%ld %g %g\n", i, grid[i], array[i]);
	fclose(f);
}

void write_array(cArray const & array, const char* filename)
{
	FILE* f = fopen(filename, "w");
	for (size_t i = 0; i < array.size(); i++)
		fprintf(f, "%ld %g %g\n", i, array[i].real(), array[i].imag());
	fclose(f);
}

void write_array(rArray const & grid, cArray const & array, const char* filename)
{
	FILE* f = fopen(filename, "w");
	for (size_t i = 0; i < array.size(); i++)
		fprintf(f, "%ld %g %g %g\n", i, grid[i], array[i].real(), array[i].imag());
	fclose(f);
}

void write_array(qArray const & array, const char* filename)
{
	FILE* f = fopen(filename, "w");
	for (size_t i = 0; i < array.size(); i++)
		fprintf(f, "%ld %Lg\n", i, array[i]);
	fclose(f);
}

// rArray logspace(double x0, double x1, unsigned N)
// {
// 	if (x0 <= 0 or x1 <= 0 or x1 < x0)
// 	{
// 		fprintf(stderr, "[logspace] It must be 0 < x1 <= x2 !\n");
// 		abort();
// 	}
// 	
// 	rArray grid(N);
// 	
// 	if (N == 1)
// 		grid[0] = x0;
// 	
// 	if (N > 1)
// 		for (unsigned i = 0; i < N; i++)
// 			grid[i] = x0 * pow(x1 / x0, i / (N - 1.));
// 	
// 	return grid;
// }

// rArray linspace(double x0, double x1, unsigned N)
// {
// 	rArray grid(N);
// 	
// 	if (N == 1)
// 		grid[0] = x0;
// 	
// 	if (N > 1)
// 		for (unsigned i = 0; i < N; i++)
// 			grid[i] = x0 + (x1 - x0) * i / (N - 1);
// 	
// 	return grid;
// }