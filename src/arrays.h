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

#ifndef HEX_ARRAYS
#define HEX_ARRAYS

#include <complex>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <typeinfo>
#include <vector>

#include <assert.h>

#ifndef NO_HDF
#include <H5Cpp.h>
#endif

#include "complex.h"
#include "misc.h"

template <typename NumberType> class Array;

/**
 * \brief Array shallow copy.
 */
template <typename NumberType> class ArrayView
{
	private:
		
		size_t N;
		NumberType * array;
		
	public:
		
		// alias
		typedef NumberType DataType;
		
		// empty constructor
		ArrayView() : N(0), array(nullptr) {}
		
		// construct view from Array non-const lvalue reference
		ArrayView(Array<NumberType> & a, size_t i = 0, size_t n = 0)
		{
			N = (n > 0) ? n : a.size();
			array = &a[0] + i;
		}
		
		// construct view from Array const lvalue reference
		ArrayView(Array<NumberType> const & a, size_t i = 0, size_t n = 0)
		{
			N = (n > 0) ? n : a.size();
			array = const_cast<NumberType*>(&a[0]) + i;
		}
		
		// destructor
		~ArrayView() {}
		
		//
		// assignments
		//
		
		virtual ArrayView<NumberType> & operator= (Array<NumberType> const & v)
		{
			if (v.size() != N)
				throw exception("[ArrayView::operator=] Cannot copy %ld elements to %ld fields!", N, v.size());
			
			for (size_t i = 0; i < N; i++)
				array[i] = v[i];
			
			return *this;
		}
		
		//
		// element-wise access (non-const)
		//
		
		virtual NumberType& operator[] (size_t i)
		{
		#ifdef NDEBUG
			return array[i];
		#else
			// bounds check
			if (i < N)
				return array[i];
			else
				throw exception("[ArrayView::operator[]] Index %ld out of bounds (size = %ld) !", i, N);
		#endif
		}
		
		//
		// element-wise access (const)
		//
		
		virtual NumberType const & operator[] (size_t i) const
		{
		#ifdef NDEBUG
			return array[i];
		#else
			if (i < N)
				return array[i];
			else
				throw exception("[ArrayView::operator[]] Index %ld out of bounds (size = %ld) !", i, N);
		#endif
		}
		
		// getters
		virtual size_t size() const { return N; }
		
		//
		// data pointer
		//
		
		virtual NumberType* data() { return array; }
		virtual const NumberType* data() const { return array; }
		
		//
		// STL-like iterator interface
		//
		
		typedef NumberType* iterator;
		typedef const NumberType* const_iterator;
		virtual iterator begin ()
				{ return array; }
		virtual const_iterator begin () const
				{ return array; }
		virtual iterator end ()
				{ return array + N; }
		virtual const_iterator end () const
				{ return array + N; }
		virtual NumberType & front (int i = 0)
				{ return *(array + i); }
		virtual NumberType const & front (int i = 0) const
				{ return *(array + i); }
		virtual NumberType & back (int i = 0)
				{ return *(array + N - 1 - i); }
		virtual NumberType const & back (int i = 0) const
				{ return *(array + N - 1 - i); }
		
		//
		// reduced arithmetic operators with other arrays
		//
		
		virtual ArrayView<NumberType>& operator += (Array<NumberType> const &  b)
		{
			// check if sizes match
 			assert(N == b.size());
			
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] += b[i];
			
			return *this;
		}
		
		virtual ArrayView<NumberType>& operator += (ArrayView<NumberType> const &  b)
		{
			// check if sizes match
 			assert(N == b.N);
			
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] += b.array[i];
			
			return *this;
		}
		
		virtual ArrayView<NumberType>& operator -= (Array<NumberType> const &  b)
		{
			// check if sizes match
 			assert(N == b.size());
			
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] -= b[i];
			
			return *this;
		}
		
		virtual ArrayView<NumberType>& operator -= (ArrayView<NumberType> const &  b)
		{
			// check if sizes match
 			assert(N == b.N);
			
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] -= b.array[i];
			
			return *this;
		}
		
		virtual ArrayView<NumberType>& operator /= (Array<NumberType> const &  b)
		{
			// check if sizes match
 			assert(N == b.size());
			
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] /= b[i];
			
			return *this;
		}
		
		virtual ArrayView<NumberType>& operator /= (ArrayView<NumberType> const &  b)
		{
			// check if sizes match
 			assert(N == b.N);
			
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] /= b.array[i];
			
			return *this;
		}
		
		// some other functions
		virtual void clear()
		{
			for (size_t i = 0; i < N; i++)
				array[i] = NumberType(0);
		}
};

/**
 * \brief A comfortable number array class.
 * 
 * Class Array is intended as a Hex's replacement for std::vector\<NumberType\>.
 * Properties:
 * - basic iterator interface (members Array::begin(), Array::end()).
 * - HDF5 interface (ability to save and load to/from HDF5 data files)
 * - a collection of overloaded arithmetic operators (sum of two arrays,
 *   difference, multiplication by a number etc.)
 */
template <typename NumberType> class Array : public ArrayView<NumberType>
{
	private:
		
		size_t N;
		NumberType * array;
		
	public:
		
		// alias
		typedef NumberType DataType;
		
		// inner product of two arrays
		template <typename NumberType1, typename NumberType2> friend auto operator | (
			Array<NumberType1> const & a, Array<NumberType2> const & b
		) -> decltype(NumberType1(0)*NumberType2(0));
				
		// default constructor, creates an empty array
		Array() : N(0), array(nullptr) {}
		
		// constructor, creates a length-n "x"-filled array
		Array(size_t n, NumberType x = 0) : N(n)
		{
			// reserve space
			array = new NumberType [N]();
					
			// set to zero
			for (size_t i = 0; i < N; i++)
				array[i] = x;
		}
		
		// constructor, copies a length-n "array
		Array(size_t n, NumberType* x) : N(n)
		{
			// reserve space
			array = new NumberType [N]();
					
			// set to zero
			for (size_t i = 0; i < N; i++)
				array[i] = x[i];
		}
		
		// copy constructor from Array const lvalue reference
		Array(Array<NumberType> const & a)
		{
			// reserve space
			N = a.N;
			array = new NumberType [N]();
	
			// run over the elements
			for (size_t i = 0; i < N; i++)
				array[i] = a.array[i];
		}
		
		// copy constructor from Array rvalue reference
		Array(Array<NumberType> && a)
		{
			// copy content
			N = a.N;
			array = a.array;
			
			// clear rvalue
			a.N = 0;
			a.array = nullptr;
		}
		
		// copy constructor from std::vector
		Array(const std::vector<NumberType>&  a)
		{
			// reserve space
			N = a.size();
			array = new NumberType [N]();
			
			// run over the elements
			for (size_t i = 0; i < N; i++)
				array[i] = a[i];
		}
		
		// copy constructor from initializer list
		Array(std::initializer_list<NumberType> a)
		{
			// reserve space
			N = a.end() - a.begin();
			array = new NumberType [N]();
			
			// run over the elements
			size_t i = 0;
			for (auto it = a.begin(); it != a.end(); it++)
				array[i++] = *it;
		}
		
		// destructor
		~Array()
		{
			if (array != nullptr)
				delete [] array;
		}
		
		//
		// storage size
		//
		
		size_t size() const { return N; }
		void resize (size_t n)
		{
			NumberType * new_array = new NumberType [n]();
			for (size_t i = 0; i < n; i++)
				new_array[i] = (i < N) ? array[i] : NumberType(0);
			delete [] array;
			N = n;
			array = new_array;
		}
		
		//
		// element-wise access (non-const)
		//
		
		inline NumberType& operator[] (size_t i)
		{
		#ifdef NDEBUG
			return array[i];
		#else
			// bounds check
			if (i < N)
				return array[i];
			else
				throw exception("[Array::operator[]] Index %ld out of bounds (size = %ld) !", i, N);
		#endif
		}
		
		//
		// element-wise access (const)
		//
		
		inline NumberType const & operator[] (size_t i) const
		{
		#ifdef NDEBUG
			return array[i];
		#else
			if (i < N)
				return array[i];
			else
				throw exception("[Array::operator[]] Index %ld out of bounds (size = %ld) !", i, N);
		#endif
		}
		
		//
		// data pointer
		//
		
		NumberType* data() { return array; }
		const NumberType* data() const { return array; }
		
		//
		// STL-like iterator interface
		//
		
		typedef NumberType* iterator;
		typedef const NumberType* const_iterator;
		iterator begin()
				{ return array; }
		const_iterator begin() const
				{ return array; }
		iterator end()
				{ return array + N; }
		const_iterator end() const
				{ return array + N; }
		NumberType & front(int i = 0)
				{ return *(array + i); }
		NumberType const & front(int i = 0) const
				{ return *(array + i); }
		NumberType & back(int i = 0)
				{ return *(array + N - 1 - i); }
		NumberType const & back(int i = 0) const
				{ return *(array + N - 1 - i); }
		
		void push_back(NumberType a)
		{
			// not very efficient... FIXME
			
			NumberType* new_array = new NumberType [N + 1]();
			for (size_t i = 0; i < N; i++)
				new_array[i] = array[i];
			new_array[N] = a;
			N++;
			delete [] array;
			array = new_array;
		}
		
		NumberType pop_back()
		{
			if (N > 0)
				return *(array + (--N));
			else
				throw exception ("Array has no element to pop!");
		}
		
		template <class InputIterator> void append (
			InputIterator first, InputIterator last
		) {
			NumberType* new_array = new NumberType [N + last - first]();
			for (size_t i = 0; i < N; i++)
				new_array[i] = array[i];
			for (InputIterator it = first; it != last; it++)
				new_array[N + it - first] = *it;
			N += last - first;
			delete [] array;
			array = new_array;
		}
		
		bool empty() const
		{
			return N == 0;
		}
		
		//
		// assignment operators
		//
		
		Array<NumberType>& operator = (Array<NumberType> const &  b)
		{
			// if we already have some allocated space, check its size,
			// so that we do not free it uselessly
			if (array != nullptr and N != b.N)
			{
				delete [] array;
				array = nullptr;
			}
			
			// set the new dimension
			N = b.N;
			
			// if necessary, reserve space
			if (array == nullptr)
				array = new NumberType [N]();
			
			// run over the elements
			for (size_t i = 0; i < N; i++)
				array[i] = b.array[i];
			
			return *this;
		}
		
		Array<NumberType>& operator = (Array<NumberType> &&  b)
		{
			// if we already have some allocated space, check its size,
			// so that we do not free it uselessly
			if (array != nullptr)
			{
				delete [] array;
				array = nullptr;
			}
			
			// move content
			N = b.N;
			array = b.array;
			
			// clear rvalue
			b.N = 0;
			b.array = nullptr;
			
			return *this;
		}
		
		//
		// reduced arithmetic operators with other arrays
		//
		
		Array<NumberType>& operator += (Array<NumberType> const &  b)
		{
			// check if sizes match
			assert(N == b.N);
			
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] += b.array[i];

			return *this;
		}
		
		Array<NumberType>& operator -= (Array<NumberType> const &  b)
		{
			// check if sizes match
			assert(N == b.N);
			
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] -= b.array[i];
			
			// return
			return *this;
		}
		
		Array<NumberType>& operator *= (Array<NumberType> const &  b)
		{
			// check if sizes match
			assert(N == b.N);
			
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] *= b.array[i];
			
			// return
			return *this;
		}
		
		Array<NumberType>& operator /= (Array<NumberType> const &  b)
		{
			// check size
			assert(b.size() == N);
			
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] /= b[i];
			
			// return
			return *this;
		}
		
		//
		// reduced arithmetic operators with complex numbers
		//
		
		template <typename NumberType2> Array<NumberType>& operator += (NumberType2 z)
		{
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] += NumberType(z);

			return *this;
		}
		
		template <typename NumberType2> Array<NumberType>& operator -= (NumberType2 z)
		{
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] -= NumberType(z);
			
			// return
			return *this;
		}
		
		template <typename NumberType2> Array<NumberType>& operator *= (NumberType2 z)
		{
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] *= NumberType(z);
			
			// return
			return *this;
		}
		
		template <typename NumberType2> Array<NumberType>& operator /= (NumberType2 z)
		{
			// run over elements
			for (size_t i = 0; i < N; i++)
				array[i] /= NumberType(z);
			
			// return
			return *this;
		}
		
		// complex conjugate
		Array<NumberType> conj() const
		{
			Array<NumberType> c = *this;
			for (size_t i = 0; i < N; i++)
			{
				Complex z = c.array[i];
				c.array[i] = Complex(z.real(), -z.imag());
			}
			return c;
		}
		
		// compute usual 2-norm
		double norm() const
		{
			double ret = 0.;
			for (size_t i = 0; i < N; i++)
			{
				Complex z = array[i];
				ret += z.real() * z.real() + z.imag() * z.imag();
			}
			return sqrt(ret);
		}
		
		// apply user transformation
		template <class Functor> auto transform(Functor f) -> Array<decltype(f(NumberType(0)))>
		{
			Array<decltype(f(NumberType(0)))> c(N);
			for (size_t i = 0; i < N; i++)
				c[i] = f(array[i]);
			return c;
		}
		
		// rrturn subarray
		Array<NumberType> slice(size_t left, size_t right)
		{
			Array<NumberType> c(right - left);
			NumberType * ptr_c = &c[0];
			
			for (size_t i = left; i < right; i++)
				*ptr_c++ = array[i];
			
			return c;
		}
		
#ifndef NO_HDF
		/**
		 * Save array to HDF file.
		 * \param name Filename.
		 * \param compress Whether to apply a trivial compression (contract the repeated zeros).
		 * \param consec Minimal consecutive occurences for compression.
		 */
		bool hdfsave(const char* name, bool compress = false, int consec = 10) const
		{
			H5::Exception::dontPrint();
			
			// compressed array
			Array<NumberType> carray;
			
			// zero blocks
			Array<int> zero_blocks;
			
			// consecutive zeros counter
			int zero_counter = 0;
			
			if (compress)
			{
				// analyze: find compressible segments
				for (size_t i = 0; i < N; i++)
				{
					if (array[i] == 0.)
					{
						// another consecutive zero
						zero_counter++;
					}
					else if (zero_counter >= consec)
					{
						// end of large zero block -> compress
						zero_blocks.push_back(i-zero_counter);
						zero_blocks.push_back(i);
						zero_counter = 0;
					}
					else
					{
						// end of tiny zero block -> do not bother with compression
						zero_counter = 0;
					}
					
				}
				if (zero_counter >= 10)
				{
					zero_blocks.push_back(N-zero_counter);
					zero_blocks.push_back(N);
				}
				
				// invert selection: get non-zero blocks
				Array<int> nonzero_blocks;
				nonzero_blocks.push_back(0);
				nonzero_blocks.append(zero_blocks.begin(), zero_blocks.end());
				nonzero_blocks.push_back(N);
				
				// compress: copy only nonzero elements
				for (size_t iblock = 0; iblock < nonzero_blocks.size()/2; iblock++)
				{
					int start = nonzero_blocks[2*iblock];
					int end = nonzero_blocks[2*iblock+1];
					carray.append (
						array + start,
						array + end
					);
				}
			}
			
			// save to HDF file
			try
			{
				H5::H5File h5file(name, H5F_ACC_TRUNC);
				
				// save data as an interleaved array
				if (compress)
				{
					hsize_t length1 = carray.size() * sizeof(NumberType) / sizeof(double);
					
					H5::DataSpace dspc1(/*rank*/1, &length1);
					H5::IntType dtype1(H5::PredType::NATIVE_DOUBLE);
					H5::DataSet dset1 = h5file.createDataSet("array", dtype1, dspc1);
					dset1.write(&carray[0], H5::PredType::NATIVE_DOUBLE);
					
					hsize_t length2 = zero_blocks.size();
					
					// make zero_blocks index doubles, not elements
					for (int& z : zero_blocks)
						z *= sizeof(NumberType)/sizeof(double);
					
					H5::DataSpace dspc2(/*rank*/1, &length2);
					H5::IntType dtype2(H5::PredType::NATIVE_INT);
					H5::DataSet dset2 = h5file.createDataSet("zero_blocks", dtype2, dspc2);
					dset2.write(&zero_blocks[0], H5::PredType::NATIVE_INT);
				}
				else
				{
					hsize_t length = N * sizeof(NumberType) / sizeof(double);
					
					H5::DataSpace dspc(/*rank*/1, &length);
					H5::IntType dtype(H5::PredType::NATIVE_DOUBLE);
					H5::DataSet dset = h5file.createDataSet("array", dtype, dspc);
					dset.write(array, H5::PredType::NATIVE_DOUBLE);
				}
				
				return true;
			}
			catch (H5::FileIException err)
			{
				return false;
			}
		}
		
		/**
		 * Load array from HDF file.
		 * \param name Filename.
		 */
		bool hdfload(const char* name)
		{
			H5::Exception::dontPrint();
			
			// number of nonzero elements
			int nnz;
			
			// array of nonzero elements
			Array<NumberType> nnz_array;
			
			// starting and (one past) ending elements of the elided zero blocks
			Array<int> zero_blocks;
			
			try
			{
				// open file
				H5::H5File h5file(name, H5F_ACC_RDONLY);
				
				// load data array (non-zero elements)
				H5::DataSet dset = h5file.openDataSet("array");
				H5::DataSpace dspc = dset.getSpace();
				nnz = dspc.getSimpleExtentNpoints() * sizeof(double) / sizeof(NumberType);
				nnz_array.resize(nnz);
				dset.read (
					reinterpret_cast<double*>(&nnz_array[0]),
					H5::PredType::NATIVE_DOUBLE,
					dspc,
					dspc
				);
				
				// load compression information (if any)
				H5::DataSet dsetz;
				H5::DataSpace dspcz;
				int n;
				try {
					dsetz = h5file.openDataSet("zero_blocks");
					dspcz = dsetz.getSpace();
					n = dspcz.getSimpleExtentNpoints();
				} catch (H5::FileIException err) {
					n = 0;
				}
				
				if (n == 0)
				{
					// no compression data
					// -> move loaded data to internal storage and return
					*this = std::move(nnz_array);
					return true;
				}
				else
				{
					// there are some compression data
					// -> load information about the zero blocks
					zero_blocks.resize(n);
					dsetz.read (
						reinterpret_cast<int*>(&zero_blocks[0]),
						H5::PredType::NATIVE_INT,
						dspcz,
						dspcz
					);
				}
			}
			catch (H5::FileIException err)
			{
				return false;
			}
			
			// make zero_blocks index elements, not doubles
			for (int& z : zero_blocks)
				z /= sizeof(NumberType)/sizeof(double);
			
			// remove previous data
			if (array != nullptr)
			{
				delete [] array;
				array = nullptr;
				N = 0;
			}
			
			// compute final size
			N = nnz;
			for (size_t i = 0; i < zero_blocks.size()/2; i++)
				N += zero_blocks[2*i+1] - zero_blocks[2*i];
			
			// resize and clean internal storage
			array = new NumberType[N];
			memset(array, 0, N * sizeof(NumberType));
			
			// copy nonzero chunks
			int this_end = 0;	// index of last updated element in "this"
			int load_end = 0;	// index of last used element in "nnz_array"
			for (size_t i = 0; i < zero_blocks.size()/2; i++)
			{
				int zero_start = zero_blocks[2*i];
				int zero_end = zero_blocks[2*i+1];
				
				// append nonzero data before this zero block
				memcpy (
					array + this_end,
					&nnz_array[0] + load_end,
					(zero_start - this_end) * sizeof(NumberType)
				);
				
				// move cursors
				load_end += zero_start - this_end;
				this_end  = zero_end;
			}
			
			// append remaining data
			memcpy (
				array + this_end,
				&nnz_array[0] + load_end,
				(N - this_end) * sizeof(NumberType)
			);
			
			return true;
		}
#endif
};

// scalar product of two arrays.
template <typename NumberType1, typename NumberType2> auto operator | (
	Array<NumberType1> const & a, Array<NumberType2> const & b
) -> decltype(NumberType1(0)*NumberType2(0)) {
	// store size
	size_t N = a.N;
	
	// check if sizes match
	assert(N == b.N);
	
	// the scalar product
	decltype(NumberType1(0)*NumberType2(0)) result = 0;
	
	// iterators
	const NumberType1* const __a = &a[0];
	const NumberType2* const __b = &b[0];
	
	// sum the products
	for (size_t i = 0; i < N; i++)
		result += __a[i] * __b[i];
	
	return result;
}

// arithmetic operators Array & Array
template <typename NumberType1, typename NumberType2> auto operator + (
	Array<NumberType1> const & a, Array<NumberType2> const & b
) -> Array<decltype(NumberType1(0) + NumberType2(0))>
{
	Array<decltype(NumberType1(0) + NumberType2(0))> c = a;
	return c += b;
}

template <typename NumberType1, typename NumberType2> auto operator - (
	Array<NumberType1> const & a, Array<NumberType2> const & b
) -> Array<decltype(NumberType1(0) - NumberType2(0))>
{
	Array<decltype(NumberType1(0) - NumberType2(0))> c = a;
	return c -= b;
}

template <typename NumberType1, typename NumberType2> auto operator * (
	Array<NumberType1> const & a, Array<NumberType2> const & b
) -> Array<decltype(NumberType1(0) * NumberType2(0))>
{
	Array<decltype(NumberType1(0) * NumberType2(0))> c = a;
	return c *= b;
}

template <typename NumberType1, typename NumberType2> auto operator / (
	Array<NumberType1> const & a, Array<NumberType2> const & b
) -> Array<decltype(NumberType1(0) / NumberType2(1))>
{
	Array<decltype(NumberType1(0) / NumberType2(1))> c = a;
	return c /= b;
}

// arithmetic operators Array & Number
template <typename NumberType1, typename NumberType2> auto operator + (
	Array<NumberType1> const & a, NumberType2 z
) -> Array<decltype(NumberType1(0) + NumberType2(0))>
{
	Array<decltype(NumberType1(0) + NumberType2(0))> c = a;
	return c += z;
}

template <typename NumberType1, typename NumberType2> auto operator - (
	Array<NumberType1> const & a, NumberType2 z
) -> Array<decltype(NumberType1(0) - NumberType2(0))>
{
	Array<decltype(NumberType1(0) - NumberType2(0))> c = a;
	return c -= z;
}

template <typename NumberType1, typename NumberType2> auto operator * (
	Array<NumberType1> const & a, NumberType2 z
) -> Array<decltype(NumberType1(0) * NumberType2(0))>
{
	Array<decltype(NumberType1(0) * NumberType2(0))> c = a;
	return c *= z;
}

template <typename NumberType1, typename NumberType2> auto operator * (
	NumberType1 z, Array<NumberType2> const & b
) -> Array<decltype(NumberType1(0) * NumberType2(0))>
{
	Array<decltype(NumberType1(0) * NumberType2(0))> c = b;
	return c *= z;
}

template <typename NumberType1, typename NumberType2> auto operator / (
	Array<NumberType1> const & a, NumberType2 z
) -> Array<decltype(NumberType1(0) / NumberType2(1))>
{
	Array<decltype(NumberType1(0) / NumberType2(1))> c = a;
	return c /= z;
}

// other vectorized functions
template <typename NumberType> Array<NumberType> sqrt (Array<NumberType> const & A)
{
	size_t N = A.size();
	Array<NumberType> B (N);

	for (size_t i = 0; i < N; i++)
		B[i] = sqrt(A[i]);

	return B;
}

inline Array<double> sqrabs (Array<Complex> const & A)
{
	size_t N = A.size();
	Array<double> B (N);

	for (size_t i = 0; i < N; i++)
		B[i] = sqrabs(A[i]);

	return B;
}

// output to text stream.
template <typename NumberType> std::ostream & operator << (std::ostream & out, Array<NumberType> const & a)
{
	out << "[";
	for (size_t i = 0; i < a.size(); i++)
	{
		if (i == 0)
			out << a[i];
		else
			out << "," << a[i];
	}
	out << "]";
	
	return out;
}

/**
 * Generate uniform grid
 * \param start Left boundary and first sample for "samples" > 0.
 * \param end Right boundary and last sample for "samples" > 1.
 * \param samples Sample count.
 */
template <typename T> Array<T> linspace(T start, T end, unsigned samples)
{
	Array<T> space(samples);
	
	if (samples == 0)
		return space;
	
	if (samples == 1)
	{
		space[0] = start;
		return space;
	}
	
	for (unsigned i = 0; i < samples; i++)
		space[i] = start + (end - start) * T(i) / T(samples - 1);
	return space;
}

/**
 * Generate logarithmic grid
 * \param x0 Left boundary and first sample for "samples" > 0.
 * \param x1 Right boundary and last sample for "samples" > 1.
 * \param N Sample count.
 */
template <typename T> Array<T> logspace(T x0, T x1, size_t N)
{
	if (x0 <= 0 or x1 <= 0 or x1 < x0)
	{
		fprintf(stderr, "[logspace] It must be 0 < x1 <= x2 !\n");
		abort();
	}
	
	Array<T> grid(N);
	
	if (N == 1)
		grid[0] = x0;
	
	if (N > 1)
		for (unsigned i = 0; i < N; i++)
			grid[i] = x0 * pow(x1 / x0, i / T(N - 1));
	
	return grid;
}

/**
 * Write array to standard output. Array will be written as a single column.
 * \param array The array to write.
 */
template <typename NumberType> void write_array(Array<NumberType> const & array)
{
	for (size_t i = 0; i < array.size(); i++) 
	{
		if (typeid(NumberType) == typeid(double))
		{
			std::cout << array[i] << std::endl;
		}
		else if (typeid(NumberType) == typeid(Complex))
		{
			std::cout << Complex(array[i]).real() << "\t" << Complex(array[i]).imag() << std::endl;
		}
		else
		{
			std::cerr << "Don't know how to write datatype with typeid " << typeid(NumberType).name() << std::endl;
		}
	}
}

/**
 * Write array to file. Array will be written as a single column into
 * an ASCII file.
 * \param array The array to write.
 * \param filename Name of the file to create/overwrite.
 */
template <typename NumberType> void write_array(Array<NumberType> const & array, const char* filename)
{
	std::ofstream fout(filename);
	for (size_t i = 0; i < array.size(); i++) 
	{
		switch (typeid(NumberType))
		{
			case typeid (double):
				fout << array[i] << std::endl;
				break;
			case typeid (Complex):
				fout << array[i].real() << "\t" << array[i].imag() << std::endl;
				break;
			default:
				std::cerr << "Don't know how to write datatype with typeid " << typeid(NumberType).name() << std::endl;
				return;
		}
	}
}

/**
 * Write array to file. There will be two columns in the resulting ASCII file.
 * One contains the elements from first supplied array (\c grid) and other
 * containing elements from the second array (\c array). Useful for named
 * plots.
 * \param grid Data labels (data for first column).
 * \param array The array to write (data for second column).
 * \param filename Name of the file to create/overwrite.
 */
template <typename NumberType1, typename NumberType2> void write_array(
	Array<NumberType1> const & grid, Array<NumberType2> const & array, const char* filename
)
{
	std::ofstream fout(filename);
	for (size_t i = 0; i < grid.size(); i++)
	{
		if (typeid(NumberType1) == typeid (double))
		{
				fout << grid[i] << "\t";
		}
		else if (typeid(NumberType1) == typeid (Complex))
		{
			fout << Complex(grid[i]).real() << "\t" << Complex(grid[i]).imag() << "\t";
		}
		else
		{
			std::cerr << "Don't know how to write datatype with typeid " << typeid(NumberType1).name() << std::endl;
			return;
		}
		
		if (typeid(NumberType2) == typeid (double))
		{
				fout << array[i] << std::endl;
		}
		else if (typeid(NumberType2) == typeid (Complex))
		{
			fout << Complex(array[i]).real() << "\t" << Complex(array[i]).imag() << std::endl;
		}
		else
		{
			std::cerr << "Don't know how to write datatype with typeid " << typeid(NumberType2).name() << std::endl;
			return;
		}
	}
}

template <typename Fetcher> bool write_1D_data (size_t m, const char* filename, Fetcher fetch)
{
	FILE* f = fopen(filename, "w");
	if (f == 0)
		return false;
	
	for (size_t i = 0; i < m; i++)
		fprintf(f, "%g\n", fetch(i));
	
	fclose(f);
	return true;
}

/**
 * Write 2D data to file. To allow maximum flexibility, only extensions
 * of the data are passed to the function and a functor that will be
 * repeatedly called with coordinate pair for new data element.
 * \param m Row count.
 * \param n Column count.
 * \param filename Filename of the file to create/overwrite.
 * \param fetch Functor with interface
 *        \code
 *             double operator() (size_t, size_t);
 *        \endcode
 * \return Write success indicator (\c false for failure).
 */
template <class Fetcher> bool write_2D_data(size_t m, size_t n, const char* filename, Fetcher fetch)
{
	FILE* f = fopen(filename, "w");
	if (f == 0)
		return false;
	
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
			fprintf(f, "%g ", fetch(i,j));
		fprintf(f, "\n");
	}
	fclose(f);
	
	return true;
}

// aliases
typedef Array<double>		rArray;
typedef Array<Complex>		cArray;
typedef Array<long double>	qArray;
typedef Array<rArray> rArrays;
typedef Array<cArray> cArrays;

typedef ArrayView<double> rArrayView;
typedef ArrayView<Complex> cArrayView;
typedef ArrayView<long double> qArrayView;

/**
 * Variadic template recurrence starter. For documentation of the function
 * itself see the other "concatenate".
 */
inline rArray concatenate()
{
	return rArray (0);
}

/**
 * Concatenate several arrays. The function template uses variadic templates
 * feature of C++, so the number of subarrays to concatenate is completely
 * arbitrary. It should be noted, though, that long concatenation list may
 * slow down the template instantiation during compilation and hence the
 * overall compilation time.
 * \param v1 First array.
 * \param ...p All other arrays.
 */
template <typename ...Params> rArray concatenate(rArray v1, Params ...p)
{
	if (sizeof...(p) == 0)
	{
		return v1;
	}
	else
	{
		rArray v2 = concatenate(p...);
		rArray v (v1.size() + v2.size());
		for (size_t i = 0; i < v1.size(); i++)
			v[i] = v1[i];
		for (size_t i = 0; i < v2.size(); i++)
			v[i + v1.size()] = v2[i];
		return v;
	}	
}

/**
 * Load / save array from a HDF5 file.
 */
#ifndef NO_HDF
bool load_array(rArray& vec, const char* name, double* pdelta = 0);
bool load_array(cArray& vec, const char* name);
bool save_array(rArray const & vec, const char* name, const double * const pdelta = 0);
bool save_array(cArray const & vec, const char* name);
#endif
/**
 * Write array to a text file.
 */
void write_array(const std::map<unsigned long long, Complex>& m, const char* filename);
void write_array(rArray const & grid, rArray const & array, const char* filename);
void write_array(rArray const & array, const char* filename);
void write_array(cArray const & array, const char* filename);
void write_array(rArray const & grid, cArray const & array, const char* filename);
void write_array(qArray const & array, const char* filename);

// return absolute values
inline rArray abs (cArray const &u)
{
	rArray v(u.size());
	
	auto iu = u.begin();
	auto iv = v.begin();
	
	while (iu != u.end())
		*(iv++) = abs(*(iu++));
	
	return v;
}

template <typename NumberType> NumberType min (Array<NumberType> const & a)
{
	NumberType z = a.front();
	for (NumberType const * it = a.begin(); it != a.end(); it++)
		if (*it < z)
			z = *it;
	return z;
}

template <typename NumberType> NumberType max (Array<NumberType> const & a)
{
	NumberType z = a.front();
	for (NumberType const * it = a.begin(); it != a.end(); it++)
		if (*it > z)
			z = *it;
	return z;
}

inline rArrays abs (cArrays const &u)
{
	rArrays v(u.size());
	
	auto iu = u.begin();
	auto iv = v.begin();
	
	while (iu != u.end())
		*(iv++) = abs(*(iu++));
	
	return v;
}

// return per-element power
template <typename T> Array<T> pow(Array<T> const &u, double e)
{
	Array<T> v(u.size());
	
	auto iu = u.begin();
	auto iv = v.begin();
	
	while (iu != u.end())
		*(iv++) = pow(*(iu++), e);
	
	return v;
}

// boolean aggregation
bool all(Array<bool> v);
bool any(Array<bool> v);

// summation
template <typename T> T sum(Array<T> v)
{
	return std::accumulate(v.begin(), v.end(), T(0));
}

// summation of nested arrays
template <typename T> Array<T> sums(Array<Array<T>> v)
{
	if (v.size() == 0)
		return Array<T>();	// empty array
	
	return std::accumulate(
		v.begin(),
		v.end(),
		Array<T> (v[0].size()),
		[](Array<T> a, Array<T> b) -> Array<T> {
			return a + b;
		}
	);
}

/**
 * Comparison of an array and a number.
 * \param u Array.
 * \param x Number.
 * \return Vector of bools for element-wise comparisons.
 */
template <typename T>
Array<bool> operator == (Array<T> u, T x)
{
	Array<bool> v(u.size());
	for (size_t i  = 0; i < u.size(); i++)
		v[i] = (u[i] == x);
	return v;
}

/**
 * Evaluation of a function over a grid.
 * \param f Functor to be evaluated using the operator() (double) interface.
 * \param grid Array-like type containing the evaluation points.
 * \param vals Array-like type to hold evaluated function on return. It is
 *             required that enough space is reserved, at least grid.size().
 */
template <typename TFunctor, typename TArray>
void eval(TFunctor f, TArray grid, TArray& vals)
{
	size_t N = grid.size();
	assert(N == vals.size());
	
	for (size_t i = 0; i < N; i++)
		vals[i] = f(grid[i]);
}

/**
 * Sum two indexed arrays if their (sorted) indices perfectly match.
 * Or, generally, create output array with elements that are sum of corresponding
 * elements of both arrays or equal to a single element in one array, if
 * that its index doesn't have a counterpart in the other array.
 * \param idx1 Sorted (!) indices of the first array.
 * \param arr1 Merge TO array.
 * \param idx2 Sorted (!) indices of the second array.
 * \param arr2 Merge FROM array.
 */
template <typename Tidx, typename Tval> void merge (
	Array<Tidx>       & idx1, Array<Tval>       & arr1,
	Array<Tidx> const & idx2, Array<Tval> const & arr2
){
	// positions in arrays
	size_t i1 = 0;
	size_t i2 = 0;
	
	// output arrays
	Array<Tidx> idx;
	Array<Tval> arr;
	
	// while there is anything to merge
	while (i1 < arr1.size() and i2 < arr2.size())
	{
		if (idx1[i1] == idx2[i2])
		{
			idx.push_back(idx1[i1]);
			arr.push_back(arr1[i1] + arr2[i2]);
			i1++;
			i2++;
		}
		else if (idx1[i1] < idx2[i2])
		{
			idx.push_back(idx1[i1]);
			arr.push_back(arr1[i1]);
			i1++;
		}
		else /* idx1[i2] > idx2[i2] */
		{
			idx.push_back(idx2[i2]);
			arr.push_back(arr2[i2]);
			i2++;
		}
	}
	
	// the rest will be done by a single copy
	if (i1 == arr1.size() and i2 < arr2.size())
	{
		idx.append(idx2.begin() + i2, idx2.end());
		arr.append(arr2.begin() + i2, arr2.end());
	}
	
	// copy to the first pair
	idx1 = idx;
	arr1 = arr;
}

#endif
