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

#ifndef HEX_ARRAYS_H
#define HEX_ARRAYS_H

#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <numeric>
#include <typeinfo>
#include <type_traits>
#include <vector>

#include <assert.h>

#ifndef NO_HDF
#include "hdffile.h"
#endif

#include "complex.h"
#include "misc.h"

class HDFFile;
template <typename DataType> class Array;
template <typename NumberType> class NumberArray;

/**
 * @brief Array shallow copy.
 */
template <typename NumberType> class ArrayView
{
private:
    
    size_t N;
    NumberType * array;
    
public:
    
    // alias
    typedef NumberType DataType;
    typedef NumberType * iterator;
    typedef NumberType const * const_iterator;
    
    // empty constructor
    ArrayView()
        : N(0), array(nullptr) {}
    
    // copy constructor
    ArrayView(ArrayView const & a, size_t i = 0, size_t n = 0)
        : N(n == 0 ? a.N : n), array(a.array + i) { assert(i + N <= a.N); }
    
    // construct view from Array non-const lvalue reference
    ArrayView(Array<NumberType> & a, size_t i = 0, size_t n = 0)
        : N(n == 0 ? a.size() : n), array(a.data() + i) { assert(i + N <= a.size()); }
    
    // construct view from Array const lvalue reference
    ArrayView(Array<NumberType> const & a, size_t i = 0, size_t n = 0)
    {
        N = (n > 0) ? n : a.size();
        array = const_cast<NumberType*>(a.data()) + i;
    }
    
    // construct from consecutive memory segment
    ArrayView(NumberType const * i, NumberType const * j)
    {
        N = j - i;
        array = const_cast<NumberType*>(&(*i));
    }
    
    // destructor
    ~ArrayView() {}
    
    //
    // assignments
    //
    
    ArrayView<NumberType> & operator= (ArrayView<NumberType> const & v)
    {
        if (v.size() != size())
            throw exception("[ArrayView::operator=] Cannot copy %ld elements to %ld fields!", v.size(), N);
        
        for (size_t i = 0; i < size(); i++)
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
        if (i < size())
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
        if (i < size())
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
    virtual NumberType const * data() const { return array; }
    
    //
    // STL-like iterator interface
    //
    
    virtual iterator begin ()
    { return array; }
    virtual const_iterator begin () const
    { return array; }
    virtual iterator end ()
    { return array + size(); }
    virtual const_iterator end () const
    { return array + size(); }
    virtual NumberType & front (int i = 0)
    { return *(array + i); }
    virtual NumberType const & front (int i = 0) const
    { return *(array + i); }
    virtual NumberType & back (int i = 0)
    { return *(array + size() - 1 - i); }
    virtual NumberType const & back (int i = 0) const
    { return *(array + size() - 1 - i); }
    
    // some other functions
    virtual void fill(NumberType x)
    {
        for (NumberType & y : *this)
            y = x;
    }
    
    bool empty() const { return size() == 0; }
    
    // 2-norm, defined only for scalar NumberType
    template <class = typename std::enable_if<is_scalar<NumberType>::value>> double norm() const
    {
        double sqrnorm = 0.;
        for (NumberType const & x : *this)
            sqrnorm += sqrabs(x);
        return sqrt(sqrnorm);
    }
};

/**
 * @brief A comfortable data array class.
 * 
 * Class Array is intended as a Hex's replacement for std::vector\<NumberType\>.
 * Properties:
 * - basic iterator interface (members Array::begin(), Array::end()).
 */
template <typename DataType> class Array : public ArrayView<DataType>
{
    
private:
    
    size_t N;
    DataType * array;
    
public:
    
    // default constructor, creates an empty array
    Array() : N(0), array(nullptr) {}
    
    // constructor, creates a length-n "x"-filled array
    Array(size_t n, DataType x = 0) : N(n)
    {
        // reserve space
        array = new DataType [N]();
        
        // set to zero
        for (size_t i = 0; i < N; i++)
            array[i] = x;
    }
    
    // constructor, copies a length-n "array
    Array(size_t n, DataType* x) : N(n)
    {
        // reserve space
        array = new DataType [N]();
        
        // set to zero
        for (size_t i = 0; i < N; i++)
            array[i] = x[i];
    }
    
    // copy constructor from Array const lvalue reference
    Array(Array<DataType> const & a)
    {
        // reserve space
        N = a.N;
        array = new DataType [N]();
        
        // run over the elements
        for (size_t i = 0; i < N; i++)
            array[i] = a.array[i];
    }
    
    // copy constructor from ArrayView const lvalue reference
    Array(ArrayView<DataType> const & a)
    {
        // reserve space
        N = a.size();
        array = new DataType [N]();
        
        // run over the elements
        for (size_t i = 0; i < N; i++)
            array[i] = a.data()[i];
    }
    
    // copy constructor from Array rvalue reference
    Array(Array<DataType> && a)
    {
        // copy content
        N = a.N;
        array = a.array;
        
        // clear rvalue
        a.N = 0;
        a.array = nullptr;
    }
    
    // copy constructor from std::vector
    Array(std::vector<DataType> const & a)
    {
        // reserve space
        N = a.size();
        array = new DataType [N]();
        
        // run over the elements
        for (size_t i = 0; i < N; i++)
            array[i] = a[i];
    }
    
    // copy constructor from initializer list
    Array(std::initializer_list<DataType> a)
    {
        // reserve space
        N = a.end() - a.begin();
        array = new DataType [N]();
        
        // run over the elements
        size_t i = 0;
        for (auto it = a.begin(); it != a.end(); it++)
            array[i++] = *it;
    }
    
    // copy constructor from two forward iterators
    template <typename ForwardIterator> Array(ForwardIterator i, ForwardIterator j)
    {
        // compute size
        N = 0;
        for (ForwardIterator k = i; k != j; k++)
            N++;
        
        // reserve space
        array = new DataType [N]();
        
        // run over the elements
        size_t n = 0;
        for (ForwardIterator k = i; k != j; k++)
            array[n++] = *k;
    }
    
    // destructor
    ~Array()
    {
        if (array != nullptr)
        {
            delete [] array;
            array = nullptr;
        }
    }
    
    //
    // storage size
    //
    
    size_t size() const { return N; }
    size_t resize (size_t n)
    {
        if (n == 0)
        {
            if (array != nullptr)
            {
                delete [] array;
                N = 0;
            }
            return 0;
        }
        
        DataType * new_array = new DataType [n]();
        for (size_t i = 0; i < n; i++)
            new_array[i] = (i < N) ? array[i] : DataType(0);
        delete [] array;
        N = n;
        array = new_array;
        return N;
    }
    
    //
    // element-wise access (non-const)
    //
    
    inline DataType& operator[] (size_t i)
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
    
    inline DataType const & operator[] (size_t i) const
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
    
    DataType* data() { return array; }
    DataType const * data() const { return array; }
    
    //
    // STL-like iterator interface
    //
    
    typedef DataType* iterator;
    typedef const DataType* const_iterator;
    iterator begin()
    { return array; }
    const_iterator begin() const
    { return array; }
    iterator end()
    { return array + N; }
    const_iterator end() const
    { return array + N; }
    DataType & front(int i = 0)
    { return *(array + i); }
    DataType const & front(int i = 0) const
    { return *(array + i); }
    DataType & back(int i = 0)
    { return *(array + N - 1 - i); }
    DataType const & back(int i = 0) const
    { return *(array + N - 1 - i); }
    
    void push_back(DataType const & a)
    {
        // not very efficient... FIXME
        
        DataType* new_array = new DataType [N + 1]();
        for (size_t i = 0; i < N; i++)
            new_array[i] = std::move(array[i]);
        new_array[N] = a;
        N++;
        delete [] array;
        array = new_array;
    }
    
    DataType pop_back()
    {
        if (N > 0)
            return *(array + (--N));
        else
            throw exception ("Array has no element to pop!");
    }
    
    template <class InputIterator> void append (
        InputIterator first, InputIterator last
    ) {
        DataType* new_array = new DataType [N + last - first]();
        for (size_t i = 0; i < N; i++)
            new_array[i] = array[i];
        for (InputIterator it = first; it != last; it++)
            new_array[N + it - first] = *it;
        N += last - first;
        delete [] array;
        array = new_array;
    }
    
    void insert (iterator it, DataType x)
    {
        // create new array (one element longer)
        DataType* new_array = new DataType [N + 1];
        
        // copy everything to the new location
        for (int i = 0; i < it - array; i++)
            new_array[i] = std::move(array[i]);
        
        // insert new element
        *(new_array + (it - array)) = std::move(x);
        
        // copy the rest
        for (int i = it - array; i < (int)N; i++)
            new_array[i+1] = std::move(array[i]);
        
        // change pointers
        delete [] array;
        N++;
        array = new_array;
    }
    
    bool empty() const
    {
        return N == 0;
    }
    
    //
    // assignment operators
    //
    
    Array<DataType>& operator = (Array<DataType> const &  b)
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
            array = new DataType [N]();
        
        // run over the elements
        for (size_t i = 0; i < N; i++)
            array[i] = b.array[i];
        
        return *this;
    }
    
    Array<DataType>& operator = (Array<DataType> &&  b)
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
};

/**
 * @brief A comfortable number array class.
 * 
 * Class NumberArray is intended as a Hex's replacement for std::vector\<NumberType\>.
 * Properties:
 * - correct memory alignment of the storage
 * - even number of allocated elements (zero padded) for the use in SSE accelerators
 * - basic iterator interface (members Array::begin(), Array::end()).
 * - HDF5 interface (ability to save and load to/from HDF5 data files)
 * - a collection of overloaded arithmetic operators (sum of two arrays,
 *   difference, multiplication by a number etc.)
 */
template <typename NumberType> class NumberArray : public Array<NumberType>
{
private:
    
    size_t N, Nres;
    NumberType * array;
    
    // allocate aligned memory
    NumberType* alloc_(size_t n)
    {
        // is there anything to allocate?
        if (n < 1)
            return nullptr;
        
        // allocate the aligned memory; make sure there will be even number of elements
        // so that we can always use pairs
        void* aligned_ptr = nullptr;
        int err = posix_memalign (
            &aligned_ptr,
            std::max(std::alignment_of<NumberType>::value, sizeof(void*)),
            (n + (n % 2)) * sizeof(NumberType)
        );
        
        // check memory allocation success
        if (err != 0)
            throw exception ("[NumberArray<T>::alloc_] Aligned memory allocation error (%d).", err);
        
        // get the number pointer
        NumberType* ptr = reinterpret_cast<NumberType*>(aligned_ptr);
        
        // clear the last element so that we may disregard it during multiplication
        *(ptr + n + (n % 2) - 1) = 0;
        
        // return the pointer
        return ptr;
    }
    
    // deallocate the memory
    void destroy_()
    {
        if (array != nullptr)
        {
            free(array);
            array = nullptr;
        }
    }
    
public:
    
    // alias
    typedef NumberType DataType;
    
    // inner product of two arrays
    template <typename NumberType1, typename NumberType2> friend auto operator | (
        NumberArray<NumberType1> const & a, NumberArray<NumberType2> const & b
    ) -> decltype(NumberType1(0)*NumberType2(0));
    
    // default constructor, creates an empty array
    NumberArray() : N(0), Nres(0), array(nullptr) {}
    
    // constructor, creates a length-n "x"-filled array
    NumberArray(size_t n, NumberType x = 0) : N(n), Nres(n)
    {
        // reserve space
        array = alloc_(Nres);
        
        // set to zero
        for (size_t i = 0; i < N; i++)
            array[i] = x;
    }
    
    // constructor, copies a length-n "array
    NumberArray(size_t n, NumberType* x) : N(n), Nres(n)
    {
        // reserve space
        array = alloc_(Nres);
        
        // set to zero
        for (size_t i = 0; i < N; i++)
            array[i] = x[i];
    }
    
    // copy constructor from Array const lvalue reference
    NumberArray(NumberArray<NumberType> const & a) : Array<NumberType>()
    {
        // reserve space
        N = Nres = a.N;
        array = alloc_(Nres);
        
        // run over the elements
        for (size_t i = 0; i < N; i++)
            array[i] = a.array[i];
    }
    
    // copy constructor from ArrayView const lvalue reference
    NumberArray(ArrayView<NumberType> const & a)
    {
        // reserve space
        N = Nres = a.size();
        array = alloc_(Nres);
        
        // run over the elements
        for (size_t i = 0; i < N; i++)
            array[i] = a.data()[i];
    }
    
    // copy constructor from Array rvalue reference
    NumberArray(NumberArray<NumberType> && a)
    {
        // copy content
        N = a.N;
        Nres = a.Nres;
        array = a.array;
        
        // clear rvalue
        a.N = 0;
        a.Nres = 0;
        a.array = nullptr;
    }
    
    // copy constructor from std::vector
    NumberArray(std::vector<NumberType> const & a)
    {
        // reserve space
        N = Nres = a.size();
        array = alloc_(Nres);
        
        // run over the elements
        for (size_t i = 0; i < N; i++)
            array[i] = a[i];
    }
    
    // copy constructor from initializer list
    NumberArray(std::initializer_list<NumberType> a)
    {
        // reserve space
        N = Nres = a.end() - a.begin();
        array = alloc_(Nres);
        
        // run over the elements
        size_t i = 0;
        for (auto it = a.begin(); it != a.end(); it++)
            array[i++] = *it;
    }
    
    // copy constructor from two forward iterators
    template <typename ForwardIterator> NumberArray(ForwardIterator i, ForwardIterator j)
    {
        // compute size
        N = 0;
        for (ForwardIterator k = i; k != j; k++)
            N++;
        
        // reserve space
        Nres = N;
        array = alloc_(Nres);
        
        // run over the elements
        size_t n = 0;
        for (ForwardIterator k = i; k != j; k++)
            array[n++] = *k;
    }
    
    // destructor
    ~NumberArray()
    {
        destroy_();
    }
    
    //
    // storage size
    //
    
    size_t size() const { return N; }
    size_t resize (size_t n)
    {
        if (n <= Nres)
        {
            N = n;
            return N;
        }
        
        Nres = n;
        NumberType * new_array = alloc_(Nres);
        
        for (size_t i = 0; i < n; i++)
            new_array[i] = (i < N) ? array[i] : NumberType(0);
        
        destroy_();
        N = n;
        
        array = new_array;
        
        return N;
    }
    size_t reserve (size_t n)
    {
        if (n > Nres)
        {
            Nres = n;
            NumberType* new_array = alloc_(Nres);
            
            if (N > 0)
            {
                memcpy(new_array, array, N * sizeof(NumberType));
                destroy_();
            }
            
            array = new_array;
        }
        
        return Nres;
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
    // data pointer to the aligned memory
    //
    
    NumberType* data()
    {
        return (NumberType *)aligned(array, std::max(alignof(NumberType),sizeof(Complex)));
    }
    NumberType const * data() const
    {
        return (NumberType * const)aligned(array, std::max(alignof(NumberType),sizeof(Complex)));
    }
    
    //
    // STL-like iterator interface
    //
    
    typedef NumberType* iterator;
    typedef const NumberType* const_iterator;
    iterator begin()
    { return array; }
    const_iterator begin() const
    { return array; }
    const_iterator cbegin() const
    { return array; }
    iterator end()
    { return array + N; }
    const_iterator end() const
    { return array + N; }
    const_iterator cend() const
    { return array + N; }
    NumberType & front(int i = 0)
    { return *(array + i); }
    NumberType const & front(int i = 0) const
    { return *(array + i); }
    NumberType & back(int i = 0)
    { return *(array + N - 1 - i); }
    NumberType const & back(int i = 0) const
    { return *(array + N - 1 - i); }
    
    void push_back(NumberType const & a)
    {
        if (N + 1 > Nres)
        {
            // double the capacity
            Nres = 2 * Nres + 1;
            
            // allocate space
            NumberType* new_array = alloc_(Nres);
            
            // copy original data
            memcpy(new_array, array, N * sizeof(NumberType));
            
            // destroy original array
            if (array != nullptr)
                destroy_();
            
            // use new array
            array = new_array;
        }
        
        // copy new element
        array[N++] = a;
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
        if (N + last - first > (int)Nres)
        {
            // raise the capacity
            Nres += last - first;
            
            // allocate space
            NumberType* new_array = alloc_(Nres);
            
            // copy original data
            memcpy(new_array, array, N * sizeof(NumberType));
            
            // destroy original array
            if (array != nullptr)
                destroy_();
            
            // use new array
            array = new_array;
        }
        
        for (InputIterator it = first; it != last; it++)
            array[N++] = *it;
    }
    
    bool empty() const
    {
        return N == 0;
    }
    
    void clear()
    {
        memset(array, 0, N * sizeof(NumberType));
    }
    
    //
    // assignment operators
    //
    
    NumberArray<NumberType>& operator = (NumberArray<NumberType> const &  b)
    {
        if (array == nullptr or Nres < b.N)
        {
            // delete insufficient storage
            destroy_();
            
            // allocate new storage
            Nres = b.N;
            array = alloc_(Nres);
        }
        
        // set the new dimension
        N = b.N;
        
        // run over the elements
        for (size_t i = 0; i < N; i++)
            array[i] = b.array[i];
        
        return *this;
    }
    
    NumberArray<NumberType>& operator = (NumberArray<NumberType> &&  b)
    {
        std::swap(N, b.N);
        std::swap(Nres, b.Nres);
        std::swap(array, b.array);
        return *this;
    }
    
    /**
     * Complex conjugate.
     */
    NumberArray<NumberType> conj() const
    {
        NumberArray<NumberType> c = *this;
        for (size_t i = 0; i < N; i++)
        {
            Complex z = c.array[i];
            c.array[i] = Complex(z.real(), -z.imag());
        }
        return c;
    }
    
    /**
     * Computes usual 2-norm.
     */
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
    
    /** 
     * Applies a user transformation.
     */
    template <class Functor> auto transform(Functor f) -> Array<decltype(f(NumberType(0)))>
    {
        Array<decltype(f(NumberType(0)))> c(N);
        for (size_t i = 0; i < N; i++)
            c[i] = f(array[i]);
        return c;
    }
    
    /**
     * Returns a subarray.
     */
    ArrayView<NumberType> slice(size_t left, size_t right) const
    {
        assert (right >= left);
        
        return ArrayView<NumberType>(array + left, array + right);
    }
    
    /**
     * Converts contents to SQL-readable BLOB (hexadecimal text format)
     * 
     * @warning The data are stored in the endianness of the current machine.
     */
    virtual std::string toBlob() const
    {
        // get byte pointer
        unsigned char const * dataptr = reinterpret_cast<unsigned char const*>(array);
        
        // get byte count
        size_t count = N * sizeof(NumberType);
        
        // resulting string
        std::ostringstream hexa;
        hexa << "x'" << std::hex << std::setfill('0');
        
        // for all bytes
        for (size_t i = 0; i < count; i++)
            hexa << std::setw(2) << static_cast<unsigned>(dataptr[i]);
        
        hexa << "'";
        return hexa.str();
    }
    
    /**
     * Decode string from SQL-readable BLOB (hexadecimal format) to correct binary array.
     * 
     * @warning The data are assumed to posess the endianness of the current machine.
     */
    virtual void fromBlob(std::string const & s)
    {
        if (array != nullptr and N != 0)
            destroy_();
        
        // the first character outght to be "x" or "X"
        // the second character outght to be "'"
        // the last character ought to be "'" as well
        if ((s[0] != 'x' and s[0] != 'X') or s[1] != '\'' or s.back() != '\'')
            throw exception ("[Array::fromBlob] Blob has wrong format, %s.", s.c_str());
        
        // create substring
        std::string ss(s.begin() + 2, s.end() - 1);
        
        // compute size
        size_t bytes = ss.size() / 2;
        N = bytes / sizeof(NumberType);
        
        // allocate space
        array = alloc_(N);
        
        // for all bytes
        for (size_t i = 0; i < bytes; i++)
        {
            unsigned byte;
            std::stringstream sst;
            
            // put two hexa digits from "ss" to "sst"
            sst << std::hex << ss.substr(2*i, 2);
            
            // read those two digits as a single byte
            sst >> byte;
            
            // store this byte
            reinterpret_cast<char*>(array)[i] = byte;
        }
    }
    
#ifndef NO_HDF
    /**
     * @brief Save array to HDF file.
     * 
     * This function will save the array data into a HDF5 file.
     * The file will then contain two datasets:
     * - "array" - the data itself, with large block of zeros omitted (if requested by
     *    the @c docompress flag
     * - "zero_blocks" - the compression information containing starting and ending
     *   positions of the zero blocks (hence the length of this dataset os always even)
     * 
     * @param name Filename.
     * @param docompress Whether to apply a trivial compression (contract the repeated zeros).
     * @param consec Minimal consecutive occurences for compression.
     */
    bool hdfsave(const char* name, bool docompress = false, int consec = 10) const
    {
        // save to HDF file
        HDFFile hdf(name, HDFFile::overwrite);
        
        if (not hdf.valid())
            return false;
        
        if (docompress)
        {
            NumberArray<int> zero_blocks;
            NumberArray<NumberType> elements;
            std::tie(zero_blocks,elements) = compress(consec);
            
            if (not zero_blocks.empty())
            {
                if (not hdf.write (
                    "zero_blocks",
                    &(zero_blocks[0]),
                    zero_blocks.size()
                )) return false;
            }
            if (not elements.empty())
            {
                if (not hdf.write (
                    "array",
                    &(elements[0]),
                    elements.size()
                )) return false;
            }
        }
        else
        {
            if (not hdf.write("array", array, N))
                return false;
        }
        
        return true;
    }
    
    /**
     * Load array from HDF file.
     * @param name Filename.
     */
    bool hdfload(std::string name)
    {
        // open the file
        HDFFile hdf(name, HDFFile::readonly);
        
        if (not hdf.valid())
            return false;
        
        NumberArray<NumberType> elements;
        NumberArray<int> zero_blocks;
        
        // read zero blocks
        if (zero_blocks.resize(hdf.size("zero_blocks")))
            if (not hdf.read("zero_blocks", &(zero_blocks[0]), zero_blocks.size()))
                return false;
            
            // get data size
            size_t size = hdf.size("array");
        
        // scale size if this is a complex array
        if (typeid(NumberType) == typeid(Complex))
            size /= 2;
        
        // read packed elements
        if (elements.resize(size))
        {
            if (not hdf.read("array", &(elements[0]), elements.size()))
                return false;
        }
        
        // remove previous data
        if (array != nullptr)
        {
            destroy_();
            array = nullptr;
            N = Nres = 0;
        }
        
        // unpack
        *this = elements.decompress(zero_blocks);
        
        return true;
    }
    
    // get compressed array
    std::tuple<NumberArray<int>,NumberArray<DataType>> compress(int consec) const
    {
        // compressed array
        NumberArray<NumberType> carray;
        carray.reserve(N);
        
        // zero blocks
        NumberArray<int> zero_blocks;
        
        // consecutive zeros counter
        int zero_counter = 0;
        
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
        NumberArray<int> nonzero_blocks;
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
        
        return std::make_tuple(zero_blocks,carray);
    }
    
    // get un-compressed array if compression info is supplied
    NumberArray<DataType> decompress(NumberArray<int> const & zero_blocks) const
    {
        if (zero_blocks.empty())
            return *this;
        
        // compute final size
        size_t final_size = N;
        for (size_t i = 0; i < zero_blocks.size()/2; i++)
            final_size += zero_blocks[2*i+1] - zero_blocks[2*i];
        
        // resize and clean internal storage
        NumberArray<DataType> unpack(final_size);
        memset(&(unpack[0]), 0, final_size * sizeof(NumberType));
        
        // copy nonzero chunks
        int this_end = 0;   // index of last updated element in "this"
        int load_end = 0;   // index of last used element in "nnz_array"
        for (size_t i = 0; i < zero_blocks.size()/2; i++)
        {
            int zero_start = zero_blocks[2*i];
            int zero_end = zero_blocks[2*i+1];
            
            // append nonzero data before this zero block
            memcpy (
                &(unpack[0]) + this_end,
                    array + load_end,
                    (zero_start - this_end) * sizeof(NumberType)
            );
            
            // move cursors
            load_end += zero_start - this_end;
            this_end  = zero_end;
        }
        
        // append remaining data
        memcpy (
            &(unpack[0]) + this_end,
                array + load_end,
                (final_size - this_end) * sizeof(NumberType)
        );
        
        return unpack;
    }
#endif
};

#include "arrithm.h"

// scalar product of two arrays.
template <typename NumberType1, typename NumberType2> auto operator | (
    NumberArray<NumberType1> const & a, NumberArray<NumberType2> const & b
) -> decltype(NumberType1(0)*NumberType2(0))
{
    // store size
    size_t N = a.N;
    
    // check if sizes match
    assert(N == b.N);
    
    // the scalar product
    decltype(NumberType1(0)*NumberType2(0)) result = 0;
    
    // iterators
    NumberType1 const * const restrict pa = &a[0];
    NumberType2 const * const restrict pb = &b[0];
    
    // sum the products
    for (size_t i = 0; i < N; i++)
        result += pa[i] * pb[i];
    
    return result;
}

template <class T1, class T2> auto outer_product (
    const ArrayView<T1> a,
    const ArrayView<T2> b
) -> NumberArray<decltype(T1(0) * T2(0))>
{
    NumberArray<decltype(T1(0) * T2(0))> c (a.size() * b.size());
    
    auto ic = c.begin();
    
    for (auto a_ : a)
    for (auto b_ : b)
        *(ic++) = a_ * b_;
    
    return c;
}

// output to text stream.
template <typename T> std::ostream & operator << (std::ostream & out, ArrayView<T> const & a)
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
 * @param start Left boundary and first sample for "samples" > 0.
 * @param end Right boundary and last sample for "samples" > 1.
 * @param samples Sample count.
 */
template <typename T> NumberArray<T> linspace(T start, T end, unsigned samples)
{
    NumberArray<T> space(samples);
    
    if (samples == 0)
        return space;
    
    if (samples == 1)
    {
        space[0] = start;
        return space;
    }
    
    for (unsigned i = 0; i < samples; i++)
        space[i] = start + ((end - start) * T(i)) / T(samples - 1);
    
    return space;
}

/**
 * Generate logarithmic grid
 * @param x0 Left boundary and first sample for "samples" > 0.
 * @param x1 Right boundary and last sample for "samples" > 1.
 * @param N Sample count.
 */
template <typename T> NumberArray<T> logspace(T x0, T x1, size_t N)
{
    if (x0 <= 0 or x1 <= 0 or x1 < x0)
    {
        std::cerr << "[logspace] It must be 0 < x1 <= x2 !\n";
        abort();
    }
    
    NumberArray<T> grid(N);
    
    if (N == 1)
        grid[0] = x0;
    
    if (N > 1)
        for (unsigned i = 0; i < N; i++)
            grid[i] = x0 * pow(x1 / x0, i / T(N - 1));
        
        return grid;
}

/**
 * Write array to file. Array will be written as a single column into
 * an ASCII file.
 * @param array The array to write.
 * @param filename Name of the file to create/overwrite.
 */
template <typename NumberType> void write_array (
    ArrayView<NumberType> array,
    const char* filename
);

/**
 * Write array to file. Array will be written as a two columns into
 * an ASCII file, first column contains grid information.
 * @param grid The 1D grid labels for the data.
 * @param array The array to write.
 * @param filename Name of the file to create/overwrite.
 */
template <typename NumberType> void write_array (
    ArrayView<double> grid,
    ArrayView<NumberType> array,
    const char* filename
);


template <typename Fetcher> bool write_1D_data (size_t m, const char* filename, Fetcher fetch)
{
    std::ofstream f(filename);
    
    if (f.bad())
        return false;
    
    for (size_t i = 0; i < m; i++)
        f << fetch(i) << "\n";
    
    return true;
}

/**
 * Write 2D data to file. To allow maximum flexibility, only extensions
 * of the data are passed to the function and a functor that will be
 * repeatedly called with coordinate pair for new data element.
 * @param m Row count.
 * @param n Column count.
 * @param filename Filename of the file to create/overwrite.
 * @param fetch Functor with interface
 *        @code
 *             double operator() (size_t, size_t);
 *        @endcode
 * @return Write success indicator (@c false for failure).
 */
template <class Fetcher> bool write_2D_data(size_t m, size_t n, const char* filename, Fetcher fetch)
{
    std::ofstream f(filename);
    
    if (f.bad())
        return false;
    
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
            f << fetch(i,j) << " ";
        f << "\n";
    }
    
    return true;
}

//
// aliases
//

typedef NumberArray<int>          iArray;
typedef NumberArray<long>         lArray;
typedef NumberArray<double>       rArray;
typedef NumberArray<long double>  qArray;
typedef NumberArray<Complex>      cArray;

typedef ArrayView<int>         iArrayView;
typedef ArrayView<long>        lArrayView;
typedef ArrayView<double>      rArrayView;
typedef ArrayView<long double> qArrayView;
typedef ArrayView<Complex>     cArrayView;

typedef Array<iArray>             iArrays;
typedef Array<lArray>             lArrays;
typedef Array<rArray>             rArrays;
typedef Array<qArray>             qArrays;
typedef Array<cArray>             cArrays;

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
 * @param v1 First array.
 * @param ...p All other arrays.
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

// return absolute values
rArray abs (cArray const &u);
rArrays abs (cArrays const &u);

// boolean aggregation
bool all(Array<bool> v);
bool any(Array<bool> v);

// minimal element
template <typename NumberType> NumberType min (ArrayView<NumberType> const & a)
{
    NumberType z = a.front();
    for (NumberType const * it = a.begin(); it != a.end(); it++)
        if (*it < z)
            z = *it;
        return z;
}

// maximal element
template <typename NumberType> NumberType max (ArrayView<NumberType> const & a)
{
    NumberType z = a.front();
    for (NumberType const * it = a.begin(); it != a.end(); it++)
        if (*it > z)
            z = *it;
        return z;
}

// return per-element power
template <typename T> Array<T> pow(ArrayView<T> const & u, double e)
{
    Array<T> v(u.size());
    
    auto iu = u.begin();
    auto iv = v.begin();
    
    while (iu != u.end())
        *(iv++) = pow(*(iu++), e);
    
    return v;
}

template <typename NumberType> NumberArray<NumberType> sqrt (ArrayView<NumberType> const & A)
{
    size_t N = A.size();
    NumberArray<NumberType> B (N);
    
    for (size_t i = 0; i < N; i++)
        B[i] = sqrt(A[i]);
    
    return B;
}

NumberArray<double> hypot (NumberArray<double> const & A, NumberArray<double> const & B);
NumberArray<double> atan2 (NumberArray<double> const & A, NumberArray<double> const & B);
NumberArray<double> sqrabs (NumberArray<Complex> const & A);
NumberArray<double> realpart (NumberArray<Complex> const & A);
NumberArray<double> imagpart (NumberArray<Complex> const & A);

// summation
template <typename T> T sum(Array<T> v)
{
    return std::accumulate(v.begin(), v.end(), T(0));
}

// summation of nested arrays
template <typename T> NumberArray<T> sums(Array<NumberArray<T>> v)
{
    if (v.size() == 0)
        return NumberArray<T>();	// empty array
        
        return std::accumulate (
            v.begin(),
            v.end(),
            NumberArray<T> (v[0].size()),
            [](NumberArray<T> a, NumberArray<T> b) -> NumberArray<T> {
                return a + b;
            }
        );
}

/**
 * Comparison of an array and a number.
 * @param u Array.
 * @param x Number.
 * @return Vector of bools for element-wise comparisons.
 */
template <typename T>
Array<bool> operator == (ArrayView<T> u, T x)
{
    Array<bool> v(u.size());
    for (size_t i  = 0; i < u.size(); i++)
        v[i] = (u[i] == x);
    return v;
}

template <typename T>
Array<bool> operator == (ArrayView<T> u, ArrayView<T> v)
{
    assert(u.size() == v.size());
    
    Array<bool> w(u.size());
    for (size_t i  = 0; i < u.size(); i++)
        w[i] = (u[i] == v[i]);
    return w;
}

inline bool all(Array<bool> B)
{
    bool ok = true;
    for (bool b : B)
        ok = ok and b;
    return ok;
}

inline bool any(Array<bool> B)
{
    bool ok = false;
    for (bool b : B)
        ok = ok or b;
    return ok;
}

/**
 * Evaluation of a function over a grid.
 * @param f Functor to be evaluated using the operator() (double) interface.
 * @param grid Array-like type containing the evaluation points.
 * @param vals Array-like type to hold evaluated function on return. It is
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
 * @param idx1 Sorted (!) indices of the first array.
 * @param arr1 Merge TO array.
 * @param idx2 Sorted (!) indices of the second array.
 * @param arr2 Merge FROM array.
 */
template <typename Tidx, typename Tval> void merge (
    NumberArray<Tidx>       & idx1, NumberArray<Tval>       & arr1,
    NumberArray<Tidx> const & idx2, NumberArray<Tval> const & arr2
){
    // positions in arrays
    size_t i1 = 0;
    size_t i2 = 0;
    
    // output arrays
    NumberArray<Tidx> idx;
    NumberArray<Tval> arr;
    
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

/**
 * Join elements from all subarrays.
 */
template <typename T> NumberArray<T> join (Array<NumberArray<T>> const & arrays)
{
    NumberArray<size_t> partial_sizes(arrays.size() + 1);
    
    // get partial sizes
    partial_sizes[0] = 0;
    for (size_t i = 0; i < arrays.size(); i++)
        partial_sizes[i+1] = partial_sizes[i] + arrays[i].size();
    
    // result array
    NumberArray<T> res(partial_sizes.back());
    
    // concatenate arrays
    for (size_t i = 0; i < arrays.size(); i++)
    {
        if (arrays[i].size() > 0)
            ArrayView<T>(res, partial_sizes[i], arrays[i].size()) = arrays[i];
    }
    
    return res;
}

/**
 * Drop all redundant repetitions from sorted array.
 */
template <class T> NumberArray<T> sorted_unique (const ArrayView<T> v, int n = 1)
{
    // create output array
    NumberArray<T> w(v.size());
    
    // iterators
    T * iw = w.begin();
    T const * iv = v.begin();
    int repeat = 0;
    
    // copy all elements, drop repetitions above "repeat"
    while (iv != v.end())
    {
        // use this element definitely if
        // - it is the first element
        // - it differs from the preceding element
        if (iw == w.begin() or *(iw - 1) != *iv)
            repeat = 0;
        
        // use the conforming element
        if (++repeat < n)
            *(iw++) = *iv;
        
        // move on to the next element
        iv++;
    }
    
    // resize and return output array
    w.resize(iw - w.begin());
    return w;
}

/**
 * @brief Running average.
 * 
 * Overwrite the given vector v = {v1, v2, v3, ..., vN} by vector
 * {v1, (v1+v2)/2, (v2+v3)/2, ..., (v_(N-1)+vN)/2}.
 */
template <class T> void smoothen (ArrayView<T> v)
{
    if (v.size() == 0)
        return;
    
    T prev = v[0], newv;
    for (int i = 1; i < v.size(); i++)
    {
        newv = 0.5 * (prev + v[i]);
        prev = v[i];
        v[i] = newv;
    }
}

/**
 * @brief Convert ArrayView<T> to a string.
 */
template <class T> std::string to_string (const ArrayView<T> v)
{
    std::ostringstream ss;
    for (T const & x : v)
        ss << x << " ";
    return ss.str();
}

#endif
