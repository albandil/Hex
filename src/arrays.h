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

template <class T> class PlainAllocator
{
    public:
        static T * alloc (size_t n) { return new T[n](); }
        static void free (T * ptr) { if (ptr != nullptr) delete [] ptr; }
};

template <class T, size_t alignment = std::alignment_of<T>::value> class AlignedAllocator
{
    public:
        // allocate aligned memory
        static T * alloc (size_t n)
        {
            // is there anything to allocate?
            if (n < 1)
                return nullptr;
            
            // allocate the aligned memory; make sure there will be even number of elements
            // so that we can always use pairs
            void* aligned_ptr = nullptr;
            int err = posix_memalign (
                &aligned_ptr,
                std::max(alignment, sizeof(void*)),
                (n + (n % 2)) * sizeof(T)
            );
            
            // check memory allocation success
            if (err != 0)
                throw exception ("[AlignedAllocator<T>::alloc] Aligned memory allocation error (%d).", err);
            
            // get the number pointer
            T* ptr = reinterpret_cast<T*>(aligned_ptr);
            
            // clear the last element so that we may disregard it during multiplication
            *(ptr + n + (n % 2) - 1) = 0;
            
            // return the pointer
            return ptr;
        }
        
        // deallocate the memory
        static void free (T * ptr)
        {
            if (ptr != nullptr)
                ::free (ptr);
        }
};

template <class T, class Alloc = PlainAllocator<T>> class Array;
template <class T, class Alloc = AlignedAllocator<T>> class NumberArray;

/**
 * @brief Array shallow copy.
 */
template <class T> class ArrayView
{
    public:
        
        // aliases
        typedef T DataType;
        typedef T * iterator;
        typedef T const * const_iterator;
        
    protected:
        
        /// Number of elements in the array.
        size_t N_;
        
        /// Pointer to the array.
        T * array_;
        
    public:
        
        // empty constructor
        ArrayView ()
            : N_(0), array_(nullptr) {}
        
        // construct from pointer and size
        ArrayView (size_t n, T * ptr)
            : N_(n), array_(ptr) {}
        
        // copy constructor
        ArrayView (ArrayView<T> const & a, size_t i = 0, size_t n = 0)
            : N_((n > 0) ? n : a.size()), array_(const_cast<T*>(a.data()) + i) { assert(i + N_ <= a.size()); }
        
        // construct view from Array const lvalue reference
        ArrayView (Array<T> const & a, size_t i = 0, size_t n = 0)
            : N_((n > 0) ? n : a.size()), array_(const_cast<T*>(a.data()) + i) { assert(i + N_ <= a.size()); }
            
        // construct view from NumberArray const lvalue reference
        ArrayView (NumberArray<T> const & a, size_t i = 0, size_t n = 0)
            : N_((n > 0) ? n : a.size()), array_(const_cast<T*>(a.data()) + i) { assert(i + N_ <= a.size()); }
        
        // construct from consecutive memory segment
        ArrayView (const_iterator i, const_iterator j)
            : N_(j - i), array_(const_cast<T*>(&(*i))) {}
        
        // construct from right-value reference
        ArrayView (ArrayView<T> && r)
            : ArrayView() { std::swap(N_, r.N_); std::swap(array_, r.array_); }
        
        // destructor
        virtual ~ArrayView () {}
    
        // assignments
        ArrayView<T> & operator = (const ArrayView<T> v)
        {
            if (v.size() != size())
                throw exception ("[ArrayView::operator=] Cannot copy %ld elements to %ld fields!", v.size(), N_);
            
            for (size_t i = 0; i < size(); i++)
                array_[i] = v[i];
            
            return *this;
        }
    
        // element-wise access (non-const)
        T & operator[] (size_t i)
        {
#ifdef NDEBUG
            return array_[i];
#else
            // bounds check
            if (i < size())
                return array_[i];
            else
                throw exception ("[ArrayView::operator[]] Index %ld out of bounds (size = %ld) !", i, N_);
#endif
        }
    
        // element-wise access (const)
        T const & operator[] (size_t i) const
        {
#ifdef NDEBUG
            return array_[i];
#else
            if (i < size())
                return array_[i];
            else
                throw exception ("[ArrayView::operator[]] Index %ld out of bounds (size = %ld) !", i, N_);
#endif
        }
        
        // getters
        size_t size () const { return N_; }
        
        // data pointer
        virtual T * data () { return array_; }
        virtual T const * data () const { return array_; }
    
        //
        // STL-like iterator interface
        //
        
        iterator begin ()
            { return data(); }
        const_iterator begin () const
            { return data(); }
        iterator end ()
            { return data() + size(); }
        const_iterator end () const
            { return data() + size(); }
        T & front (int i = 0)
            { return *(data() + i); }
        T const & front (int i = 0) const
            { return *(data() + i); }
        T & back (int i = 0)
            { return *(data() + size() - 1 - i); }
        T const & back (int i = 0) const
            { return *(data() + size() - 1 - i); }
        
        // some other functions
        void fill (T x) { for (T & y : *this) y = x; }
        bool empty () const { return size() == 0; }
        
        // 2-norm, defined only for scalar NumberType
        template <class = typename std::enable_if<is_scalar<T>::value>> double norm () const
        {
            double sqrnorm = 0.;
            for (T const & x : *this)
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
template <class T, class Alloc> class Array : public ArrayView<T>
{
    public:
        
        // aliases
        typedef T DataType;
        typedef T * iterator;
        typedef T const * const_iterator;
    
    public:
    
        // default constructor, creates an empty array
        Array ()
            : ArrayView<T>() {}
        
        // constructor, creates a length-n "x"-filled array
        Array (size_t n, T x = 0)
            : ArrayView<T>(n, Alloc::alloc(n)) { for (size_t i = 0; i < size(); i++) (*this)[i] = x; }
        
        // constructor, copies a length-n "array
        Array (size_t n, T const * const x)
            : ArrayView<T>(n, Alloc::alloc(n)) { for (size_t i = 0; i < size(); i++) (*this)[i] = x[i]; }
        
        // copy constructor from Array const lvalue reference
        Array (Array<T> const & a)
            : ArrayView<T>(a.size(), Alloc::alloc(a.size())) { for (size_t i = 0; i < size(); i++) (*this)[i] = a[i]; }
        
        // copy constructor from const ArrayView
        Array (const ArrayView<T> a)
            : ArrayView<T>(a.size(), Alloc::alloc(a.size())) { for (size_t i = 0; i < size(); i++) (*this)[i] = a[i]; }
        
        // copy constructor from Array rvalue reference
        Array (Array<T> && a)
            : ArrayView<T>() { std::swap (ArrayView<T>::N_, a.ArrayView<T>::N_); std::swap (ArrayView<T>::array_, a.ArrayView<T>::array_); }
        
        // copy constructor from std::vector
        Array (std::vector<T> const & a)
            : ArrayView<T>(a.size(), Alloc::alloc(a.size())) { for (size_t i = 0; i < size(); i++) (*this)[i] = a[i]; }
        
        // copy constructor from initializer list
        Array (std::initializer_list<T> a)
            : ArrayView<T>(a.end()-a.begin(), Alloc::alloc(a.end()-a.begin())) { size_t i = 0; for (auto it = a.begin(); it != a.end(); it++) (*this)[i++] = *it; }
        
        // copy constructor from two forward iterators
        template <typename ForwardIterator> Array (ForwardIterator i, ForwardIterator j)
            : ArrayView<T>(std::distance(i,j), Alloc::alloc(std::distance(i,j)))
        {
            size_t n = 0;
            for (ForwardIterator k = i; k != j; k++)
                (*this)[n++] = *k;
        }
        
        // destructor
        ~Array ()
        {
            if (ArrayView<T>::array_ != nullptr)
            {
                Alloc::free (ArrayView<T>::array_);
                ArrayView<T>::array_ = nullptr;
            }
        }
        
        // storage size
        size_t size () const { return ArrayView<T>::size(); }
        virtual size_t resize (size_t n)
        {
            if (n == 0)
            {
                if (ArrayView<T>::array_ != nullptr)
                {
                    Alloc::free (ArrayView<T>::array_);
                    ArrayView<T>::array_ = nullptr;
                    ArrayView<T>::N_ = 0;
                }
                return 0;
            }
            
            T * new_array = Alloc::alloc(n);
            for (size_t i = 0; i < n; i++)
                new_array[i] = (i < size()) ? (*this)[i] : T(0);
            
            Alloc::free (ArrayView<T>::array_);
            
            ArrayView<T>::N_ = n;
            ArrayView<T>::array_ = new_array;
            
            return size();
        }
        
        // element-wise access (non-const)
        inline T & operator[] (size_t i) { return ArrayView<T>::operator[](i); }

        // element-wise access (const)
        inline T const & operator[] (size_t i) const { return ArrayView<T>::operator[](i); }
    
        // data pointer
        virtual T* data () { return ArrayView<T>::data(); }
        virtual T const * data () const { return ArrayView<T>::data(); }
    
        //
        // STL-like iterator interface
        //
        
        iterator begin ()                  { return ArrayView<T>::begin(); }
        const_iterator begin () const      { return ArrayView<T>::begin(); }
        iterator end ()                    { return ArrayView<T>::end(); }
        const_iterator end () const        { return ArrayView<T>::end(); }
        T & front (int i = 0)              { return ArrayView<T>::front(i); }
        T const & front (int i = 0) const  { return ArrayView<T>::front(i); }
        T & back (int i = 0)               { return ArrayView<T>::back(i); }
        T const & back (int i = 0) const   { return ArrayView<T>::back(i); }
        
        virtual void push_back (T const & a)
        {
            // not very efficient... FIXME
            
            T* new_array = Alloc::alloc (size() + 1);
            for (size_t i = 0; i < size(); i++)
                new_array[i] = std::move(ArrayView<T>::array_[i]);
            new_array[size()] = a;
            ArrayView<T>::N_++;
            Alloc::free (ArrayView<T>::array_);
            ArrayView<T>::array_ = new_array;
        }
        
        virtual T pop_back ()
        {
            if (size() > 0)
                return *(data() + (--ArrayView<T>::N_));
            else
                throw exception ("Array has no element to pop!");
        }
        
        template <class InputIterator> void append (
            InputIterator first, InputIterator last
        ) {
            T* new_array = Alloc::alloc(size() + last - first);
            for (size_t i = 0; i < size(); i++)
                new_array[i] = (*this)[i];
            for (InputIterator it = first; it != last; it++)
                new_array[size() + it - first] = *it;
            ArrayView<T>::N_ += last - first;
            Alloc::free (ArrayView<T>::array_);
            ArrayView<T>::array_ = new_array;
        }
        
        void insert (iterator it, T x)
        {
            // create new array (one element longer)
            T* new_array = Alloc::alloc(size() + 1);
            
            // copy everything to the new location
            for (int i = 0; i < it - begin(); i++)
                new_array[i] = std::move((*this)[i]);
            
            // insert new element
            *(new_array + (it - begin())) = std::move(x);
            
            // copy the rest
            for (int i = it - begin(); i < (int)size(); i++)
                new_array[i+1] = std::move((*this)[i]);
            
            // change pointers
            Alloc::free (ArrayView<T>::array_);
            ArrayView<T>::N_++;
            ArrayView<T>::array_ = new_array;
        }
        
        bool empty () const { return size() == 0; } 
    
        //
        // assignment operators
        //
        
        Array<T>& operator = (Array<T> const &  b)
        {
            // if we already have some allocated space, check its size,
            // so that we do not free it uselessly
            if (data() != nullptr and size() != b.size())
            {
                Alloc::free (ArrayView<T>::array_);
                ArrayView<T>::array_ = nullptr;
            }
            
            // set the new dimension
            ArrayView<T>::N_ = b.size();
            
            // if necessary, reserve space
            if (data() == nullptr)
                ArrayView<T>::array_ = Alloc::alloc(size());
            
            // run over the elements
            for (size_t i = 0; i < size(); i++)
                (*this)[i] = b[i];
            
            return *this;
        }
        
        Array<T>& operator = (Array<T> &&  b)
        {
            // if we already have some allocated space, check its size,
            // so that we do not free it uselessly
            if (data() != nullptr)
            {
                Alloc::free (ArrayView<T>::array_);
                ArrayView<T>::array_ = nullptr;
                ArrayView<T>::N_ = 0;
            }
            
            // swap content
            std::swap (ArrayView<T>::N_, b.ArrayView<T>::N_);
            std::swap (ArrayView<T>::array_, b.ArrayView<T>::array_);
            
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
template <class T, class Alloc> class NumberArray : public Array<T, Alloc>
{
    public:
        
        // aliases
        typedef T DataType;
        typedef T * iterator;
        typedef T const * const_iterator;
        
    protected:
        
        /// Allocated memory.
        size_t Nres_;
        
    public:
        
        // default constructor, creates an empty array
        NumberArray ()
            : Array<T,Alloc>(), Nres_(size()) {}
        
        // constructor, creates a length-n "x"-filled array
        NumberArray (size_t n, T x = 0)
            : Array<T,Alloc>(n, x), Nres_(size()) {}
        
        // constructor, copies a length-n "array
        NumberArray (size_t n, T const * x)
            : Array<T,Alloc>(n, x), Nres_(size()) {}
        
        // copy constructor from ArrayView const lvalue reference
        NumberArray (ArrayView<T> const & a)
            : Array<T,Alloc>(a), Nres_(size()) {}
        
        // copy constructor from Array const lvalue reference
        NumberArray (NumberArray<T> const & a)
            : Array<T,Alloc>((ArrayView<T> const &)a), Nres_(size()) {}
        
        // copy constructor from std::vector
        NumberArray (std::vector<T> const & a)
            : Array<T,Alloc>(a), Nres_(size()) {}
        
        // copy constructor from initializer list
        NumberArray (std::initializer_list<T> a)
            : Array<T,Alloc>(a), Nres_(size()) {}
        
        // copy constructor from two forward iterators
        template <typename ForwardIterator> NumberArray (ForwardIterator i, ForwardIterator j)
            : Array<T,Alloc>(i,j), Nres_(size()) {}
        
        // copy constructor from Array rvalue reference
        NumberArray (NumberArray<T> && a)
            : Array<T,Alloc>()
        {
            std::swap(ArrayView<T>::N_, a.ArrayView<T>::N_);
            std::swap(ArrayView<T>::array_, a.ArrayView<T>::array_);
            std::swap(Nres_,a.Nres_);
        }
        
        // destructor
        ~NumberArray ()
        {
            // Array<T> will do the destruction
            // ... do nothing here
        }
        
        //
        // storage size
        //
        
        size_t size () const { return ArrayView<T>::size(); }
        virtual size_t resize (size_t n)
        {
            if (n <= Nres_)
            {
                ArrayView<T>::N_ = n;
                return n;
            }
            
            Nres_ = n;
            T * new_array = Alloc::alloc(Nres_);
            
            for (size_t i = 0; i < n; i++)
                new_array[i] = (i < size()) ? (*this)[i] : T(0);
            
            Alloc::free (ArrayView<T>::array_);
            ArrayView<T>::N_ = n;
            ArrayView<T>::array_ = new_array;
            
            return size();
        }
        size_t reserve (size_t n)
        {
            if (n > Nres_)
            {
                Nres_ = n;
                T* new_array = Alloc::alloc(Nres_);
                
                if (size() > 0)
                {
                    memcpy (new_array, data(), size() * sizeof(T));
                    Alloc::free (ArrayView<T>::array_);
                }
                
                ArrayView<T>::array_ = new_array;
            }
            
            return Nres_;
        }
        
        // element-wise access (non-const)
        inline T & operator[] (size_t i) { return ArrayView<T>::operator[](i); }
        
        // element-wise access (const)
        inline T const & operator[] (size_t i) const { return ArrayView<T>::operator[](i); }
        
        //
        // data pointer to the aligned memory
        //
        
        virtual T * data ()
        {
            return (T *)aligned_ptr(ArrayView<T>::array_, std::max(alignof(T),sizeof(Complex)));
        }
        virtual T const * data () const
        {
            return (T * const)aligned_ptr(ArrayView<T>::array_, std::max(alignof(T),sizeof(Complex)));
        }
        
        //
        // STL-like iterator interface
        //
        
        iterator begin ()                  { return ArrayView<T>::begin(); }
        const_iterator begin () const      { return ArrayView<T>::begin(); }
        const_iterator cbegin () const     { return ArrayView<T>::begin(); }
        iterator end ()                    { return ArrayView<T>::end(); }
        const_iterator end () const        { return ArrayView<T>::end(); }
        const_iterator cend () const       { return ArrayView<T>::end(); }
        T & front (int i = 0)              { return ArrayView<T>::front(i); }
        T const & front (int i = 0) const  { return ArrayView<T>::front(i); }
        T & back (int i = 0)               { return ArrayView<T>::back(i); }
        T const & back (int i = 0) const   { return ArrayView<T>::back(i); }
        
        virtual void push_back (T const & a)
        {
            if (size() + 1 > Nres_)
            {
                // double the capacity
                Nres_ = 2 * Nres_ + 1;
                
                // allocate space
                T* new_array = Alloc::alloc(Nres_);
                
                // copy original data
                memcpy(new_array, data(), size() * sizeof(T));
                
                // destroy original array
                if (data() != nullptr)
                    Alloc::free (ArrayView<T>::array_);
                
                // use new array
                ArrayView<T>::array_ = new_array;
            }
            
            // copy new element
            (*this)[ArrayView<T>::N_++] = a;
        }
        
        T pop_back()
        {
            return Array<T,Alloc>::pop_back();
        }
        
        template <class InputIterator> void append (
            InputIterator first, InputIterator last
        ) {
            if (size() + last - first > (int)Nres_)
            {
                // raise the capacity
                Nres_ += last - first;
                
                // allocate space
                T* new_array = Alloc::alloc(Nres_);
                
                // copy original data
                memcpy(new_array, data(), size() * sizeof(T));
                
                // destroy original array
                if (data() != nullptr)
                    Alloc::free (ArrayView<T>::array_);
                
                // use new array
                ArrayView<T>::array_ = new_array;
            }
            
            for (InputIterator it = first; it != last; it++)
                (*this)[ArrayView<T>::N_++] = *it;
        }
        
        bool empty () const
        {
            return size() == 0;
        }
        
        void clear ()
        {
            memset(ArrayView<T>::array_, 0, size() * sizeof(T));
        }
        
        //
        // assignment operators
        //
        
        NumberArray<T>& operator = (NumberArray<T> const &  b)
        {
            if (data() == nullptr or Nres_ < b.size())
            {
                // delete insufficient storage
                Alloc::free (ArrayView<T>::array_);
                
                // allocate new storage
                Nres_ = b.size();
                ArrayView<T>::array_ = Alloc::alloc(Nres_);
            }
            
            // set the new dimension
            ArrayView<T>::N_ = b.size();
            
            // run over the elements
            for (size_t i = 0; i < size(); i++)
                (*this)[i] = b[i];
            
            return *this;
        }
        
        NumberArray<T>& operator = (NumberArray<T> &&  b)
        {
            std::swap(ArrayView<T>::N_, b.ArrayView<T>::N_);
            std::swap(ArrayView<T>::array_, b.ArrayView<T>::array_);
            std::swap(Nres_, b.Nres_);
            return *this;
        }
        
        /**
        * Complex conjugate.
        */
        NumberArray<T> conj () const
        {
            NumberArray<T> c = *this;
            for (size_t i = 0; i < size(); i++)
            {
                Complex z = c[i];
                c[i] = Complex(z.real(), -z.imag());
            }
            return c;
        }
        
        /**
        * Computes usual 2-norm.
        */
        double norm () const
        {
            double ret = 0.;
            for (size_t i = 0; i < size(); i++)
            {
                Complex z = (*this)[i];
                ret += z.real() * z.real() + z.imag() * z.imag();
            }
            return sqrt(ret);
        }
        
        /** 
        * Applies a user transformation.
        */
        template <class Functor> auto transform (Functor f) -> NumberArray<decltype(f(T(0)))>
        {
            NumberArray<decltype(f(T(0)))> c(size());
            for (size_t i = 0; i < size(); i++)
                c[i] = f((*this)[i]);
            return c;
        }
        
        /**
        * Returns a subarray.
        */
        ArrayView<T> slice (size_t left, size_t right) const
        {
            assert (right >= left);
            
            return ArrayView<T>(begin() + left, begin() + right);
        }
        
        /**
        * Converts contents to SQL-readable BLOB (hexadecimal text format)
        * 
        * @warning The data are stored in the endianness of the current machine.
        */
        std::string toBlob () const
        {
            // get byte pointer
            unsigned char const * dataptr = reinterpret_cast<unsigned char const*>(data());
            
            // get byte count
            size_t count = size() * sizeof(T);
            
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
        void fromBlob (std::string const & s)
        {
            if (data() != nullptr and size() != 0)
                Alloc::free (ArrayView<T>::array_);
            
            // the first character outght to be "x" or "X"
            // the second character outght to be "'"
            // the last character ought to be "'" as well
            if ((s[0] != 'x' and s[0] != 'X') or s[1] != '\'' or s.back() != '\'')
                throw exception ("[NumberArray::fromBlob] Blob has wrong format, %s.", s.c_str());
            
            // create substring
            std::string ss (s.begin() + 2, s.end() - 1);
            
            // compute size
            size_t bytes = ss.size() / 2;
            ArrayView<T>::N_ = bytes / sizeof(T);
            
            // allocate space
            ArrayView<T>::array_ = Alloc::alloc(size());
            
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
                reinterpret_cast<char*>(data())[i] = byte;
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
        bool hdfsave (const char* name, bool docompress = false, int consec = 10) const
        {
            // save to HDF file
            HDFFile hdf(name, HDFFile::overwrite);
            
            if (not hdf.valid())
                return false;
            
            if (docompress)
            {
                NumberArray<int> zero_blocks;
                NumberArray<T> elements;
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
                if (not hdf.write("array", data(), size()))
                    return false;
            }
            
            return true;
        }
        
        /**
        * Load array from HDF file.
        * @param name Filename.
        */
        bool hdfload (std::string name)
        {
            // open the file
            HDFFile hdf(name, HDFFile::readonly);
            
            if (not hdf.valid())
                return false;
            
            NumberArray<T> elements;
            NumberArray<int> zero_blocks;
            
            // read zero blocks
            if (zero_blocks.resize(hdf.size("zero_blocks")))
                if (not hdf.read("zero_blocks", &(zero_blocks[0]), zero_blocks.size()))
                    return false;
                
                // get data size
                size_t size = hdf.size("array");
            
            // scale size if this is a complex array
            if (typeid(T) == typeid(Complex))
                size /= 2;
            
            // read packed elements
            if (elements.resize(size))
            {
                if (not hdf.read("array", &(elements[0]), elements.size()))
                    return false;
            }
            
            // remove previous data
            if (data() != nullptr)
            {
                Alloc::free (ArrayView<T>::array_);
                ArrayView<T>::array_ = nullptr;
                ArrayView<T>::N_ = Nres_ = 0;
            }
            
            // unpack
            *this = elements.decompress(zero_blocks);
            
            return true;
        }
        
        // get compressed array
        std::tuple<NumberArray<int>,NumberArray<T>> compress (int consec) const
        {
            // compressed array
            NumberArray<T> carray;
            carray.reserve(size());
            
            // zero blocks
            NumberArray<int> zero_blocks;
            
            // consecutive zeros counter
            int zero_counter = 0;
            
            // analyze: find compressible segments
            for (size_t i = 0; i < size(); i++)
            {
                if ((*this)[i] == 0.)
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
                zero_blocks.push_back(size()-zero_counter);
                zero_blocks.push_back(size());
            }
            
            // invert selection: get non-zero blocks
            NumberArray<int> nonzero_blocks;
            nonzero_blocks.push_back(0);
            nonzero_blocks.append(zero_blocks.begin(), zero_blocks.end());
            nonzero_blocks.push_back(size());
            
            // compress: copy only nonzero elements
            for (size_t iblock = 0; iblock < nonzero_blocks.size()/2; iblock++)
            {
                int start = nonzero_blocks[2*iblock];
                int end = nonzero_blocks[2*iblock+1];
                carray.append (
                    begin() + start,
                    begin() + end
                );
            }
            
            return std::make_tuple(zero_blocks,carray);
        }
        
        // get un-compressed array if compression info is supplied
        NumberArray<T> decompress (NumberArray<int> const & zero_blocks) const
        {
            if (zero_blocks.empty())
                return *this;
            
            // compute final size
            size_t final_size = size();
            for (size_t i = 0; i < zero_blocks.size()/2; i++)
                final_size += zero_blocks[2*i+1] - zero_blocks[2*i];
            
            // resize and clean internal storage
            NumberArray<DataType> unpack(final_size);
            memset(&(unpack[0]), 0, final_size * sizeof(T));
            
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
                    begin() + load_end,
                    (zero_start - this_end) * sizeof(T)
                );
                
                // move cursors
                load_end += zero_start - this_end;
                this_end  = zero_end;
            }
            
            // append remaining data
            memcpy (
                &(unpack[0]) + this_end,
                begin() + load_end,
                (final_size - this_end) * sizeof(T)
            );
            
            return unpack;
        }
#endif
};

#include "arrithm.h"

// scalar product of two arrays.
template <class T> T operator | (const ArrayView<T> a, const ArrayView<T> b)
{
    // get size; check if sizes match
    size_t N = a.size();
    assert(N == b.size());
    
    // the scalar product
    T result = 0;
    
    // iterators
    T const * const restrict pa = a.data();
    T const * const restrict pb = b.data();
    
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
template <typename T> NumberArray<T> linspace (T start, T end, unsigned samples)
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
template <typename T> NumberArray<T> logspace (T x0, T x1, size_t N)
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
    NumberArray<NumberType> array,
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
    NumberArray<double> grid,
    NumberArray<NumberType> array,
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
template <class Fetcher> bool write_2D_data (size_t m, size_t n, const char* filename, Fetcher fetch)
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
typedef NumberArray<Complex>      cArray;

typedef ArrayView<int>         iArrayView;
typedef ArrayView<long>        lArrayView;
typedef ArrayView<double>      rArrayView;
typedef ArrayView<Complex>     cArrayView;

typedef Array<iArray>             iArrays;
typedef Array<lArray>             lArrays;
typedef Array<rArray>             rArrays;
typedef Array<cArray>             cArrays;

/**
 * Variadic template recurrence starter. For documentation of the function
 * itself see the other "concatenate".
 */
inline rArray concatenate ()
{
    return rArray();
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
template <typename ...Params> rArray concatenate (rArray const & v1, Params ...p)
{
    if (sizeof...(p) == 0)
    {
        return v1;
    }
    else
    {
        rArray v2 = concatenate (p...);
        rArray v (v1.size() + v2.size());
        for (size_t i = 0; i < v1.size(); i++)
            v[i] = v1[i];
        for (size_t i = 0; i < v2.size(); i++)
            v[i + v1.size()] = v2[i];
        return v;
    }
}

// return absolute values
rArray abs (const cArrayView u);
rArrays abs (cArrays const &u);

// boolean aggregation
bool all (const ArrayView<bool> v);
bool any (const ArrayView<bool> v);

// minimal element
template <typename T> T min (const ArrayView<T> a)
{
    T z = a.front();
    for (T const * it = a.begin(); it != a.end(); it++)
        if (*it < z)
            z = *it;
        return z;
}

// maximal element
template <typename T> T max (const ArrayView<T> a)
{
    T z = a.front();
    for (T const * it = a.begin(); it != a.end(); it++)
        if (*it > z)
            z = *it;
        return z;
}

// return per-element power
template <typename T> NumberArray<T> pow (NumberArray<T> const & u, double e)
{
    NumberArray<T> v(u.size());
    
    auto iu = u.begin();
    auto iv = v.begin();
    
    while (iu != u.end())
        *(iv++) = pow(*(iu++), e);
    
    return v;
}
template <typename T> Array<T> pow (Array<T> const & u, double e)
{
    Array<T> v(u.size());
    
    auto iu = u.begin();
    auto iv = v.begin();
    
    while (iu != u.end())
        *(iv++) = pow(*(iu++), e);
    
    return v;
}
template <class T> NumberArray<T> sqrt (NumberArray<T> const & A)
{
    size_t N = A.size();
    NumberArray<T> B (N);
    
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
template <typename T> T sum (const ArrayView<T> v)
{
    return std::accumulate(v.begin(), v.end(), T(0));
}

// summation of nested arrays
template <typename T> NumberArray<T> sums (const ArrayView<NumberArray<T>> v)
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
Array<bool> operator == (const ArrayView<T> u, T x)
{
    Array<bool> v(u.size());
    for (size_t i  = 0; i < u.size(); i++)
        v[i] = (u[i] == x);
    return v;
}

template <typename T>
Array<bool> operator == (const ArrayView<T> u, const ArrayView<T> v)
{
    assert(u.size() == v.size());
    
    Array<bool> w(u.size());
    for (size_t i  = 0; i < u.size(); i++)
        w[i] = (u[i] == v[i]);
    return w;
}

inline bool all(const ArrayView<bool> B)
{
    bool ok = true;
    for (bool b : B)
        ok = ok and b;
    return ok;
}

inline bool any(const ArrayView<bool> B)
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
void eval (TFunctor f, TArray grid, TArray& vals)
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
template <typename T> NumberArray<T> join (const ArrayView<NumberArray<T>> arrays)
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
        if (++repeat <= n)
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
