//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2018, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#ifndef HEX_ARRAYS
#define HEX_ARRAYS

// --------------------------------------------------------------------------------- //

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <numeric>
#include <typeinfo>
#include <type_traits>
#include <vector>

// --------------------------------------------------------------------------------- //

#include "hex-hdffile.h"
#include "hex-memory.h"
#include "hex-misc.h"
#include "hex-numbers.h"

// --------------------------------------------------------------------------------- //

// Forward declaration of Array (unaligned array of items).
template < class T, class Alloc_ = PlainAllocator<T> > class Array;

// Forward declaration of NumberArray (aligned array of numbers -> SIMD).
template < class T, class Alloc_ = AlignedAllocator<T,SIMD_VECTOR_BYTES> > class NumberArray;

// --------------------------------------------------------------------------------- //

/**
 * @brief Array view.
 * 
 * This class holds a shallow copy of an array of items of type T.
 * The view is represented by a pointer to the first item and by the total
 * number of elements. It is offers a comfortable means for work with
 * subarrays. Also, it is a base class for derived data types Array
 * and CLArrayView, that add some advanced functionality to the plain
 * pointer-and-size storage.
 * 
 * ArrayView neither allocates nor deallocates any memory.
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
        std::size_t N_;

        /// Pointer to the array.
        T * array_;

    public:

        // constructor
        ArrayView ()
            : N_(0), array_(nullptr) {}

        // construct from pointer and size
        ArrayView (std::size_t n, T * ptr)
            : N_(n), array_(ptr) {}

        // copy constructor
        ArrayView (ArrayView<T> const & a, std::size_t i = 0, std::size_t n = 0)
            : N_((n > 0) ? n : a.size()), array_(const_cast<T*>(a.data()) + i) { assert(i + N_ <= a.size()); }

        // construct view from Array const lvalue reference
        ArrayView (Array<T> const & a, std::size_t i = 0, std::size_t n = 0)
            : N_((n > 0) ? n : a.size()), array_(const_cast<T*>(a.data()) + i) { assert(i + N_ <= a.size()); }

        // construct view from NumberArray const lvalue reference
        ArrayView (NumberArray<T> const & a, std::size_t i = 0, std::size_t n = 0)
            : N_((n > 0) ? n : a.size()), array_(const_cast<T*>(a.data()) + i) { assert(i + N_ <= a.size()); }

        // construct from consecutive memory segment
        ArrayView (const_iterator i, const_iterator j)
            : N_(j - i), array_(const_cast<T*>(&(*i))) {}

        // construct from right-value reference
        ArrayView (ArrayView<T> && r)
            : ArrayView() { std::swap(N_, r.N_); std::swap(array_, r.array_); }

        // destructor
        ~ArrayView () {}

        /// Change view.
        void reset (std::size_t n, T const * ptr)
        {
            N_ = n;
            array_ = const_cast<T*>(ptr);
        }

        /// Assignment operator.
        ArrayView<T> & operator = (const ArrayView<T> v)
        {
            if (v.size() != size())
                HexException ("[ArrayView::operator=] Cannot copy %ld elements to %ld fields!", v.size(), N_);

            for (std::size_t i = 0; i < size(); i++)
                array_[i] = v[i];

            return *this;
        }

        /// Element-wise access (non-const).
        T & operator[] (std::size_t i)
        {
#ifdef NDEBUG
            return array_[i];
#else
            // bounds check
            if (i < size())
                return array_[i];
            else
                HexException("[ArrayView::operator[]] Index %ld out of bounds (size = %ld) !", i, N_);
#endif
        }

        /// Element-wise access (const).
        T const & operator[] (std::size_t i) const
        {
#ifdef NDEBUG
            return array_[i];
#else
            if (i < size())
                return array_[i];
            else
                HexException("[ArrayView::operator[]] Index %ld out of bounds (size = %ld) !", i, N_);
#endif
        }

        /// Length of the array (number of elements).
        std::size_t size () const { return N_; }

        /// Pointer to the data.
        //@{
        virtual T * data () { return array_; }
        virtual T const * data () const { return array_; }
        //@}

        //
        // STL-like iterator interface
        //

        iterator begin ()                   { return data(); }
        const_iterator begin () const       { return data(); }
        iterator end ()                     { return data() + size(); }
        const_iterator end () const         { return data() + size(); }
        T & front (int i = 0)               { return *(data() + i); }
        T const & front (int i = 0) const   { return *(data() + i); }
        T & back (int i = 0)                { return *(data() + size() - 1 - i); }
        T const & back (int i = 0) const    { return *(data() + size() - 1 - i); }

        /// Fill the array with a value.
        void fill (T x) { for (T & y : *this) y = x; }

        /// Check whether the size is equal to zero.
        bool empty () const { return size() == 0; }

        /// Square of the two-norm (defined only for scalar data type).
        template <class = typename std::enable_if<is_scalar<T>::value>> double sqrnorm () const
        {
            double sqrnorm = 0.;
            for (T const & x : *this)
                sqrnorm += sqrabs(x);
            return sqrnorm;
        }

        /// Two-norm (defined only for scalar data type).
        template <class = typename std::enable_if<is_scalar<T>::value>> double norm () const
        {
            return std::sqrt(sqrnorm());
        }

        /// Return sub-view.
        ArrayView<T> slice (std::size_t i, std::size_t j) const
        {
            return ArrayView<T>(begin() + i, begin() + j);
        }

        /**
         * @brief Output to a set of strings.
         * 
         * This will print the contents of the array to a set of strings
         * of given maximal length. Elements are delimited by commas.
         */
        std::vector<std::string> lines (unsigned len) const
        {
            std::stringstream line, tmp;
            std::vector<std::string> output;

            for (std::size_t i = 0; i < size(); i++)
            {
                // if line is empty, always add the next element
                if (line.str().empty())
                    line << '[' << (*this)[i];

                // otherwise check length first
                else
                {
                    // use temporary stream
                    tmp.str("");
                    tmp << line.str() << (*this)[i];

                    // if the new element fits within the maximal length, add it
                    if (tmp.str().size() + 1 <= len)
                        line << (*this)[i];

                    // otherwise flush the line and start a new one
                    else
                    {
                        output.push_back(line.str());
                        line.str("");
                        line << (*this)[i];
                    }
                }

                // add comma, if this is not the last element in the array
                if (i != size() - 1)
                    line << ',';

                // otherwise flush this last line
                else
                {
                    line << ']';
                    output.push_back(line.str());
                }
            }

            return output;
        }
};

// --------------------------------------------------------------------------------- //

/**
 * @brief A comfortable data array class.
 * 
 * Class Array is intended as a Hex's replacement for std::vector\<T\>.
 * Properties:
 * - User specified allocator (given as the second template argument
 *   of the class). See @ref PlainAllocator and @ref AlignedAllocator
 *   for examples of such allocators. Allocator is expected to have two
 *   static methods: void* alloc(size_t) that returns a pointer to new memory
 *   chunk and void free() which deallocates the pointer.
 * - Basic iterator interface simillar to STL containers -- methods begin(),
 *   end(), and thus also the ability to appear in the range-based for loops.
 */
template <class T, class Alloc_> class Array : public ArrayView<T>
{
    public:

        // aliases
        typedef T DataType;
        typedef T * iterator;
        typedef T const * const_iterator;
        typedef Alloc_ Alloc;

    public:

        // constructors, creates an empty array
        Array ()
            : ArrayView<T>() {}

        // constructor, creates a length-n "x"-filled array
        explicit Array (std::size_t n, T x = T(0))
            : ArrayView<T>(n, Alloc::alloc(n)) { for (std::size_t i = 0; i < size(); i++) (*this)[i] = x; }

        // constructor, copies a length-n "array
        explicit Array (std::size_t n, T const * const x)
            : ArrayView<T>(n, Alloc::alloc(n)) { for (std::size_t i = 0; i < size(); i++) (*this)[i] = x[i]; }

        // copy constructor from Array const lvalue reference
        Array (Array<T> const & a)
            : ArrayView<T>(a.size(), Alloc::alloc(a.size())) { for (std::size_t i = 0; i < size(); i++) (*this)[i] = a[i]; }

        // copy constructor from const ArrayView
        Array (const ArrayView<T> a)
            : ArrayView<T>(a.size(), Alloc::alloc(a.size())) { for (std::size_t i = 0; i < size(); i++) (*this)[i] = a[i]; }

        // copy constructor from Array rvalue reference
        Array (Array<T> && a)
            : ArrayView<T>() { std::swap (ArrayView<T>::N_, a.ArrayView<T>::N_); std::swap (ArrayView<T>::array_, a.ArrayView<T>::array_); }

        // copy constructor from std::vector
        Array (std::vector<T> const & a)
            : ArrayView<T>(a.size(), Alloc::alloc(a.size())) { for (std::size_t i = 0; i < size(); i++) (*this)[i] = a[i]; }

        // copy constructor from initializer list
        Array (std::initializer_list<T> a)
            : ArrayView<T>(a.end()-a.begin(), Alloc::alloc(a.end()-a.begin())) { std::size_t i = 0; for (auto it = a.begin(); it != a.end(); it++) (*this)[i++] = *it; }

        // copy constructor from two forward iterators
        template <typename ForwardIterator> Array (ForwardIterator i, ForwardIterator j)
            : ArrayView<T>(std::distance(i,j), Alloc::alloc(std::distance(i,j)))
        {
            std::size_t n = 0;
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

        /// Return number of elements.
        std::size_t size () const { return ArrayView<T>::size(); }

        /**
         * @brief Resize array.
         * 
         * The method will change the length of the array, most
         * probably by reallocating the storage. The possible new
         * elements will be initialized to T(0). The function is
         * virtual, so that is can be safely overridden in the
         * derived classes.
         */
        virtual std::size_t resize (std::size_t n)
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
            for (std::size_t i = 0; i < n; i++)
                new_array[i] = (i < size()) ? (*this)[i] : T(0);

            Alloc::free (ArrayView<T>::array_);

            ArrayView<T>::N_ = n;
            ArrayView<T>::array_ = new_array;

            return size();
        }

        /// Element-wise access (non-const).
        inline T & operator[] (std::size_t i) { return ArrayView<T>::operator[](i); }

        /// Element-wise access (const).
        inline T const & operator[] (std::size_t i) const { return ArrayView<T>::operator[](i); }

        /// Data pointer.
        //@{
        virtual T* data () { return ArrayView<T>::data(); }
        virtual T const * data () const { return ArrayView<T>::data(); }
        //@}

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

        /**
         * @brief Append an element to the end.
         * 
         * The function will append the element "a" to the end of the array.
         * This will require allocation of a longer array and copy of the
         * original elements to the new array (using the rvalue reference).
         * The function is declared as virtual so that it can be safely
         * overridden in derived classes.
         */
        virtual void push_back (T const & a)
        {
            T* new_array = Alloc::alloc (size() + 1);
            for (std::size_t i = 0; i < size(); i++)
                new_array[i] = std::move(ArrayView<T>::array_[i]);
            new_array[size()] = a;
            ArrayView<T>::N_++;
            Alloc::free (ArrayView<T>::array_);
            ArrayView<T>::array_ = new_array;
        }

        /**
         * @brief Remove the last element from the array.
         * 
         * The last element of the array will be ignored, which is
         * achieved by decrementing array length. The memory will not be
         * deallocated to save processor time. For this reason, if an array
         * is "erased" by subsequent calls to pop_back, it will still occupy
         * the same memory as before. The deallocation will take place only
         * on resize.
         */
        virtual T pop_back ()
        {
            if (size() > 0)
                return *(data() + (--ArrayView<T>::N_));
            else
                HexException("Array has no element to pop!");
        }

        /**
         * @brief Remove the first element from the array.
         * 
         * Moves all elements one element back.
         */
        virtual T pop_front ()
        {
            if (size() > 0)
            {
                T elem = std::move(data()[0]);

                --ArrayView<T>::N_;
                for (std::size_t i = 0; i < ArrayView<T>::N_; i++)
                    data()[i] = std::move(data()[i+1]);

                return elem;
            }
            else
            {
                HexException("Array has no element to pop!");
            }
        }

        /**
         * @brief Append more items.
         * 
         * The function appends a range of items to the end of the
         * array. A reallocation takes place. The range is specified
         * using input iterators.
         */
        template <class InputIterator> void append
        (
            InputIterator first, InputIterator last
        )
        {
            T* new_array = Alloc::alloc(size() + last - first);
            for (std::size_t i = 0; i < size(); i++)
                new_array[i] = (*this)[i];
            for (InputIterator it = first; it != last; it++)
                new_array[size() + it - first] = *it;
            ArrayView<T>::N_ += last - first;
            Alloc::free (ArrayView<T>::array_);
            ArrayView<T>::array_ = new_array;
        }

        /**
         * @brief Insert item.
         * 
         * Function "insert" inserts new item "x" to the array
         * at position pointed to by the iterator "it". Reallocation
         * always taked place. If the iterator points to the end
         * of the original array, the effect is the same as push_back.
         * Otherwise, the items behind the iterator are shifted to
         * make place for the new item.
         */
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

        /// Check that size equals to zero.
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
            for (std::size_t i = 0; i < size(); i++)
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
            std::swap(ArrayView<T>::N_, b.ArrayView<T>::N_);
            std::swap(ArrayView<T>::array_, b.ArrayView<T>::array_);

            return *this;
        }

        Array<T>& operator &= (Array<T> const & b)
        {
            append(b.begin(), b.end());
        }
};

// --------------------------------------------------------------------------------- //

/**
 * @brief A comfortable number array class.
 * 
 * Class NumberArray is intended as a Hex's replacement for std::vector\<T\> for
 * number type T.
 * Properties:
 * - correct memory alignment of the storage
 * - even number of allocated elements (zero padded) for the use in SSE accelerators
 * - reserved storage reducing necessary reallocations
 * - basic iterator interface (members Array::begin(), Array::end()).
 * - datafile interface (ability to save and load to/from scratch data files)
 * - a collection of overloaded arithmetic operators (sum of two arrays,
 *   difference, multiplication by a number etc.)
 */
template <class T, class Alloc_> class NumberArray : public Array<T, Alloc_>
{
    public:

        // aliases
        typedef T DataType;
        typedef T * iterator;
        typedef T const * const_iterator;
        typedef Alloc_ Alloc;

    protected:

        /// Allocated memory.
        std::size_t Nres_;

        /// HDF file to store to / load from.
        std::string name_;

    public:

        // default constructor, creates an empty array
        NumberArray ()
            : Array<T,Alloc>(), Nres_(size()), name_() {}

        // constructor, creates a length-n "x"-filled array
        explicit NumberArray (std::size_t n, T x = T(0))
            : Array<T,Alloc>(n, x), Nres_(size()), name_() {}

        // constructor, copies a length-n "array
        explicit NumberArray (std::size_t n, T const * x)
            : Array<T,Alloc>(n, x), Nres_(size()), name_() {}

        // copy constructor from ArrayView const lvalue reference
        NumberArray (ArrayView<T> const & a)
            : Array<T,Alloc>(a), Nres_(size()), name_() {}

        // copy constructor from Array const lvalue reference
        NumberArray (NumberArray<T> const & a)
            : Array<T,Alloc>((ArrayView<T> const &)a), Nres_(size()), name_() {}

        // copy constructor from std::vector
        NumberArray (std::vector<T> const & a)
            : Array<T,Alloc>(a), Nres_(size()), name_() {}

        // copy constructor from initializer list
        NumberArray (std::initializer_list<T> a)
            : Array<T,Alloc>(a), Nres_(size()), name_() {}

        // copy constructor from two forward iterators
        template <typename ForwardIterator> NumberArray (ForwardIterator i, ForwardIterator j)
            : Array<T,Alloc>(i,j), Nres_(size()), name_() {}

        // copy constructor from Array rvalue reference
        NumberArray (NumberArray<T> && a)
            : NumberArray<T,Alloc>()
        {
            std::swap(ArrayView<T>::N_, a.ArrayView<T>::N_);
            std::swap(ArrayView<T>::array_, a.ArrayView<T>::array_);
            std::swap(Nres_,a.Nres_);
            std::swap(name_, a.name_);
        }

        // destructor
        ~NumberArray ()
        {
            // Array<T> will do the destruction
            // ... do nothing here
        }

        /// Item count.
        std::size_t size () const { return ArrayView<T>::size(); }

        /**
         * @brief Resize array.
         * 
         * Set size of the array. If the new size is less than the
         * reserved memory, the end of the array is just appropriately
         * shifted. The "new" elements are not initialized.
         * If the size is bigger than the reserve, the reallocation
         * takes place and new elements are initialized to T(0).
         * @return New size of the array.
         */
        virtual std::size_t resize (std::size_t n)
        {
            if (n <= Nres_)
            {
                ArrayView<T>::N_ = n;
                return n;
            }

            Nres_ = n;
            T * new_array = Alloc::alloc(Nres_);

            for (std::size_t i = 0; i < n; i++)
                new_array[i] = (i < size()) ? (*this)[i] : T(0);

            Alloc::free (ArrayView<T>::array_);
            ArrayView<T>::N_ = n;
            ArrayView<T>::array_ = new_array;

            return size();
        }

        /**
         * @brief Reserve memory.
         * 
         * If the requested reserve size is greater than than the current
         * reserve size, the memory is reallocated. Reserved memory allows
         * fast push_back-s and similar operations, because no reallocation
         * is necessary then.
         * @return New reserve size.
         */
        std::size_t reserve (std::size_t n)
        {
            if (n > Nres_)
            {
                Nres_ = n;
                T* new_array = Alloc::alloc(Nres_);

                if (size() > 0)
                {
                    std::memmove(new_array, data(), size() * sizeof(T));
                    Alloc::free (ArrayView<T>::array_);
                }

                ArrayView<T>::array_ = new_array;
            }

            return Nres_;
        }

        /// Element-wise access (non-const).
        inline T & operator[] (std::size_t i) { return ArrayView<T>::operator[](i); }

        /// Element-wise access (const).
        inline T const & operator[] (std::size_t i) const { return ArrayView<T>::operator[](i); }

        /**
         * @brief Data pointer.
         * 
         * Gets pointer to the (possibly aligned) memory.
         */
        //@{
        virtual T * data ()
        {
            return ArrayView<T>::array_;
        }
        virtual T const * data () const
        {
            return ArrayView<T>::array_;
        }
        //@}

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

        /**
         * @brief Add element to beginning.
         * 
         * Prepends a new element to the start of the array. If the reserved
         * storage has sufficient size, no reallocation takes place.
         * The whole array has to be shifted by one element in any case.
         */
        virtual void push_front (T const & a)
        {
            // can we just whift the whole data?
            if (size() + 1 <= Nres_)
            {
                // shift data by one element
                for (std::size_t i = size(); i > 0; i--)
                    (*this)[i] = (*this)[i-1];

                // set the new front element
                (*this)[0] = a;

                // update size of the array
                ArrayView<T>::N_++;
            }
            else
            {
                // reset storage size
                Nres_ = 2 * std::max<std::size_t>(1, Nres_);

                // reallocate
                T* new_array = Alloc::alloc(Nres_);

                // check if we have any previous data
                if (size() > 0 and data() != nullptr)
                {
                    // copy original data
                    std::memmove(new_array + 1, data(), size() * sizeof(T));

                    // destroy original array
                    Alloc::free (ArrayView<T>::array_);
                }

                // use new array
                ArrayView<T>::array_ = new_array;

                // set the new front element
                (*this)[0] = a;
            }
        }

        /**
         * @brief Add element to end.
         * 
         * Appends a new element to the end of the array. If the reserved
         * storage has sufficient size, no reallocation takes place.
         */
        virtual void push_back (T const & a)
        {
            if (size() + 1 > Nres_)
            {
                // double the capacity
                Nres_ = 2 * Nres_ + 1;

               // allocate space
                T* new_array = Alloc::alloc(Nres_);

                // check if we have any previous data
                if (size() > 0 and data() != nullptr)
                {
                    // copy original data
                    std::memmove(new_array, data(), size() * sizeof(T));

                    // destroy original array
                    Alloc::free (ArrayView<T>::array_);
                }

                // use new array
                ArrayView<T>::array_ = new_array;
            }

            // copy new element
            (*this)[ArrayView<T>::N_++] = a;
        }

        /**
         * @brief Remove the last element.
         * 
         * Calls Array::pop_back without any modification.
         * @return Value of the removed element.
         */
        T pop_back()
        {
            return Array<T,Alloc>::pop_back();
        }

        /**
         * @brief Append a range of values at end.
         * 
         * Takes the values specified by two iterators and appends them
         * to the end of the array. If the reserved storage is large enough,
         * no reallocation will take place.
         */
        template <class InputIterator> void append (InputIterator first, InputIterator last)
        {
            if (size() + std::distance(first, last) > Nres_)
            {
                // raise the capacity
                Nres_ += last - first;

                // allocate space
                T* new_array = Alloc::alloc(Nres_);

                // check if we have any previous data
                if (size() > 0 and data() != nullptr)
                {
                    // copy original data
                    std::memmove(new_array, data(), size() * sizeof(T));

                    // destroy original array
                    Alloc::free (ArrayView<T>::array_);
                }

                // use new array
                ArrayView<T>::array_ = new_array;
            }

            // copy new elements
            for (InputIterator it = first; it != last; it++)
                (*this)[ArrayView<T>::N_++] = *it;
        }
        void append (const ArrayView<T> a)
        {
            append(a.begin(), a.end());
        }

        /// Check that size equals to zero.
        bool empty () const
        {
            return size() == 0;
        }

        /// Fill array with zeros.
        void clear ()
        {
            this->fill(0);
        }

        /// Reset array: deallocate everything, resize to zero.
        void drop ()
        {
            Nres_ = ArrayView<T>::N_ = 0;
            Alloc::free (ArrayView<T>::array_);
            ArrayView<T>::array_ = nullptr;
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
            for (std::size_t i = 0; i < size(); i++)
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

        /// Return complex conjugated array.
        NumberArray<T> conj () const
        {
            NumberArray<T> c = *this;
            for (std::size_t i = 0; i < size(); i++)
            {
                Complex z = c[i];
                c[i] = Complex(z.real(), -z.imag());
            }
            return c;
        }

        /// Compute (square of) the usual 2-norm.
        double sqrnorm () const
        {
            double ret = 0.;
            for (T z : *this)
                ret += sqrabs(z);
            return ret;
        }

        /// Compute the usual 2-norm.
        double norm () const
        {
            return std::sqrt(sqrnorm());
        }

        /** 
         * @brief Apply a user transformation.
         * 
         * The functor "f" will be applied on every item and the resulting
         * array is returned. It is expected that the return value of the functor
         * is a number type, so that NumberArray can be used.
         */
        template <class Functor> auto transform (Functor f) const -> NumberArray<decltype(f(T(0)))>
        {
            std::size_t n = size();
            NumberArray<decltype(f(T(0)))> c(n);

            for (std::size_t i = 0; i < n; i++)
                c[i] = f((*this)[i]);

            return c;
        }

        /// Return a subarray using ArrayView.
        ArrayView<T> slice (std::size_t left, std::size_t right) const
        {
            assert (right >= left);

            return ArrayView<T>(begin() + left, begin() + right);
        }

        /**
         * @brief Convert to SQL BLOB.
         * 
         * Converts contents to SQL-readable BLOB (hexadecimal text format)
         * 
         * @warning The data are stored in the endianness of the current machine.
         */
        std::string toBlob () const
        {
            // get byte pointer
            unsigned char const * dataptr = reinterpret_cast<unsigned char const*>(data());

            // get byte count
            std::size_t count = size() * sizeof(T);

            // resulting string
            std::ostringstream hexa;
            hexa << "x'" << std::hex << std::setfill('0');

            // for all bytes
            for (std::size_t i = 0; i < count; i++)
                hexa << std::setw(2) << static_cast<unsigned>(dataptr[i]);

            hexa << "'";
            return hexa.str();
        }

        /**
         * @brief Convert from SQL BLOB.
         * 
         * Decode string from SQL-readable BLOB (hexadecimal format) to correct binary array.
         * 
         * @warning The data are assumed to possess the endianness of the current machine.
         */
        void fromBlob (std::string const & s)
        {
            if (data() != nullptr and size() != 0)
                Alloc::free (ArrayView<T>::array_);

            // the first character outght to be "x" or "X"
            // the second character outght to be "'"
            // the last character ought to be "'" as well
            if ((s[0] != 'x' and s[0] != 'X') or s[1] != '\'' or s.back() != '\'')
                HexException("[NumberArray::fromBlob] Blob has wrong format, %s.", s.c_str());

            // create substring
            std::string ss (s.begin() + 2, s.end() - 1);

            // compute size
            std::size_t bytes = ss.size() / 2;
            ArrayView<T>::N_ = bytes / sizeof(T);

            // allocate space
            ArrayView<T>::array_ = Alloc::alloc(size());

            // for all bytes
            for (std::size_t i = 0; i < bytes; i++)
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

        /**
         * @brief Link to HDF file.
         * 
         * In order to avoid repetitious specifying of the HDF filename,
         * it is possible to link the array to a single file and then use
         * @ref hdfload and @ref hdfsave functions without arguments.
         */
        void link (std::string name)
        {
            name_ = name;
        }

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
        //@{
        bool hdfsave (std::string name, bool docompress = false, int consec = 10) const
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

                if (not zero_blocks.empty() and not hdf.write("zero_blocks", &(zero_blocks[0]), zero_blocks.size()))
                    return false;

                if (not elements.empty() and not hdf.write("array", &(elements[0]), elements.size()))
                    return false;
            }
            else
            {
                if (not hdf.write("array", data(), size()))
                    return false;
            }

            return true;
        }
        bool hdfsave () const
        {
            return hdfsave (name_);
        }
        //@}

        /**
         * @brief Load array from HDF file.
         * 
         * See @ref hdfsave for the expected structure of the HDF file.
         * @param name Filename.
         */
        //@{
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
            std::size_t size = hdf.size("array") / typeinfo<T>::ncmpt;

            // read packed elements
            if (elements.resize(size) > 0 and not hdf.read("array", &(elements[0]), elements.size()))
                return false;

            // remove previous data
            if (data() != nullptr)
            {
                Alloc::free (ArrayView<T>::array_);
                ArrayView<T>::array_ = nullptr;
                ArrayView<T>::N_ = Nres_ = 0;
            }

            // unpack (if necessary)
            if (zero_blocks.size() > 0)
                *this = std::move(elements.decompress(zero_blocks));
            else
                *this = std::move(elements);

            return true;
        }
        bool hdfload ()
        {
            return hdfload (name_);
        }
        //@}

        std::string hdfname () const
        {
            return name_;
        }

        /**
         * @brief Get compressed array.
         * 
         * Omit all consecutive occurences of zero (threshold number of consecutive
         * occurences is specified in "consec") and save the start and (one after) end position
         * of the zero blocks in the array to a new array conventionally called "zero_blocks".
         * 
         * @return Pair "zero_blocks", "compressed_array".
         */
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
            for (std::size_t i = 0; i < size(); i++)
            {
                if ((*this)[i] == T(0))
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
            for (std::size_t iblock = 0; iblock < nonzero_blocks.size()/2; iblock++)
            {
                int start = nonzero_blocks[2*iblock];
                int end = nonzero_blocks[2*iblock+1];
                carray.append(begin() + start, begin() + end);
            }

            return std::make_tuple(zero_blocks,carray);
        }

        /**
         * @brief Decompress array.
         * 
         * This function will create a new array using self and the supplied
         * compression info "zero_blocks". The parameter "zero_blocks" is expected
         * to have an even count of items; every pair of items then specifies start
         * and (one after) end position of a block of zeros in the decompressed
         * array. The decompressed array is returned and the original array left
         * intact.
         */
        NumberArray<T> decompress (NumberArray<int> const & zero_blocks) const
        {
            if (zero_blocks.empty())
                return *this;

            // compute final size
            std::size_t final_size = size();
            for (std::size_t i = 0; i < zero_blocks.size()/2; i++)
                final_size += zero_blocks[2*i+1] - zero_blocks[2*i];

            // resize and clean internal storage
            NumberArray<DataType> unpack(final_size);

            // copy nonzero chunks
            int this_end = 0;   // index of last updated element in "this"
            int load_end = 0;   // index of last used element in "nnz_array"
            for (std::size_t i = 0; i < zero_blocks.size()/2; i++)
            {
                int zero_start = zero_blocks[2*i];
                int zero_end = zero_blocks[2*i+1];

                // append nonzero data before this zero block
                std::memmove
                (
                    &(unpack[0]) + this_end,
                    begin() + load_end,
                    (zero_start - this_end) * sizeof(T)
                );

                // move cursors
                load_end += zero_start - this_end;
                this_end  = zero_end;
            }

            // append remaining data
            std::memmove
            (
                &(unpack[0]) + this_end,
                begin() + load_end,
                (final_size - this_end) * sizeof(T)
            );

            return unpack;
        }
};


// --------------------------------------------------------------------------------- //

/**
 * @brief Helper class returning either reference or new object.
 * 
 * This class is a wrapper around an object or a reference to the object.
 * It serves as a useful return type of functions that either return
 * new data, or just reference existing, so that no duplicate memory
 * needs to be allocated.
 */
template <class T> class TmpNumberArray
{
    public:

        // delete default constructor
        TmpNumberArray () = delete;

        // delete default copy constructor
        TmpNumberArray (TmpNumberArray<T> const & array) = delete;

        // initialize as array reference
        TmpNumberArray (const ArrayView<T> view)
            : owner_(false), view_(view)
        {}

        // initialize as allocated array
        TmpNumberArray (NumberArray<T> && array)
            : owner_(true), array_(array)
        {}

        // move constructor
        TmpNumberArray (TmpNumberArray<T> && array)
            : owner_(array.owner_), view_(array.view_), array_(std::move(array.array_))
        {}

        // access view of the underlying data
        ArrayView<T> operator() () const
        {
            return owner_ ? array_ : view_;
        }

        // access pointer to the underlying data
        T const * ptr () const
        {
            return owner_ ? array_.data() : view_.data();
        }

        // delete default copy operator
        TmpNumberArray & operator= (TmpNumberArray const & other) = delete;

    private:

        bool owner_;
        ArrayView<T> view_;
        NumberArray<T> array_;
};

// --------------------------------------------------------------------------------- //

template <class T> class BlockArray
{
    public:

        BlockArray ()
            : arrays_(), inmemory_(true), filename_() {}
        BlockArray (std::size_t nblocks, bool inmemory = true, std::string filename = "")
            : arrays_(nblocks), inmemory_(inmemory), filename_(filename) {}
        BlockArray (BlockArray const & x)
            : arrays_(x.arrays_), inmemory_(x.inmemory_), filename_(x.filename_) {}
        BlockArray (BlockArray && x)
            : arrays_(std::move(x.arrays_)), inmemory_(std::move(x.inmemory_)), filename_(std::move(x.filename_)) {}

        ~BlockArray () {}

        BlockArray<T> & operator= (BlockArray<T> const & x)
        {
            arrays_ = x.arrays_;
            inmemory_ = x.inmemory_;
            filename_ = x.filename_;
            return *this;
        }

        BlockArray<T> & operator= (BlockArray<T> && x)
        {
            arrays_ = std::move(x.arrays_);
            inmemory_ = std::move(x.inmemory_);
            filename_ = std::move(x.filename_);
            return *this;
        }

        std::size_t size () const
        {
            return arrays_.size();
        }

        std::size_t size (std::size_t iblock) const
        {
            assert(iblock < arrays_.size());

            if (inmemory_)
            {
                return arrays_[iblock].size();
            }
            else
            {
                // returns length of the "array" dataset within the on-disk file
                // - assumes no compression
                HDFFile hdf (subname(iblock), HDFFile::readonly);
                return hdf.size("array") / typeinfo<T>::ncmpt;
            }
        }

        void drop (std::size_t iblock)
        {
            assert(iblock < arrays_.size());
            arrays_[iblock].drop();
        }

        NumberArray<T> * begin () { return arrays_.begin(); }
        NumberArray<T> * end () { return arrays_.end(); }

        NumberArray<T> const * begin () const { return arrays_.begin(); }
        NumberArray<T> const * end () const { return arrays_.end(); }

        NumberArray<T> const * cbegin () const { return arrays_.cbegin(); }
        NumberArray<T> const * cend () const { return arrays_.cend(); }

        NumberArray<T> & operator[] (std::size_t i)
        {
            assert(i < arrays_.size());
            return arrays_[i];
        }

        NumberArray<T> const & operator[] (std::size_t i) const
        {
            assert(i < arrays_.size());
            return arrays_[i];
        }

        /**
         * @brief Retrieve a segment of the block array.
         * 
         * The returned object contains either a reference to the requested data, if those are
         * present in the memory, or a copy of the requested data, if those are kept on disk.
         * 
         * @param iblock Index of the block.
         * @param seg Position of the segment within the block.
         * @param nseg Number of elements in the segment.
         */
        TmpNumberArray<T> segment (std::size_t iblock, std::size_t seg, std::size_t nseg) const
        {
            // return a view to existing data ...
            if (inmemory_)
                return TmpNumberArray<T>(ArrayView<T>(arrays_[iblock], seg, nseg));

            // ... or load data from disk
            else
                return TmpNumberArray<T>(std::move(hdfread(subname(iblock),seg, nseg)));
        }

        void setSegment (std::size_t iblock, std::size_t offset, std::size_t n, const ArrayView<T> data) const
        {
            if (data.size() != n)
                HexException("Incompatible dimensions %ld != %ld.", n, data.size());

            // update a view of existing data ...
            if (inmemory_)
                ArrayView<T>(arrays_[iblock], offset, n) = data;

            // ... or save data from disk
            else
                hdfwrite(subname(iblock), offset, n, data.data());
        }

        bool hdfcheck (std::size_t iblock)
        {
            return HDFFile(subname(iblock), HDFFile::readonly).valid();
        }

        bool hdfload (std::size_t iblock)
        {
            assert(iblock < arrays_.size());
            return arrays_[iblock].hdfload(subname(iblock));
        }

        bool hdfsave (std::size_t iblock) const
        {
            assert(iblock < arrays_.size());
            if (not arrays_[iblock].hdfsave(subname(iblock)))
                return false;
            return true;
        }

        bool inmemory () const
        {
            return inmemory_;
        }

        std::string subname (std::size_t iblock) const
        {
            return format("%s-%ld.ooc", filename_.c_str(), iblock);
        }

        /**
         * @brief Read part of an array from HDF file.
         * 
         * This function does not support compression.
         * The file must exist.
         */
        static NumberArray<T> hdfread (std::string const & name, std::size_t offset, std::size_t n)
        {
            NumberArray<T> elements(n);

            // open the file
            HDFFile hdf (name, HDFFile::readonly);

            if (not hdf.valid())
                HexException("File \"%s\" can't be opened for reading.", name.c_str());

            // read elements
            if (not hdf.read("array", &(elements[0]), n, offset))
                HexException("Can't read dataset \"array\" from file \"%s\".", name.c_str());

            return elements;
        }

        /**
         * @brief Write part of an array into HDF file.
         * 
         * This function does not support compression.
         * The file must exist.
         */
        static void hdfwrite (std::string const & name, std::size_t offset, std::size_t n, T const * data)
        {
            // open the file
            HDFFile hdf (name, HDFFile::readwrite);

            if (not hdf.valid())
                HexException("File \"%s\" can't be opened for writing.", name.c_str());

            // read elements
            if (not hdf.write("array", data, n, offset))
                HexException("Can't write to dataset \"array\" in file \"%s\".", name.c_str());
        }

    private:

        /// Array blocks.
        Array<NumberArray<T>> arrays_;

        /// Whether all data are always present in memory (or offloaded to disk).
        bool inmemory_;

        /// Scratch disk filename prefix ("-%d.hdf" will be added).
        std::string filename_;
};

// --------------------------------------------------------------------------------- //

typedef BlockArray<Complex> cBlockArray;

// --------------------------------------------------------------------------------- //

// load array arithmetic operators
#include "hex-arrithm.h"

// --------------------------------------------------------------------------------- //

/// Scalar product of two arrays.
template <class T> T operator | (const ArrayView<T> a, const ArrayView<T> b)
{
    // get size; check if sizes match
    std::size_t N = a.size();
    assert(N == b.size());

    // the scalar product
    T result = 0;

    // iterators
    T const * const restrict pa = a.data();
    T const * const restrict pb = b.data();

    // sum the products
    for (std::size_t i = 0; i < N; i++)
        result += pa[i] * pb[i];

    return result;
}

/**
 * @brief Outer product of two arrays.
 * 
 * Returns a new array with the following values:
 * @f[
 * a_1 b_1, a_1 b_2, \dots, a_1, b_n, a_2 b_1, \dots, a_m b_n
 * @f]
 */
template <class T> NumberArray<T> operator ^
(
    const ArrayView<T> a,
    const ArrayView<T> b
)
{
    NumberArray<T> c (a.size() * b.size());

    auto ic = c.begin();

    for (auto a_ : a)
    for (auto b_ : b)
        *(ic++) = a_ * b_;

    return c;
}

/// Output to text stream.
template <typename T> std::ostream & operator << (std::ostream & out, ArrayView<T> const & a)
{
    out << "[";
    for (std::size_t i = 0; i < a.size(); i++)
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
 * @brief Generate uniform grid.
 * 
 * Return a uniform array
 * @f[
 *      a_1, a_2, a_3, \dots, a_n,
 * @f]
 * where @f$ a_1 @f$ is equal to "start", @f$ a_n @f$ is equal to "end"
 * and @f$ n @f$ is equal to "samples". For consecutive elements it holds
 * @f[
 *      a_{k+1} - a_k = \frac{a_n - a_1}{n - 1} \ .
 * @f]
 * 
 * @param start Left boundary and first sample for "samples" > 0.
 * @param end Right boundary and last sample for "samples" > 1.
 * @param samples Sample count.
 */
template <typename T> NumberArray<T> linspace (T start, T end, unsigned samples)
{
    NumberArray<T> space(samples);

    if (samples == 0)
    {
        return space;
    }

    if (samples == 1)
    {
        space[0] = start;
        return space;
    }

    for (unsigned i = 0; i < samples - 1; i++)
        space[i] = start + ((end - start) * T(i)) / T(samples - 1);
    space[samples - 1] = end;

    return space;
}

/**
 * @brief Generate logarithmic grid.
 * 
 * Return a uniformly diverging array
 * @f[
 *      a_1, a_2, a_3, \dots, a_n,
 * @f]
 * where @f$ a_1 @f$ is equal to "start", @f$ a_n @f$ is equal to "end"
 * and @f$ n @f$ is equal to "samples". For consecutive elements it holds
 * @f[
 *      \frac{a_{k+1}}{a_k} = q = \left(\frac{a_n}{a_1}\right)^{1/(n-1)} \ .
 * @f]
 * 
 * @param x0 Left boundary and first sample for "samples" > 0.
 * @param x1 Right boundary and last sample for "samples" > 1.
 * @param samples Sample count.
 */
template <typename T> NumberArray<T> logspace (T x0, T x1, std::size_t samples)
{
    if (x0 <= 0 or x1 <= 0 or x1 < x0)
        HexException("[logspace] It must be 0 < x1 <= x2 !");

    NumberArray<T> grid(samples);

    if (samples == 0)
    {
        return grid;
    }

    if (samples == 1)
    {
        grid[0] = x0;
        return grid;
    }

    for (unsigned i = 0; i < samples; i++)
    {
        grid[i] = x0 * std::pow(x1 / x0, i / T(samples - 1));
    }

    return grid;
}

/**
 * @brief Generate geometric grid.
 * 
 * Return a geometrically increasing array
 * @f[
 *      a_1, a_2, a_3, \dots, a_n,
 * @f]
 * where @f$ a_1 @f$ is equal to "start", @f$ a_n @f$ is equal to "end"
 * and @f$ n @f$ is equal to "samples". For consecutive elements it holds
 * @f[
 *      \frac{a_{k+1} - a_{k}}{a_k - a_{k-1}} = q \ .
 * @f]
 * Altogether is
 * @f[
 *      a_k = a_1 + (a_n - a_1) \frac{1-q^{k-1}}{1-q^{n-1}}
 * @f]
 */
template <typename T> NumberArray<T> geomspace (T x0, T x1, std::size_t samples, double q)
{
    NumberArray<T> grid(samples);

    if (samples == 0)
    {
        return grid;
    }

    if (samples == 1)
    {
        grid[0] = x0;
        return grid;
    }

    for (unsigned i = 0; i < samples; i++)
    {
        grid[i] = x0 + (x1 - x0) * (1. - std::pow(q,i)) / (1. - std::pow(q,samples-1));
    }

    return grid;
}

/**
 * Write array to file. Array will be written as a single column into
 * an ASCII file.
 * @param array The array to write.
 * @param filename Name of the file to create/overwrite.
 */
template <class T> void write_array
(
    const ArrayView<T> array,
    std::string filename
)
{
    std::ofstream fout(filename);
    for (std::size_t i = 0; i < array.size(); i++)
        fout << array[i] << "\n";
    fout.close();
}

template<> void write_array (const ArrayView<Complex> array, std::string filename);

/**
 * Write array to file. Array will be written as a two columns into
 * an ASCII file, first column contains grid information.
 * @param grid The 1D grid labels for the data.
 * @param array The array to write.
 * @param filename Name of the file to create/overwrite.
 */
template <class T1, class T2> void write_array
(
    const ArrayView<T1> grid,
    const ArrayView<T2> array,
    std::string filename
);

/**
 * @brief Write elements.
 * 
 * For all integers from 0 to m-1 call fetch(i) and save the results to
 * a text file as a single column.
 */
template <typename Fetcher> bool write_1D_data (std::size_t m, std::string filename, Fetcher fetch)
{
    std::ofstream f(filename);

    if (f.bad())
        return false;

    for (std::size_t i = 0; i < m; i++)
        f << fetch(i) << std::endl;

    return true;
}

/**
 * @brief Write 2D data to file.
 * 
 * Write 2D data to file. To allow maximum flexibility, only extensions
 * of the data are passed to the function and a functor that will be
 * repeatedly called with coordinate pair for new data element.
 * @param m Row count.
 * @param n Column count.
 * @param filename Filename of the file to create/overwrite.
 * @param fetch Functor with interface
          @code
               double operator() (size_t, size_t);
          @endcode
 * @return Write success indicator (@c false for failure).
 */
template <class Fetcher> bool write_2D_data (std::size_t m, std::size_t n, std::string filename, Fetcher fetch)
{
    std::ofstream f(filename);

    if (f.bad())
        return false;

    for (std::size_t i = 0; i < m; i++)
    {
        for (std::size_t j = 0; j < n; j++)
            f << fetch(i,j) << " ";
        f << std::endl;
    }

    return true;
}

//
// aliases
//

typedef NumberArray<int>          iArray;
typedef NumberArray<int64_t>      lArray;
typedef NumberArray<Real>         rArray;
typedef NumberArray<Complex>      cArray;

typedef ArrayView<int>         iArrayView;
typedef ArrayView<int64_t>     lArrayView;
typedef ArrayView<Real>        rArrayView;
typedef ArrayView<Complex>     cArrayView;

typedef Array<iArray>             iArrays;
typedef Array<lArray>             lArrays;
typedef Array<rArray>             rArrays;
typedef Array<cArray>             cArrays;

/**
 * Variadic template recurrence starter. For documentation of the function
 * itself see the other "concatenate".
 */
template <class T> NumberArray<T> concatenate ()
{
    return NumberArray<T>();
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
template <class T, class ...Params> NumberArray<T> concatenate (NumberArray<T> const & v1, Params ...p)
{
    if (sizeof...(p) == 0)
    {
        return v1;
    }
    else
    {
        NumberArray<T> v2 = concatenate<T>(p...);
        NumberArray<T> v (v1.size() + v2.size());
        for (std::size_t i = 0; i < v1.size(); i++)
            v[i] = v1[i];
        for (std::size_t i = 0; i < v2.size(); i++)
            v[i + v1.size()] = v2[i];
        return v;
    }
}

// return absolute values
rArray abs (const cArrayView u);
rArrays abs (cArrays const &u);

/// Minimal element of array.
template <typename T> T min (const ArrayView<T> a)
{
    T z = a.front();
    for (T const * it = a.begin(); it != a.end(); it++)
        if (*it < z)
            z = *it;
    return z;
}

/// Maximal element of array.
template <typename T> T max (const ArrayView<T> a)
{
    T z = a.front();
    for (T const * it = a.begin(); it != a.end(); it++)
        if (*it > z)
            z = *it;
    return z;
}

#define DEFINE_FUN_1ARR(fun,kern)           \
template <class T> NumberArray<T> fun       \
(                                           \
    const ArrayView<T> in                   \
)                                           \
{                                           \
    std::size_t N = in.size();              \
    NumberArray<T> out(N);                  \
    for (std::size_t i = 0; i < N; i++)     \
        out[i] = kern(in[i]);               \
    return out;                             \
}

#define DEFINE_FUN_1ARR_1DBL(fun,kern)      \
template <class T> NumberArray<T> fun       \
(                                           \
    const ArrayView<T> in,                  \
    double y                                \
)                                           \
{                                           \
    std::size_t N = in.size();              \
    NumberArray<T> out(N);                  \
    for (std::size_t i = 0; i < N; i++)     \
        out[i] = kern(in[i],y);             \
    return out;                             \
}

#define DEFINE_FUN_2ARR(fun,kern)           \
template <class T> NumberArray<T> fun       \
(                                           \
    const ArrayView<T> a,                   \
    const ArrayView<T> b                    \
)                                           \
{                                           \
    assert(a.size() == b.size());           \
    std::size_t N = a.size();               \
    NumberArray<T> out(N);                  \
    for (std::size_t i = 0; i < N; i++)     \
        out[i] = kern(a[i],b[i]);           \
    return out;                             \
}

DEFINE_FUN_1ARR(exp,std::exp)
DEFINE_FUN_1ARR(sin,std::sin)
DEFINE_FUN_1ARR(cos,std::cos)
DEFINE_FUN_1ARR(asin,std::asin)
DEFINE_FUN_1ARR(sqrt,std::sqrt)
DEFINE_FUN_2ARR(atan2,std::atan2)
DEFINE_FUN_1ARR_1DBL(pow,std::pow)

/// Return per-element hypot.
NumberArray<Real> hypot (NumberArray<Real> const & A, NumberArray<Real> const & B);
/// Return per-element square of absolute value.
NumberArray<Real> sqrabs (NumberArray<Complex> const & A);
/// Return per-element real part.
NumberArray<Real> realpart (NumberArray<Complex> const & A);
/// Return per-element imag part.
NumberArray<Real> imagpart (NumberArray<Complex> const & A);

/// Sum elements in array.
template <typename T> T sum (const ArrayView<T> v)
{
    return std::accumulate(v.begin(), v.end(), T(0));
}

/// Sum arrays.
template <typename T> NumberArray<T> sums (const ArrayView<NumberArray<T>> v)
{
    if (v.size() == 0)
        return NumberArray<T>();    // empty array

    return std::accumulate
    (
        v.begin(),
        v.end(),
        NumberArray<T> (v[0].size()),
        [](NumberArray<T> a, NumberArray<T> b) -> NumberArray<T>
        {
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
    for (std::size_t i  = 0; i < u.size(); i++)
        v[i] = (u[i] == x);
    return v;
}

/// Comparison of two arrays.
template <typename T>
Array<bool> operator == (const ArrayView<T> u, const ArrayView<T> v)
{
    assert(u.size() == v.size());

    Array<bool> w(u.size());
    for (std::size_t i  = 0; i < u.size(); i++)
        w[i] = (u[i] == v[i]);
    return w;
}

/// Check that all values are "true".
inline bool all (const ArrayView<bool> B)
{
    bool ok = true;
    for (bool b : B)
        ok = ok and b;
    return ok;
}

/// Check that all values are "true".
inline bool all (const ArrayView<int> B)
{
    bool ok = true;
    for (bool b : B)
        ok = ok and b;
    return ok;
}

/// Check if any value is "true".
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
    std::size_t N = grid.size();
    assert(N == vals.size());

    for (std::size_t i = 0; i < N; i++)
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
template <typename Tidx, typename Tval> void merge
(
    NumberArray<Tidx>       & idx1, NumberArray<Tval>       & arr1,
    NumberArray<Tidx> const & idx2, NumberArray<Tval> const & arr2
)
{
    // positions in arrays
    std::size_t i1 = 0;
    std::size_t i2 = 0;

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

/// Join elements from all subarrays.
template <typename T> NumberArray<T> join (const ArrayView<NumberArray<T>> arrays)
{
    NumberArray<std::size_t> partial_sizes (arrays.size() + 1);

    // get partial sizes
    partial_sizes[0] = 0;
    for (std::size_t i = 0; i < arrays.size(); i++)
        partial_sizes[i+1] = partial_sizes[i] + arrays[i].size();

    // result array
    NumberArray<T> res (partial_sizes.back());

    // concatenate arrays
    for (std::size_t i = 0; i < arrays.size(); i++)
    {
        if (arrays[i].size() > 0)
            ArrayView<T>(res, partial_sizes[i], arrays[i].size()) = arrays[i];
    }

    return res;
}

/// Drop all redundant repetitions from sorted array.
template <class T> NumberArray<T> sorted_unique (const ArrayView<T> v, int n = 1)
{
    // create output array
    NumberArray<T> w (v.size());

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

/// Convert ArrayView to a string.
template <class T> std::string to_string (const ArrayView<T> v, char sep = ' ')
{
    std::ostringstream ss;
    for (std::size_t i = 0; i < v.size(); i++)
    {
        if (i == 0)
            ss << v[i];
        else
            ss << sep << v[i];
    }
    return ss.str();
}

/// Drop small elements of array (replace by zero).
rArray threshold (const rArrayView a, double eps);

/// Combine real and imaginary parts.
cArray interleave (const rArrayView re, const rArrayView im);

/// Join two Arrays.
template <class T> Array<T> operator & (Array<T> const & a, Array<T> const & b)
{
    Array<T> c = a;
    c &= b;
    return c;
}

/// Join two NumberArrays.
template <class T> NumberArray<T> operator & (NumberArray<T> const & a, NumberArray<T> const & b)
{
    NumberArray<T> c = a;
    c &= b;
    return c;
}

/// Return array generated by std::iota
template <class T> NumberArray<T> identity (std::size_t N)
{
    NumberArray<T> A (N);
    std::iota(A.begin(), A.end(), T(0));
    return A;
}

// --------------------------------------------------------------------------------- //

#endif // HEX_ARRAYS
