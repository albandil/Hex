//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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

#ifndef HEX_SPARSE_GRID
#define HEX_SPARSE_GRID

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ndcube.h"

namespace spgrid
{

/**
 * @brief VTK cell shapes for various dimensions.
 * 
 * This list contains VTK cell shapes that are used when exporting the sparse grid
 * into VTK file. The types used are:
 * 
 * dimension | VTK id | meaning
 * --------- | ------ | ----------
 * 0         | 1      | point
 * 1         | 3      | line
 * 2         | 9      | quad
 * 3         | 12     | hexahedron
 */
extern const std::vector<unsigned> vtk_cell;

/**
 * @brief Get further point of a cube along given axis.
 * 
 * Given a subdivision history of a sparse grid point, return a the nearest shifted point along
 * the specified axis.
 * 
 * This function template accepts a "subdivision history" of a point in a sparse grid.
 * The history is a set of bitsets, number of bits in every bitset is equal to the dimension
 * of the grid; every bit of a specific bitset belongs to a single coordinate. The string
 * of axis-specific bits across the bitsets makes a big-endian fractional binary number
 * between 0 and 1 that encodes the position of a point within the (unit) grid.
 * 
 * Example
 * -------
 * 
 * Let the history of a 2-D point looks like this:
 @verbatim
   00, 10, 11
 @endverbatim
 * Then its position is x = 0.011 and y = 0.001 (note once again that these are binary numbers!).
 * The encoded subdivision history of the near cells is then:
 <pre>
   1.0 ┌───────────────┬───────────────┐
       │               │               │
       │               │               │
       │               │               │
       │       ?       │       ?       │
       │               │               │
       │               │               │
       │               │               │
   0.1 ├───────┬───────┼───────────────┤
       │       │       │               │
       │   ?   │   ?   │               │
       │       │       │               │
  0.01 ├───────┼───┬───┤       ?       │
       │       │   │   │               │
 0.001 │   ?   ├───●───┤               │
       │       │   │   │               │
   0.0 └───────┴───┴───┴───────────────┘
      0.0      ↑   ↑  0.1             1.0
               ↑   ↑
               .01 .011
 </pre>
 * The x-shifted point would have coordinate x=0.100, whereas the y-shifted coordinate would be y = 0.010.
 * Note that the number of digits matters; x=0.1 and x=0.100 are the same point but the subdivision
 * history is different.
 */
template <int dim> std::vector<char> shift_fwd (char axis, std::vector<char> pt)
{
    // get bit mask corresponding to the chosen axis
    char mask = special::pow2(axis);
    
    // find the latest history point where the point can be shifted forward
    for (auto bit = pt.rbegin(); bit != pt.rend(); bit++)
    {
        if ((*bit & mask) != 0)
        {
            *bit &= ~mask;  // set bit to 0 (and "carry" to the next binary "digit")
            continue;
        }
        else
        {
            *bit |= mask;  // set bit to 1 (no "carry", done)
            return pt;
        }
    }
    
    HexException("Point cannot be shifted forward.");
}

/**
 * @brief Get previous point of a cube along given axis.
 * 
 * Get the previous (back-shifted) point. See @ref shift_fwd for details on the storage.
 * The back-shifted point from the example in @ref shift_fwd would be x = 0.010 (if the x-axis
 * were chosen) or y = 0.000 (if the y-axis were chosen).
 */
template <int dim> std::vector<char> shift_bck (char axis, std::vector<char> pt)
{
    // get bit mask corresponding to the chosen axis
    char mask = special::pow2(axis);
    
    // find the latest history point where the point can be shifted back
    for (auto bit = pt.rbegin(); bit != pt.rend(); bit++)
    {
        if ((*bit & mask) != 0)
        {
            *bit &= ~mask;  // set bit to 0 (no "carry", done)
            return pt;
        }
        else
        {
            *bit |= mask;  // set bit to 1 (and "carry" to the next binary "digit")
            continue;
        }
    }
    
    HexException("Point cannot be shifted back.");
}

/**
 * @brief Blind template.
 * 
 * This template should never be called (it will throw an exception otherwise).
 * Its explicit specializations for 1-D, 2-D and 3-D are expected to be called instead.
 */
template <int dim> std::vector<std::vector<char>> extrude_point_to_regular_VTK_cell (std::vector<char> const & pt)
{
    HexException("Regular extrusion of point into more than 3 dimensions is not implemented. Please do not use VTK grid export here.");
    return std::vector<std::vector<char>>(); // this is only for the compiler not to complain
}

/**
 * @brief Make line segment from origin point.
 * 
 * The function will return set of points (in the form of subdivision history, see @ref shift_fwd for details)
 * making a line, starting from given origin.
 */
template<> inline std::vector<std::vector<char>> extrude_point_to_regular_VTK_cell<1> (std::vector<char> const & pt)
{
    // order required by VTK: (xmin) (xmax)
    std::vector<std::vector<char>> cell_points = { pt };
    cell_points.push_back(shift_fwd<1>(0,pt));
    return cell_points;
}

/**
 * @brief Make square from origin point.
 * 
 * The function will return set of points (in the form of subdivision history, see @ref shift_fwd for details)
 * making a square, starting from given origin. The ordering of the vertices is chosen such that it will
 * make a valid VTK_QUAD specification.
 */
template<> inline std::vector<std::vector<char>> extrude_point_to_regular_VTK_cell<2> (std::vector<char> const & pt)
{
    // order required by VTK: (xmin,ymin) (xmax,ymin) (xmax,ymax) (xmin,ymax)
    std::vector<std::vector<char>> cell_points = { pt };
    cell_points.push_back(shift_fwd<2>(0,pt));
    cell_points.push_back(shift_fwd<2>(1,cell_points.back()));
    cell_points.push_back(shift_bck<2>(0,cell_points.back()));
    return cell_points;
}

/**
 * @brief Make cube from origin point.
 * 
 * The function will return set of points (in the form of subdivision history, see @ref shift_fwd for details)
 * making a cube, starting from given origin. The ordering of the vertices is chosen such that it will
 * make a valid VTK_HEXAHEDRON specification.
 */
template<> inline std::vector<std::vector<char>> extrude_point_to_regular_VTK_cell<3> (std::vector<char> const & pt)
{
    // order required by VTK: (xmin,ymin,zmin) (xmax,ymin,zmin) (xmax,ymax,zmin) (xmin,ymax,zmin)
    //                        (xmin,ymin,zmax) (xmax,ymin,zmax) (xmax,ymax,zmax) (xmin,ymax,zmax)
    std::vector<std::vector<char>> cell_points = { pt };
    cell_points.push_back(shift_fwd<3>(0,pt));
    cell_points.push_back(shift_fwd<3>(1,cell_points.back()));
    cell_points.push_back(shift_bck<3>(0,cell_points.back()));
    cell_points.push_back(shift_fwd<3>(2,pt));
    cell_points.push_back(shift_fwd<3>(0,cell_points.back()));
    cell_points.push_back(shift_fwd<3>(1,cell_points.back()));
    cell_points.push_back(shift_bck<3>(0,cell_points.back()));
    return cell_points;
}

/**
 * @brief Available sparse grid types.
 * 
 * The sparse grids can be downloaded from the page http://sparse-grids.de/#Nodes.
 */
typedef enum
{
    /* dim = 1 */
    d1l1n1, d1l2n3, d1l3n3, d1l4n7, d1l5n7, /*d1l6n7, d1l7n15, d1l8n15,
    d1l9n15, d1l10n15, d1l11n15, d1l12n15, d1l13n31, d1l14n31, d1l15n31,
    d1l16n31, d1l17n31, d1l18n31, d1l19n31, d1l20n31, d1l21n31, d1l22n31,
    d1l23n31, d1l24n31, d1l25n63,*/
    
    /* dim = 2 */
    d2l1n1, d2l2n5, d2l3n9, d2l4n17, d2l5n33, d2l6n33, d2l7n65, /*d2l8n97,
    d2l9n97, d2l10n161, d2l11n161, d2l12n161, d2l13n257, d2l14n321,
    d2l15n321, d2l16n449, d2l17n449, d2l18n449, d2l19n705, d2l20n705,
    d2l21n705, d2l22n705, d2l23n705, d2l24n705, d2l25n1025,*/
    
    /* dim = 3 */
    d3l1n1, d3l2n7, d3l3n19, d3l4n39, d3l5n87, d3l6n135, /*d3l7n207, d3l8n399,
    d3l9n495, d3l10n751, d3l11n1135, d3l12n1135, d3l13n1759, d3l14n2335,
    d3l15n2527, d3l16n3679, d3l17n4447, d3l18n4447, d3l19n6495, d3l20n8031,
    d3l21n8031, d3l22n11103, d3l23n11103, d3l24n11103, d3l25n15039,*/
    
    /* dim = 4 */
    d4l1n1, d4l2n9, d4l3n33, d4l4n81, d4l5n193, /*d4l6n385, d4l7n641, d4l8n1217,
    d4l9n1985, d4l10n2881, d4l11n4929, d4l12n6465, d4l13n8705, d4l14n13697,
    d4l15n16001, d4l16n22401, d4l17n31617, d4l18n34689, d4l19n47489,
    d4l20n63873,*/
    
    /* dim = 5 */
    d5l1n1, d5l2n11, d5l3n51, d5l4n151, d5l5n391, /*d5l6n903, d5l7n1743,
    d5l8n3343, d5l9n6223, d5l10n10063, d5l11n17103, d5l12n27343, d5l13n38303,
    d5l14n60703,*/
    
    /* dim = 6 */
    d6l1n1, d6l2n13, d6l3n73, d6l4n257, d6l5n737,
}
SparseGridId;

/**
 * @brief Sparse grid nodes.
 * 
 * Structure containing nodes and weights for a given dimension
 * and point count.
 */
typedef struct
{
    int dim;
    int n;
    std::vector<double> x;
    std::vector<double> w;
}
SparseGridNodes;

/// Sparse grid nodes.
extern const std::map<SparseGridId,SparseGridNodes> nodes;

template <class T, int dim> class Domains
{
    public:
        
        // default constructor
        Domains (int lvl = 0) : lvl_(lvl), N_(0), file_(nullptr)
        {
#ifdef _OPENMP
            omp_init_lock(&lock_);
#endif
        }
        
        // move constructor
        Domains (Domains && D) : Domains()
        {
            std::swap(lvl_,      D.lvl_);
            std::swap(N_,        D.N_);
            std::swap(nbytes_,   D.nbytes_);
            std::swap(data_,     D.data_);
            std::swap(file_,     D.file_);
            std::swap(filename_, D.filename_);
        }
        
        // destructor
        ~Domains ()
        {
#ifdef _OPENMP
            omp_destroy_lock(&lock_);
#endif
            if (file_ != nullptr)
            {
                std::fclose(file_);
                std::remove(filename_);
            }
        }
        
        // interchange data structure
        typedef struct
        {
            bool set;
            T val;
            geom::ndCube<dim> cube;
        } Domain;
        
        // subdivision level
        int level () const { return lvl_; }
        
        // number of domains
        std::size_t size () const { return N_; }
        
        // add new subdomain (OpenMP thread safe)
        void add (Domain const & dom)
        {
            // reset cube parameters
            if (N_ == 0)
            {
                lvl_ = dom.cube.level();
                nbytes_ = sizeof(bool) + sizeof(T) + lvl_ * sizeof(char);
            }
            
            assert (dom.cube.level() == lvl_);
            
            // serialize the item to add
            std::vector<char> data (nbytes_);
            pack_(dom, data.data());
            
#ifdef _OPENMP
            // enable only a single thread at once to access the data
            omp_set_lock(&lock_);
#endif
            
            // check if there is already an open scratch file
            if (file_ != nullptr)
            {
                // write data to the end of scratch file
                std::fseek(file_, 0, SEEK_END);
                std::fwrite(data.data(), nbytes_, 1, file_);
            }
            else
            {
                try
                {
                    // try to add new element to memory
                    data_.insert(data_.end(), data.data(), data.data() + nbytes_);
                }
                catch (std::bad_alloc const & ex)
                {
                    // -> not enough memory, we need to create a disk file!
                    
                    // create a unique scratch file name
                    do
                    {
                        static const char templt[] = "pwba2-XXXXXXXXX.tmp";
                        static const char alphnm[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
                        
                        std::srand(std::time(nullptr));
                        std::strncpy(filename_, templt, sizeof(filename_));
                        for (unsigned i = 0; i < std::strlen(filename_); i++)
                            if (filename_[i] == 'X')
                                filename_[i] = alphnm[std::rand() % std::strlen(alphnm)];
                    }
                    while (std::ifstream(filename_).good());
                    
                    // open the file
                    file_ = std::fopen(filename_, "wb+");
                    
                    // check that the file has been created successfully
                    if (file_ == nullptr)
                    {
                        std::perror(format("Error while opening \"%s\": ", filename_).c_str());
                        std::exit(EXIT_FAILURE);
                    }
                    
                    // write existing data to the file
                    std::fwrite(data_.data(), nbytes_, N_, file_);
                    
                    // destroy memory data
                    data_.clear();
                    
                    // finally add the new element
                    std::fwrite(data.data(), nbytes_, 1, file_);
                }
            }
            
            // update subdomain count
            N_++;
            
#ifdef _OPENMP
            // allow other threads to access the data
            omp_unset_lock(&lock_);
#endif
        }
        
        // change values of the domain attributes (OpenMP thread safe)
        void update (std::size_t i, Domain const & dom)
        {
            // serialize the item to write
            std::vector<char> data (nbytes_);
            pack_(dom, data.data());
            
#ifdef _OPENMP
            // enable only a single thread at once to access the data
            omp_set_lock(&lock_);
#endif
            
            // check if there is already an open scratch file
            if (file_ != nullptr)
            {
                std::fseek(file_, i * nbytes_, SEEK_SET);
                std::fwrite(data.data(), nbytes_, 1, file_);
            }
            else
            {
                std::memcpy(data_.data() + i * nbytes_, data.data(), nbytes_);
            }
            
#ifdef _OPENMP
            // allow other threads to access the data
            omp_unset_lock(&lock_);
#endif
        }
        
        // access subdomain attributes
        Domain operator[] (std::size_t i)
        {
            Domain ret;
            
#ifdef _OPENMP
            // enable only a single thread at once to access the data
            omp_set_lock(&lock_);
#endif
            // packed data
            std::vector<char> data (nbytes_);
            
            // check if there is already an open scratch file
            if (file_ != nullptr)
            {
                std::fseek(file_, i * nbytes_, SEEK_SET);
                std::fread(data.data(), nbytes_, 1, file_);
            }
            else
            {
                std::memcpy(data.data(), data_.data() + i * nbytes_, nbytes_);
            }
            
#ifdef _OPENMP
            // allow other threads to access the data
            omp_unset_lock(&lock_);
#endif
            
            // deserialize
            unpack_(ret, data.data());
            
            // return
            return ret;
        }
        
    private:
        
        // current subdivision level
        int lvl_;
        
        // number of integration sub-domains currently stored
        std::size_t N_;
        
        // size of serialized interchange structure
        std::size_t nbytes_;
        
        // data storage in memory
        std::vector<char> data_;
        
        // data storage on disk
        FILE* file_;
        
        // name of the disk file
        char filename_[20];
        
#ifdef _OPENMP
        // concurrency lock
        omp_lock_t lock_;
#endif
        
        void pack_ (Domain const & dom, char * data)
        {
            *reinterpret_cast<bool*>(data) = dom.set;
            *reinterpret_cast<T*>(data + sizeof(bool)) = dom.val;
            if (lvl_ > 0)
                std::memcpy(data + sizeof(bool) + sizeof(T), dom.cube.history(), lvl_);
        }
        
        void unpack_ (Domain & dom, char const * data)
        {
            dom.set = *reinterpret_cast<bool const *>(data);
            dom.val = *reinterpret_cast<T const*>(data + sizeof(bool));
            
            geom::ndCube<dim> c(lvl_, data + sizeof(bool) + sizeof(T));
            dom.cube = c;
        }
};

/**
 * @brief Globally adaptive sparse-grid-based integrator.
 * 
 * This class implements a multidimensional integator. The main routine
 * is @ref integrate_adapt, which accepts a functor to integrate
 * (of signature T (*)(int, double const*)) and a pair of sparse grids
 * (two different-order rules) which will be used to estimate the
 * integration error.
 * 
 * The integrator has to be used for unit cube integration domain.
 * 
 * The method works in the following way: The integral estimate
 * is computed using two different-order rules to estimate also the
 * error. If the error satisfies convergence criteria, integration
 * terminates. Otherwise all non-converged domains will be exploded
 * into 2^dim sub-cubes. The converged cells' estimates are summed
 * and only the non-converged cells are kept for further processing.
 * 
 * See also description of the function @ref integrate_adapt.
 */
template <class T> class SparseGrid
{
    private:
        
        /// Result of the last integration.
        T result_;
        
        /// Error of the last integration.
        T error_;
        
        /// Number of evaluations of the last integration.
        std::size_t neval_;
        
        /// Number of cells of all levels used in the last integration.
        std::size_t ncells_;
        
        /// Absolute tolerance for low- and high-order rule difference.
        double epsabs_;
        
        /// Relative tolerance for low- and high-order rule difference.
        double epsrel_;
        
        /// Absolute tolerance for the contribution of a sub-domain.
        double global_epsabs_;
        
        /// Relative tolerance for the contribution of a sub-domain.
        double global_epsrel_;
        
        /// Minimal subdivision level for integration termination.
        int minlevel_;
        
        /// Maxima subdivision level for integration termination.
        int maxlevel_;
        
        /// Status of the last integration.
        bool ok_;
        
        /// Whether do output diagnosti information (not used at the moment).
        bool verbose_;
        
        /// Verbose outpu prefix.
        std::string prefix_;
        
        /// Whether to write final grid to VTK (only usable when dim <= 3).
        bool write_vtk_;
        
        /// VTK output file name.
        std::string vtkfilename_;
        
        /// Whether to process cells in parallel.
        bool parallel_;
    
    public:
        
        /// Constructor.
        SparseGrid ()
            : epsabs_(1e-8), epsrel_(1e-8), global_epsabs_(1e-8),
              global_epsrel_(1e-8), minlevel_(0), maxlevel_(5),
              ok_(true), verbose_(false), prefix_(""), write_vtk_(false),
              vtkfilename_("spgrid.vtk"), parallel_(false)
        {}
        
        /**
         * @brief Integrand interface.
         * 
         * @param n Number of evaluation points.
         * @param dim Number of coordinates for every evaluation point.
         * @param origin Coordinates of the corner of the integration area that is closest to the origin.
         * @param range Size of integration area.
         * @param scale Fractional positions of the evaluation points within the integration area.
         * @param eval Output array for evaluations.
         * 
         * Evaluates the integrated function at points x[ι] = origin[i] + range * scale[i].
         */
        typedef void (*integrand_t) (int n, int dim, double const * origin, double range, const double * scale, T * eval);
        
        /// Return the result value computed before.
        T result () const { return result_; }
        
        /// Return the error estimate computed before (not implemented yet).
        T error () const { return error_; }
        
        /// Get evaluation count from the last integration.
        std::size_t evalcount () const { return neval_; }
        
        /// Get number of cells of all subdivision levels that have been used in the last integration.
        std::size_t cellcount () const { return ncells_; }
        
        /// Check that the last integration has been successful.
        bool ok () const { return ok_; }
        
        /**
         * @brief Set local absolute tolerance.
         * 
         * If the absolute change of integral estimates between the rough and fine
         * quadrature rule is smaller than the local relative tolerance set by this function,
         * the cell is considered converged.
         */
        void setLocEpsAbs (double epsabs) { epsabs_ = epsabs; }
        
        /**
         * @brief Set local relative tolerance.
         * 
         * If the relative change of integral estimates between the rough and fine
         * quadrature rule is smaller than the local relative tolerance set by this function,
         * the cell is considered converged.
         */
        void setLocEpsRel (double epsrel) { epsrel_ = epsrel; }
        
        /**
         * @brief Set global absolute tolerance
         * 
         * Every cell will be considered converged (or here actually 'irrelevant') if the
         * extrapolation of its integral over the whole integration domain yields a number
         * that is less than the number set by this function.
         */
        void setGlobEpsAbs (double epsabs) { global_epsabs_ = epsabs; }
        
        /**
         * @brief Set global relative tolerance.
         * 
         * Every cell will be considered converged (or here actually 'irrelevant') if the
         * extrapolation of its integral over the whole integration domain yields
         * a fraction of the total estimate that is less than the number set by this function.
         */
        void setGlobEpsRel (double epsrel) { global_epsrel_ = epsrel; }
        
        /**
         * @brief Set minimal subdivision level.
         * 
         * If 'n' is passed as the argument of this function, the integration domain
         * will be subdivided into 2^n sub-cells before the convergence criteria start
         * to be tested.
         */
        void setMinLevel (int minlevel) { minlevel_ = minlevel; }
        
        /**
         * @brief Set maximal subdivision level.
         * 
         * Inhibits further subdivision of small enough cells, even though they
         * will be still considered non-converged.
         */
        void setMaxLevel (int maxlevel) { maxlevel_ = maxlevel; }
        
        /**
         * @brief Set verbosity level.
         * 
         * If set to true, the integrator will print a large amount of diagnostic
         * information.
         */
        void setVerbose (bool verbose) { verbose_ = verbose; }
        
        /**
         * @brief Set verbose output prefix.
         * 
         * Every line of the verbose output (if enabled by @ref setVerbose) will
         * be prefixed by the string passed to this function. For example a string
         * of spaces can be used to align the verbose output among the other output
         * of a program.
         */
        void setPrefix (std::string prefix) { prefix_ = std::string(prefix); }
        
        /**
         * @brief Set whether to write final grid to VTK.
         * 
         * If the write-out is enabled, every final integration grid will be written
         * to a VTK file so that it can be viewed e.g. in ParaView. This function
         * can be used only for one- to three-dimensional grids, because VTK file
         * format does not support higher dimensions.
         */
        void setWriteVTK (bool write_vtk, std::string vtkfilename = "spgrid.vtk")
        {
            write_vtk_ = write_vtk;
            vtkfilename_ = vtkfilename;
        }
        
        /**
         * @brief Set whether to process cells in parallel.
         * 
         * If the parallelization is enabled, integral estimates for several cells
         * will be evaluated simultaneously. OpenMP is used here, so the number
         * of simultaneous threads is controlled by the environment variable
         * OMP_NUM_THREADS.
         */
        void setParallel (bool parallel) { parallel_ = parallel; }
        
        /**
         * @brief Fixed-order spare-grid quadrature.
         * 
         * This function will integrate the function 'F' over a 'dim'-dimensional cube.
         * The specific grid is specified by the third parameter. Function 'F' is
         * expected to accept two parameters:
         *  - 'n': the integer number of dimensions
         *  - 'ptr': pointer to a field of doubles of length 'n'
         * The return value can be of arbitrary arithmetic type (but consistent with
         * the template parameter of the class).
         * 
         * @param F Function to integrate. Expected signature:
         *            void (*) (int npt, int dim, const double* origin, double range, const double* scale, T* eval)
         *          The parameter 'npt' is number of points to evaluate, 'dim' is the dimension of the integration
         *          space (and hence number of coordinates of every point), 'origin[i] + range * scale[i]' are
         *          the evaluation points and 'eval[i]' are the evaluated values to fill. The number of evaluation
         *          points is always equal to the quadrature's number of points.
         * @param i Sparse grid type.
         * @param domain N-dimensional cube (integration volume) which is a subdivision of the unit cube.
         *               This argument is optional; default is a unit hyper-cube.
         */
        template <int dim, class Functor> T integrate_fixed
        (
            Functor F,
            SparseGridId i,
            geom::ndCube<dim> const & domain = geom::ndCube<dim>()
        ) const
        {
            // check that the chosen sparse grid is available
            assert(nodes.find(i) != nodes.end());
            
            // check dimedions
            assert(dim == nodes.at(i).dim);
            
            // shorthands
            int N = nodes.at(i).n;
            const double * px = nodes.at(i).x.data();
            const double * pw = nodes.at(i).w.data();
            
            // evaluate function at the nodes
            NumberArray<T> eval(N);
            F(N, dim, domain.origin().data(), domain.edge(), px, eval.data());
            
            // evaluate the quadrature rule
            T res = 0;
            double orphaned_weights = 0;
            for (int j = 0; j < N; j++)
            {
                // check that the result's magnitude is finite
                if (not std::isfinite(std::abs(eval[j])))
                {
                    // NO : omit this value as we may have hit an (otherwise integrable) singularity
                    orphaned_weights += pw[j];
                    
                    // print warning
                    Debug << "Warning: Non-numerical evaluation in " << domain << std::endl;
                }
                else
                {
                    // YES : update the estimate
                    res += pw[j] * eval[j];
                }
            }
            
            // return the result (extrapolate omitted nodes)
            return res / (1.0 - orphaned_weights) * domain.volume();
        }
        
        /**
         * @brief Globally adaptive sparse-grid-based quadrature.
         * 
         * This function will integrate the function 'F' a 'dim'-dimensional unit cube.
         * Two sparse grids have to be specified. For every integration sub-cell
         * the integral of the function 'F' is estimated by these two quadrature
         * rules and whenever the change is within the tolerances, the result of the
         * second grid is used as the final estimate. Otherwise when the evaluations
         * differ too much, the volume is further subdivided (up to the limit).
         * The subdivision bisects simultaneously in every dimension making always
         * 2^dim new sub-cells from every original cell.
         * 
         * Function 'F' is expected to accept the following parameters:
         *  void (*) (int npt, int dim, const double* origin, double range, const double* scale, T* eval)
         * The parameter 'npt' is number of points to evaluate, 'dim' is the dimension of the integration
         * space (and hence number of coordinates of every point), 'origin[i] + range * scale[i]' are
         * the evaluation points and 'eval[i]' are the evaluated values to fill. The number of evaluation
         * points is always equal to the quadrature's number of points.
         * 
         * This function uses the function @ref integrate_fixed to evaluate the
         * integral estimates on the two specified grids.
         * 
         * @param F Function to integrate.
         * @param ruleLow Low-order rule for error estimation.
         * @param ruleHigh High-order rule for error estimation.
         */
        template <int dim, class Functor> bool integrate_adapt
        (
            Functor F, SparseGridId ruleLow, SparseGridId ruleHigh
        )
        {
            // synonym for 'dim'-dimensional cube
            typedef geom::ndCube<dim> dCube;
            
            // high-order integral estimate over processed cells
            T finalEstimate = 0;
            
            // low-order integral estimate over to-be-yet-processed cells
            T globalEstimate = integrate_fixed(F, ruleLow, dCube());
            
            // declare root domain
            typename Domains<T,dim>::Domain root;
            root.set = false;
            root.val = globalEstimate;
            root.cube = dCube();
            
            // unprocessed subdomains of the integration domain
            Domains<T,dim> domains;
            domains.add(root);
            
            // VTK debug data (only used when write_vtk_ == true)
            std::map<std::vector<char>,T> vtk;
            
            // initialize evaluations count and number of cells
            neval_ = nodes.at(ruleLow).n;
            ncells_ = 1;
            
            // for all subdivision levels
            for (int level = 0; level <= maxlevel_; level++)
            {
                // also compute the volume (as a fraction of root volume) for use in convergence checking
                double levelVolume = domains.size() * std::pow(0.5,level*dim);
                
                // debug info
                if (verbose_)
                {
                    std::cout << std::endl;
                    std::cout << prefix_ << "Level " << level << std::endl;
                    std::cout << prefix_ << "  initial estimate : " << globalEstimate << std::endl;
                    std::cout << prefix_ << "  cells to process : " << domains.size() << " (";
                    
                    double frac = levelVolume * 100;
                    if (frac == 0)
                        std::cout << "= 0";
                    else if (frac < 1)
                        std::cout << "< 1";
                    else
                        std::cout << "~ " << (int)std::ceil(frac);
                    
                    std::cout << "% of orig. volume)" << std::endl;
                }
                
                // auxiliary variables
                std::size_t abslocal = 0, rellocal = 0, absglobal = 0, relglobal = 0, pos = 0, i_pos = 0;
                std::vector<std::size_t> bisectIDs;
                Timer timer;
                
                // check convergence for all not yet converged sub-domains
                # pragma omp parallel if (parallel_) \
                    default (shared) \
                    firstprivate(i_pos) \
                    shared (pos) \
                    reduction (+:abslocal,absglobal,rellocal,relglobal)
                {
                    // per-thread data
                    T thread_finalEstimate = 0;
                    std::vector<std::size_t> thread_bisectIDs;
                    
                    // for all domains (consecutive access)
                    # pragma omp for schedule (static)
                    for (std::size_t i = 0; i < domains.size(); i++)
                    {
                        // get next comain ID
                        # pragma omp critical
                        i_pos = pos++;
                        
                        // get current domain data
                        typename Domains<T,dim>::Domain dom = domains[i_pos];
                        
                        // compute the fine estimate for this domain
                        T fineEstimate = integrate_fixed(F, ruleHigh, dom.cube);
                        
                        // check absolute local convergence
                        if (std::abs(dom.val - fineEstimate) < epsabs_)
                        {
                            dom.set = true;
                            abslocal++;
                        }
                        
                        // check relative local convergence
                        if (std::abs(dom.val - fineEstimate) / std::abs(fineEstimate) < epsrel_)
                        {
                            dom.set = true;
                            rellocal++;
                        }
                        
                        // check absolute global (extrapolated) convergence
                        if (std::abs(dom.val - fineEstimate) * (levelVolume / dom.cube.volume()) < global_epsabs_)
                        {
                            dom.set = true;
                            absglobal++;
                        }
                        
                        // check relative global (extrapolated) convergence
                        if (std::abs(dom.val - fineEstimate) * (levelVolume / dom.cube.volume()) < global_epsrel_ * std::abs(globalEstimate))
                        {
                            dom.set = true;
                            relglobal++;
                        }
                        
                        // store the better estimate
                        dom.val = fineEstimate;
                        
                        // is the cell converged or was max level reached?
                        if (level >= minlevel_ and (dom.set or level == maxlevel_))
                        {
                            // YES : use the current value as final
                            thread_finalEstimate += dom.val;
                            
                            // YES : save data for optional VTK debug output
                            if (write_vtk_ and dim <= 3)
                            {
                                // compose map key (the history)
                                std::vector<char> key (dom.cube.history(), dom.cube.history() + dom.cube.level());
                                
                                // store the pair "(key,average)" to the map
                                # pragma omp critical
                                vtk[key] = dom.val / dom.cube.volume();
                            }
                        }
                        else
                        {
                            // NO : add domain ID to the list of domains needing bisection
                            thread_bisectIDs.push_back(i_pos);
                        }
                        
                        // update domain in the storage
                        domains.update(i_pos,dom);
                    }
                    
                    # pragma omp critical
                    {
                        // update shared variables by contributions from this thread
                        finalEstimate += thread_finalEstimate;
                        bisectIDs.insert(bisectIDs.end(), thread_bisectIDs.begin(), thread_bisectIDs.end());
                    }                
                }
                
                // update evaluation count
                neval_ += domains.size() * nodes.at(ruleHigh).n;
                
                // debug info
                if (verbose_)
                {
                    std::cout << prefix_ << "  processing time : " << timer.nice_time() << std::endl;
                    
                    // print how many cells are converged (in every convergence category)
                        std::cout << prefix_ << "  converged cells due to" << std::endl;
                    if (epsabs_ > 0)
                        std::cout << prefix_ << "    absolute local change between rules : " << abslocal << std::endl;
                    if (epsrel_ > 0)
                        std::cout << prefix_ << "    relative local change between rules : " << rellocal << std::endl;
                    if (global_epsabs_ > 0)
                        std::cout << prefix_ << "    absolute global threshold : " << absglobal << std::endl;
                    if (global_epsrel_ > 0)
                        std::cout << prefix_ << "    relative global threshold : " << relglobal << std::endl;
                    
                    // print how many cells will need further bisection
                    double frac = bisectIDs.size() * 100. / domains.size();
                    std::cout << prefix_ << "  cells to bisect : " << bisectIDs.size() << " (";
                    if (0 < frac and frac < 1)
                        std::cout << "< 1";
                    else
                        std::cout << "~ " << (int)std::ceil(frac);
                    std::cout << "% of processed cells)" << std::endl;
                    
                    // warn if we reached the max level (no further subdivision will be done)
                    if (level == maxlevel_)
                        std::cout << prefix_ << "  maximal level " << maxlevel_ << " reached" << std::endl;
                }
                
                // move current domains to temporary array
                Domains<T,dim> oldDomains = std::move(domains);
                
                // further subdivide all not yet converged sub-domains (consecutive access)
                timer.reset(); globalEstimate = finalEstimate; pos = 0; i_pos = 0;
                # pragma omp parallel for schedule (static) if (parallel_) \
                    default (shared) \
                    shared (pos,globalEstimate) \
                    firstprivate (i_pos)
                for (unsigned i = 0; i < bisectIDs.size(); i++)
                {
                    // get next domain ID
                    # pragma omp critical
                    i_pos = pos++;
                    
                    // get old domain
                    typename Domains<T,dim>::Domain olddom = oldDomains[bisectIDs[i_pos]];
                    
                    // subdivide this domain into 2^dim subcubes
                    for (dCube const & c : olddom.cube.subdivide())
                    {
                        // compute the rough estimate of integral in this new cell 'c'
                        T val = integrate_fixed(F, ruleLow, c);
                        
                        // update global estimate over to-be-yet-processed cells
                        # pragma omp critical
                        globalEstimate += val;
                        
                        // create new sub-domain
                        typename Domains<T,dim>::Domain dom;
                        dom.set = false;
                        dom.val = val;
                        dom.cube = c;
                        
                        // add cube and data
                        domains.add(dom);
                    }
                }
                
                // update cell const and evaluation count
                ncells_ += domains.size();
                neval_ += domains.size() * nodes.at(ruleLow).n;
                
                // print some more statistics
                if (verbose_)
                {
                    std::cout << prefix_ << "  bisection time : " << timer.nice_time() << std::endl;
                    std::cout << prefix_ << "  function evaluations up to now : " << neval_ << std::endl;
                }
                
                // stop if no domain needs subdivision
                if (domains.size() == 0)
                    break;
            }
            
            // sum contributions to the integral
            result_ = globalEstimate;
            
            if (verbose_)
            {
                if (domains.size() == 0)
                    std::cout << prefix_ << "  cell list successfully emptied" << std::endl;
                else
                    std::cout << prefix_ << "  integral did not converge in " << domains.size() << " cells" << std::endl;
                std::cout << prefix_ << "  final estimate : " << result_ << std::endl;
            }
            
            if (write_vtk_ and dim <= 3)
            {
                if (verbose_)
                {
                    std::cout << std::endl;
                    std::cout << prefix_ << "Writing VTK to file \"" << vtkfilename_ << "\"." << std::endl;
                }
                
                // auxiliary data structures needed for points and cells
                std::map<std::vector<char>,std::size_t> point_histories_map;
                std::vector<std::pair<std::array<std::size_t,special::pow2(dim)>,T>> cell_data;
                
                // for every cell in grid: extrude its origin to the whole cell
                for (auto data : vtk)
                {
                    // start with the root corner of the cube (use its subdivision history, prepend zero)
                    std::vector<char> history = { 0 };
                    history.insert(history.end(), data.first.begin(), data.first.end());
                    std::vector<std::vector<char>> pthistory = extrude_point_to_regular_VTK_cell<dim>(history);
                    
                    // point indices for this cell
                    std::array<std::size_t,special::pow2(dim)> pts;
                    
                    // for all points of the current cell
                    for (std::size_t i = 0; i < special::pow2(dim); i++)
                    {
                        // check whether this point exists in the database
                        if (point_histories_map.find(pthistory[i]) == point_histories_map.end())
                        {
                            // if it doesn't, add it with a new consecutive id
                            // NOTE : The following two lines need to be split, because the first std::map access
                            //        must precede the second (and '=' is not a sequence point). Otherwise a new
                            //        element could be created in the map, which we do not want to count.
                            std::size_t npoints = point_histories_map.size();
                            point_histories_map[pthistory[i]] = npoints;
                        }
                        
                        // save index of this point
                        pts[i] = point_histories_map[pthistory[i]];
                    }
                    
                    // add cell point data to database
                    cell_data.push_back(std::make_pair(pts, data.second));
                }
                
                // sort points so that their order matches their indices
                std::vector<std::vector<char>> point_histories(point_histories_map.size());
                for (auto data : point_histories_map)
                    point_histories[data.second] = data.first;
                
                // write data file header
                std::ofstream vtkf(vtkfilename_.c_str());
                vtkf << "# vtk DataFile Version 3.0" << std::endl;
                vtkf << "Sparse grid dump from the Hex package." << std::endl;
                vtkf << "ASCII" << std::endl;
                vtkf << "DATASET UNSTRUCTURED_GRID" << std::endl;
                vtkf << "POINTS " << point_histories.size() << " float" << std::endl;
                
                // decode subdivision histories into final coordinates (= translate the binary number)
                for (auto ptsubs : point_histories)
                {
                    // compute coordinates of this point given its subdivision history
                    std::vector<double> coords(dim);
                    for (unsigned sub = 0; sub < ptsubs.size(); sub++)
                    for (int axis = 0; axis < dim; axis++)
                    if ((ptsubs[sub] & special::pow2(axis)) != 0)
                        coords[axis] += std::pow(0.5,sub);
                    
                    // write exactly three coordinates, even in the lower-dimensional case (VTK wants it this way)
                    for (unsigned axis = 0; axis < dim; axis++)
                        vtkf << coords[axis] << " ";
                    for (int axis = dim; axis < 3; axis++)
                        vtkf << 0 << " ";
                    
                    vtkf << std::endl;
                }
                
                // write cell information header
                vtkf << "CELLS " << cell_data.size() << " " << cell_data.size() * (special::pow2(dim) + 1) << std::endl;
                
                // write cell vertex list in the correct order
                for (auto cell : cell_data)
                {
                    vtkf << special::pow2(dim) << " ";
                    for (unsigned ptid : cell.first)
                        vtkf << ptid << " ";
                    vtkf << std::endl;
                }
                
                // write cell types and field header
                vtkf << "CELL_TYPES " << cell_data.size() << std::endl;
                for (std::size_t i = 0; i < cell_data.size(); i++)
                    vtkf << vtk_cell[dim] << std::endl;
                vtkf << "CELL_DATA " << cell_data.size() << std::endl;
                vtkf << "FIELD EvaluatedIntegrand " << typeinfo<T>::ncmpt << std::endl;
                
                // for all components of the data type 'T' (i.e. one for real number or two for a complex number)
                for (std::size_t i = 0; i < typeinfo<T>::ncmpt; i++)
                {
                    // write header for this component
                    vtkf << "Component" << i << " 1 " << cell_data.size() << " float" << std::endl;
                    
                    // write cell data
                    for (auto celldata : cell_data)
                        vtkf << typeinfo<T>::cmpt(i,celldata.second) << std::endl;
                }
                
                // close the file, we are done
                vtkf.close();
            }
            
            // terminate
            return ok_ = true;
        }

}; // class SparseGrid

} // namespace spgrid

#endif // HEX_SPARSE_GRID
