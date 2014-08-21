/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2014                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_SPARSE_GRID
#define HEX_SPARSE_GRID

#include <algorithm>
#include <array>
#include <cassert>
#include <map>
#include <numeric>
#include <vector>

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
 * @brief Ordering for 'std::bitset' type.
 * 
 * This template defines a "less than" relation for the 'std::bitset' type.
 * The relation is necessary for use of 'std::bitset' as the indexing key 
 * of tree storage ('std::map').
 * 
 * The comparison is done by transforming the bitset to an unsigned integral
 * number (if its size is sufficient) or to text string of zeros and ones
 * (if more space is necessary). It is not really astonishingly fast method,
 * but it is used only in debug output of the VTK file, which makes it
 * acceptable.
 */
template <int d> class order_bitsets
{
    public:
        
        bool operator() (std::bitset<d> const & u, std::bitset<d> const & v)
        {
            // define ordering for 'std::bitset<d>' variables
            return d < 8 * sizeof(unsigned long) ?
                    u.to_ulong() < v.to_ulong() : u.to_string() < v.to_string() ;
        }
        
        
};

/**
 * @brief Ordering for vectors of bitsets.
 * 
 * This template defined a "less than" relation for vectors of bitsets, or to
 * be exact for the type std::vector&lt;std::bitset&lt;T,dim&gt;&gt;. The
 * normal lexicographical comparison library call is used.
 * 
 * This function relies on the existence of @ref order_bitsets, which defines
 * the ordering of the 'std::bitset' type.
 */
template <int d> class order_vector_of_bitsets
{
    public:
        
        bool operator()
        (
            std::vector<std::bitset<d>> const & a,
            std::vector<std::bitset<d>> const & b
        )
        {
            // compare two vectors of bitsets
            return std::lexicographical_compare
            (
                a.begin(), a.end(),
                b.begin(), b.end(),
                order_bitsets<d>()
            );
        }
};

/**
 * @brief Get further point of a cube along given axis.
 * 
 * Given a subdivision history of a sparse grid point, return a the neares shifted point along
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
      ┌───────────────┬───────────────┐
      │               │               │
      │               │               │
      │               │               │
      │       ?       │       ?       │
      │               │               │
      │               │               │
      │               │               │
  .1  ├───────┬───────┼───────────────┤
      │       │       │               │
      │   ?   │   ?   │               │
      │       │       │               │
  .01 ├───────┼───┬───┤       ?       │
      │       │   │   │               │
  .001│   ?   ├───●───┤               │
      │       │   │   │               │
  .0  └───────┴───┴───┴───────────────┘
      .0      ↑   ↑   .1           
              ↑   ↑
              .01 .011
 </pre>
 * The x-shifted point would have coordinate x=0.100, whereas the y-shifted coordinate would be y = 0.010.
 * Note that the number of digits matters; x=0.1 and x=0.100 are the same point but the subdivision
 * history is different.
 */
template <int dim> std::vector<std::bitset<dim>> shift_fwd (unsigned axis, std::vector<std::bitset<dim>> pt)
{
    // find the latest history point where the point can be shifted forward
    for (auto bit = pt.rbegin(); bit != pt.rend(); bit++)
    {
        if ((*bit)[axis])
        {
            (*bit)[axis] = 0;  // set bit to 0 (and "carry" to the next binary "digit")
            continue;
        }
        else
        {
            (*bit)[axis] = 1;  // set bit to 1 (no "carry", done)
            return pt;
        }
    }
    
    throw exception ("Cannot shift foward.");
}

/**
 * @brief Get previous point of a cube along given axis.
 * 
 * Get the previous (back-shifted) point. See @ref shift_fwd for details on the storage.
 * The back-shifted point from the example in @ref shift_fwd would be x = 0.010 (if the x-axis
 * were chosen) or y = 0.000 (if the y-axis were chosen).
 */
template <int dim> std::vector<std::bitset<dim>> shift_bck (unsigned axis, std::vector<std::bitset<dim>> pt)
{
    // find the latest history point where the point can be shifted back
    for (auto bit = pt.rbegin(); bit != pt.rend(); bit++)
    {
        if ((*bit)[axis])
        {
            (*bit)[axis] = 0;  // set bit to 0 (no "carry", done)
            return pt;
        }
        else
        {
            (*bit)[axis] = 1;  // set bit to 1 (and "carry" to the next binary "digit")
            continue;
        }
    }
    
    throw exception ("Cannot shift backward.");
}

/**
 * @brief Blind template.
 * 
 * This template should never be called (it will throw an exception otherwise).
 * Its explicit specializations for 1-D, 2-D and 3-D are expected to be called instead.
 */
template <int dim> std::vector<std::vector<std::bitset<dim>>> extrude_point_to_regular_VTK_cell (std::vector<std::bitset<dim>> const & pt)
{
    throw exception ("Regular extrusion of point into more than 3 dimensions is not implemented. Please do not use VTK grid export here.");
}

/**
 * @brief Make line segment from origin point.
 * 
 * The function will return set of points (in the form of subdivision history, see @ref shift_fwd for details)
 * making a line, starting from given origin.
 */
template<> inline std::vector<std::vector<std::bitset<1>>> extrude_point_to_regular_VTK_cell<1> (std::vector<std::bitset<1>> const & pt)
{
    // order required by VTK: (xmin) (xmax)
    std::vector<std::vector<std::bitset<1>>> cell_points = { pt };
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
template<> inline std::vector<std::vector<std::bitset<2>>> extrude_point_to_regular_VTK_cell<2> (std::vector<std::bitset<2>> const & pt)
{
    // order required by VTK: (xmin,ymin) (xmax,ymin) (xmax,ymax) (xmin,ymax)
    std::vector<std::vector<std::bitset<2>>> cell_points = { pt };
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
template<> inline std::vector<std::vector<std::bitset<3>>> extrude_point_to_regular_VTK_cell<3> (std::vector<std::bitset<3>> const & pt)
{
    // order required by VTK: (xmin,ymin,zmin) (xmax,ymin,zmin) (xmax,ymax,zmin) (xmin,ymax,zmin)
    //                        (xmin,ymin,zmax) (xmax,ymin,zmax) (xmax,ymax,zmax) (xmin,ymax,zmax)
    std::vector<std::vector<std::bitset<3>>> cell_points = { pt };
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
 * Not all of these sparse grids are implemented in Hex. But they can be
 * downloaded from the page http://sparse-grids.de/#Nodes.
 */
typedef enum
{
    /* dim = 1 */
    d1l1n1, d1l2n3, d1l3n3, d1l4n7, d1l5n7, d1l6n7, d1l7n15, d1l8n15,
    d1l9n15, d1l10n15, d1l11n15, d1l12n15, d1l13n31, d1l14n31, d1l15n31,
    d1l16n31, d1l17n31, d1l18n31, d1l19n31, d1l20n31, d1l21n31, d1l22n31,
    d1l23n31, d1l24n31, d1l25n63,
    
    /* dim = 2 */
    d2l1n1, d2l2n5, d2l3n9, d2l4n17, d2l5n33, d2l6n33, d2l7n65, d2l8n97,
    d2l9n97, d2l10n161, d2l11n161, d2l12n161, d2l13n257, d2l14n321,
    d2l15n321, d2l16n449, d2l17n449, d2l18n449, d2l19n705, d2l20n705,
    d2l21n705, d2l22n705, d2l23n705, d2l24n705, d2l25n1025,
    
    /* dim = 3 */
    d3l1n1, d3l2n7, d3l3n19, d3l4n39, d3l5n87, d3l6n135, d3l7n207, d3l8n399,
    d3l9n495, d3l10n751, d3l11n1135, d3l12n1135, d3l13n1759, d3l14n2335,
    d3l15n2527, d3l16n3679, d3l17n4447, d3l18n4447, d3l19n6495, d3l20n8031,
    d3l21n8031, d3l22n11103, d3l23n11103, d3l24n11103, d3l25n15039,
    
    /* dim = 4 */
    d4l1n1, d4l2n9, d4l3n33, d4l4n81, d4l5n193, d4l6n385, d4l7n641, d4l8n1217,
    d4l9n1985, d4l10n2881, d4l11n4929, d4l12n6465, d4l13n8705, d4l14n13697,
    d4l15n16001, d4l16n22401, d4l17n31617, d4l18n34689, d4l19n47489,
    d4l20n63873,
    
    /* dim = 5 */
    d5l1n1, d5l2n11, d5l3n51, d5l4n151, d5l5n391, d5l6n903, d5l7n1743,
    d5l8n3343, d5l9n6223, d5l10n10063, d5l11n17103, d5l12n27343, d5l13n38303,
    d5l14n60703,
    
    /* dim = 6 */
    d6l1n1, d6l2n13, d6l3n73, d6l4n257, d6l5n737, d6l6n1889, d6l7n4161,
    d6l8n8481, d6l9n16929, d6l10n30689, d6l11n53729, d6l12n93665,

    /* dim = 7 */
    d7l1n1, d7l2n15, d7l3n99, d7l4n407, d7l5n1303, d7l6n3655, d7l7n8975,
    d7l8n19855, d7l9n42031, d7l10n83247,
    
    /* dim = 8 */
    d8l1n1, d8l2n17, d8l3n129, d8l4n609, d8l5n2177, d8l6n6657, d8l7n17921,
    d8l8n43137, d8l9n97153,
    
    /* dim = 9 */
    d9l1n1, d9l2n19, d9l3n163, d9l4n871, d9l5n3463, d9l6n11527, d9l7n33679,
    d9l8n87823,
    
    /* dim = 10 */
    d10l1n1, d10l2n21, d10l3n201, d10l4n1201, d10l5n5281, d10l6n19105,
    d10l7n60225,
    
    /* dim = 11 */
    d11l1n1, d11l2n23, d11l3n243, d11l4n1607, d11l5n7767, d11l6n30471,
    
    /* dim = 12 */
    d12l1n1, d12l2n25, d12l3n289, d12l4n2097, d12l5n11073, d12l6n46977,
    
    /* dim = 13 */
    d13l1n1, d13l2n27, d13l3n339, d13l4n2679, d13l5n15367, d13l6n70279,
    
    /* dim = 14 */
    d14l1n1, d14l2n29, d14l3n393, d14l4n3361, d14l5n20833, d14l6n102369,
    
    /* dim = 15 */
    d15l1n1, d15l2n31, d15l3n451, d15l4n4151, d15l5n27671,
    
    /* dim = 16 */
    d16l1n1, d16l2n33, d16l3n513, d16l4n5057, d16l5n36097,
    
    /* dim = 17 */
    d17l1n1, d17l2n35, d17l3n579, d17l4n6087, d17l5n46343,
    
    /* dim = 18 */
    d18l1n1, d18l2n37, d18l3n649, d18l4n7249, d18l5n58657,
    
    /* dim = 19 */
    d19l1n1, d19l2n39, d19l3n723, d19l4n8551, d19l5n73303,
    
    /* dim = 20 */
    d20l1n1, d20l2n41, d20l3n801, d20l4n10001, d20l5n90561
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
         * @param F Function to integrate.
         * @param domain N-dimensional cube (integration volume).
         * @param i Sparse grid type.
         */
        template <class Functor, int dim> T integrate_fixed (Functor F, geom::ndCube<dim> const & domain, SparseGridId i) const
        {
            // check that the chosen sparse grid is available
            assert(nodes.find(i) != nodes.end());
            
            // check dimedions
            assert(dim == nodes.at(i).dim);
            
            // shorthands
            int N = nodes.at(i).n;
            std::vector<double>::const_iterator px = nodes.at(i).x.begin();
            std::vector<double>::const_iterator pw = nodes.at(i).w.begin();
            
            // compute fixed-point quadrature
            T res = 0;
            for (int j = 0; j < N; j++)
            {
                // compute translated and scaled node coordinates
                double x[dim];
                for (int k = 0; k < dim; k++)
                    x[k] = domain.origin()[k] + (*px++) * domain.edge();
                
                // evaluate function at the node
                res += (*pw++) * F(dim,x);
            }
            
            // return the result
            return res * domain.volume();
        };
        
        /**
         * @brief Adaptive sparse-grid quadrature.
         * 
         * This function will integrate the function 'F' a 'dim'-dimensional cube.
         * Two sparse grids have to be specified. For every integration sub-cell
         * the integral of the function 'F' is estimated by these two quadrature
         * rules and whenever the change is within the tolerances, the result of the
         * second grid is used as the final estimate. Otherwise when the evaluations
         * differ too much, the volume is further subdivided (up to the limit).
         * The subdivision bisects simultaneously in every dimension making always
         * 2^dim new sub-cells from every original cell.
         * 
         * Function 'F' is expected to accept two parameters:
         *  - 'n': the integer number of dimensions
         *  - 'ptr': pointer to a field of doubles of length 'n'
         * The return value can be of arbitrary arithmetic type (but consistent with
         * the template parameter of the class).
         * 
         * This function uses the function @ref integrate_fixed to evaluate the
         * integral estimates on the two specified grids.
         * 
         * @param F Function to integrate.
         * @param root N-dimensional cube (integration volume).
         * @param ruleLow Low-order rule for error estimation.
         * @param ruleHigh High-order rule for error estimation.
         */
        template <int dim, class Functor> bool integrate_adapt
        (
            Functor F, geom::ndCube<dim> const & root,
            SparseGridId ruleLow, SparseGridId ruleHigh
        )
        {
            // synonym for 'dim'-dimensional cube
            typedef geom::ndCube<dim> dCube;
            
            // integration domain info
            typedef struct sDomain
            {
                dCube cube;
                bool set;
                T val;
            } tDomain;
            
            // integral over processed cells
            T finalEstimate = 0;
            
            // unprocessed subdomains of the integration domain
            std::vector<tDomain> domains = {{ root, false, integrate_fixed(F, root, ruleLow) }};
            
            // VTK debug data (only used when write_vtk_ == true)
            std::map<std::vector<std::bitset<dim>>,T,order_vector_of_bitsets<dim>> vtk;
            
            // initialize evaluations count and number of cells
            neval_ = nodes.at(ruleLow).n;
            ncells_ = 1;
            
            // for all subdivision levels
            for (int level = 0; level <= maxlevel_; level++)
            {
                // first of all, sum the current estimates for use in the global convergence check
                T globalEstimate = finalEstimate;
                for (tDomain const & dom : domains)
                    globalEstimate += dom.val;
                
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
                
                // debug variables
                std::size_t abslocal = 0, rellocal = 0, absglobal = 0, relglobal = 0;
                std::vector<std::size_t> bisectIDs;
                Timer timer;
                
                // check convergence for all not yet converged sub-domains
                if (level >= minlevel_)
                {
                    # pragma omp parallel if (parallel_) \
                        reduction (+:abslocal,absglobal,rellocal,relglobal)
                    {
                        // per-thread data
                        std::size_t thread_neval = 0;
                        T thread_finalEstimate = 0;
                        std::vector<std::size_t> thread_bisectIDs;
                        
                        // for all domains
                        # pragma omp for schedule (static)
                        for (auto dom = domains.begin(); dom < domains.end(); dom++)
                        {
                            // compute the fine estimate for this domain
                            T fineEstimate = integrate_fixed(F, dom->cube, ruleHigh);
                            thread_neval += nodes.at(ruleHigh).n;
                            
                            // check absolute local convergence
                            if (std::abs(dom->val - fineEstimate) < epsabs_)
                            {
                                dom->set = true;
                                abslocal++;
                            }
                            
                            // check relative local convergence
                            if (std::abs(dom->val - fineEstimate) / std::abs(fineEstimate) < epsrel_)
                            {
                                dom->set = true;
                                rellocal++;
                            }
                            
                            // check absolute global (extrapolated) convergence
                            if (std::abs(dom->val - fineEstimate) * (levelVolume / dom->cube.volume()) < global_epsabs_)
                            {
                                dom->set = true;
                                absglobal++;
                            }
                            
                            // check relative global (extrapolated) convergence
                            if (std::abs(dom->val - fineEstimate) * (levelVolume / dom->cube.volume()) < global_epsrel_ * std::abs(globalEstimate))
                            {
                                dom->set = true;
                                relglobal++;
                            }
                            
                            // store the better estimate
                            dom->val = fineEstimate;
                            
                            // is the cell converged or was max level reached?
                            if (dom->set or level == maxlevel_)
                            {
                                // YES : use the current value as final
                                thread_finalEstimate += dom->val;
                                
                                // YES : save data for optional VTK debug output
                                if (write_vtk_ and dim <= 3)
                                {
                                    # pragma omp critical
                                    vtk[dom->cube.history()] = dom->val / dom->cube.volume();
                                }
                            }
                            else
                            {
                                // NO : add domain ID to the list of domains needing bisection
                                thread_bisectIDs.push_back(dom - domains.begin());
                            }
                        }
                        
                        # pragma omp critical
                        {
                            // update shared variables by contributions from this thread
                            neval_ += thread_neval;
                            finalEstimate += thread_finalEstimate;
                            bisectIDs.insert(bisectIDs.end(), thread_bisectIDs.begin(), thread_bisectIDs.end());
                        }
                    }
                    
                    if (verbose_)
                        std::cout << prefix_ << "  processing time : " << timer.nice_time() << std::endl;
                }
                
                // debug info
                if (verbose_)
                {
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
                std::vector<tDomain> oldDomains = std::move(domains);
                
                // further subdivide all not yet converged sub-domains
                timer.reset();
                # pragma omp parallel if (parallel_)
                {
                    // per-thread data
                    std::size_t Nlow = nodes.at(ruleLow).n;
                    std::size_t thread_neval = 0, thread_ncells = 0;
                    std::vector<tDomain> thread_domains;
                    
                    // for all domains that need subdivision
                    # pragma omp for schedule (static)
                    for (auto itID = bisectIDs.begin(); itID < bisectIDs.end(); itID++)
                    {
                        // subdivide this domain
                        for (dCube const & c : oldDomains[*itID].cube.subdivide())
                        {
                            // compute also the rough estimate of integral in this new cell 'c'
                            thread_domains.push_back({ c, false, integrate_fixed(F, c, ruleLow) });
                            
                            // update (per-thread) evaluation count and cell count
                            thread_neval += Nlow;
                            thread_ncells++;
                        }
                    }
                    
                    # pragma omp critical
                    {
                        // update shared variables by contributions from this thread
                        ncells_ += thread_ncells;
                        neval_ += thread_neval;
                        domains.insert(domains.end(), thread_domains.begin(), thread_domains.end());
                    }
                }
                
                if (verbose_)
                {
                    // print some more statistics
                    std::cout << prefix_ << "  bisection time : " << timer.nice_time() << std::endl;
                    std::cout << prefix_ << "  function evaluations up to now : " << neval_ << std::endl;
                }
                
                // stop if no domain needs subdivision
                if (domains.empty())
                    break;
            }
            
            // sum contributions to the integral
            result_ = finalEstimate;
            for (tDomain const & dom : domains)
                result_ += dom.val;
            
            if (verbose_)
            {
                if (domains.empty())
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
                std::map<std::vector<std::bitset<dim>>,std::size_t,order_vector_of_bitsets<dim>> point_histories_map;
                std::vector<std::pair<std::array<std::size_t,special::pow2(dim)>,T>> cell_data;
                
                // for every cell in grid: extrude its origin to the whole cell
                for (auto data : vtk)
                {
                    // start with the root corner of the cube (use its subdivision history)
                    std::vector<std::vector<std::bitset<dim>>> pthistory = extrude_point_to_regular_VTK_cell<dim>(data.first);
                    
                    // point indices for this cell
                    std::array<std::size_t,special::pow2(dim)> pts;
                    
                    // for all points of the current cell
                    for (std::size_t i = 0; i < special::pow2(dim); i++)
                    {
                        // check whether this points exists in the database
                        if (point_histories_map.find(pthistory[i]) == point_histories_map.end())
                        {
                            // if it doesn't, add it with a new consecutive id
                            // NOTE : The followig two lines need to be split, because the first std::map access
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
                std::vector<std::vector<std::bitset<dim>>> point_histories(point_histories_map.size());
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
                    std::vector<double> coords = root.origin();
                    for (unsigned sub = 0; sub < ptsubs.size(); sub++)
                    for (int axis = 0; axis < dim; axis++)
                    if (ptsubs[sub][axis])
                        coords[axis] += root.edge() * std::pow(0.5,sub);
                    
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
