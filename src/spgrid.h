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

static std::vector<unsigned> vtk_cell = {
     1, // 0-D: vertex
     3, // 1-D: line
     9, // 2-D: quad
    12  // 3-D: hexahedron
};

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

template <int dim> std::vector<std::vector<std::bitset<dim>>> extrude_point_to_regular_VTK_cell (std::vector<std::bitset<dim>> const & pt)
{
    throw exception ("Regular extrusion of point into more than 3 dimensions is not implemented. Please do not use VTK grid export here.");
}

template<> inline std::vector<std::vector<std::bitset<1>>> extrude_point_to_regular_VTK_cell<1> (std::vector<std::bitset<1>> const & pt)
{
    // order required by VTK: (xmin) (xmax)
    std::vector<std::vector<std::bitset<1>>> cell_points = { pt };
    cell_points.push_back(shift_fwd<1>(0,pt));
    return cell_points;
}

template<> inline std::vector<std::vector<std::bitset<2>>> extrude_point_to_regular_VTK_cell<2> (std::vector<std::bitset<2>> const & pt)
{
    // order required by VTK: (xmin,ymin) (xmax,ymin) (xmax,ymax) (xmin,ymax)
    std::vector<std::vector<std::bitset<2>>> cell_points = { pt };
    cell_points.push_back(shift_fwd<2>(0,pt));
    cell_points.push_back(shift_fwd<2>(1,cell_points.back()));
    cell_points.push_back(shift_bck<2>(0,cell_points.back()));
    return cell_points;
}

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
    
    public:
        
        /// Constructor.
        SparseGrid ()
            : epsabs_(1e-8), epsrel_(1e-8), global_epsabs_(1e-8),
              global_epsrel_(1e-8), minlevel_(0), maxlevel_(5),
              ok_(true), verbose_(false), prefix_(""), write_vtk_(false),
              vtkfilename_("spgrid.vtk")
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
        
        /// Set absolute tolerance.
        void setLocEpsAbs (double epsabs) { epsabs_ = epsabs; }
        
        /// Set relative tolerance.
        void setLocEpsRel (double epsrel) { epsrel_ = epsrel; }
        
        /// Set global absolute tolerance.
        void setGlobEpsAbs (double epsabs) { global_epsabs_ = epsabs; }
        
        /// Set global relative tolerance.
        void setGlobEpsRel (double epsrel) { global_epsrel_ = epsrel; }
        
        /// Set minimal subdivision level.
        void setMinLevel (int minlevel) { minlevel_ = minlevel; }
        
        /// Set maximal subdivision level.
        void setMaxLevel (int maxlevel) { maxlevel_ = maxlevel; }
        
        /// Set verbosity level.
        void setVerbose (bool verbose) { verbose_ = verbose; }
        
        /// Set verbose output prefix.
        void setPrefix (std::string prefix) { prefix_ = std::string(prefix); }
        
        /// Set whether to write final grid to VTK (can be used only when dim <= 3).
        void setWriteVTK (bool write_vtk, std::string vtkfilename = "spgrid.vtk")
        {
            write_vtk_ = write_vtk;
            vtkfilename_ = vtkfilename;
        }
        
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
         * @param domain N-dimensional cube (integration volume).
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
            
            // short form for absolute value type of the template data type 'T'
            typedef decltype(std::abs(T(0))) AbsT;
            
            /**
             * @brief Integration domain info.
             * 
             * The structure 'sDomain' contains integration domain data, namely
             * the subdivision level of the domain (i.e. how many parents it has),
             * the position and extent of the domain (via the @ref ndCube type) and
             * index of the parent domain. The variables 'set' and 'val' are used
             * by the integrator to keep the values of the integral estimates.
             */
            typedef struct sDomain
            {
                int level;
                dCube cube;
                int parent;
                bool set;
                T val;
            } tDomain;
            
            const int NO_PARENT = -1;
            
            // subdomains of the integration domain (hypercube) that are yet to be processed
            std::vector<tDomain> domains = { { 0, root, NO_PARENT, false, T(0) } };
            
            // initialize subdivision count
            ncells_ = 1;
            
            // initialize evaluations count
            neval_ = 0;
            
            // there is always something to do
            while (true)
            {
                // get last added domain info
                tDomain & Dom = domains.back();
                
                // check if we already have an estimate in this subdomain
                if (Dom.set)
                {
                    // check if this is the root domain
                    if (Dom.parent == NO_PARENT)
                    {
                        ok_ = true;
                        result_ = Dom.val;
                        return ok_;
                    }
                    
                    // otherwise add this estimate to the parent estimate
                    domains[Dom.parent].set = true;
                    domains[Dom.parent].val += Dom.val;
                    
                    // drop this subdomain
                    domains.pop_back();
                    
                    // loop next
                    continue;
                }
                
                // evaluate estimates, check convergence and possibly subdivide
                {
                    // integral estimates (rough and fine)
                    T estimateLow = integrate_fixed(F, Dom.cube, ruleLow);
                    T estimateHigh = integrate_fixed(F, Dom.cube, ruleHigh);
                    
                    // update evaluation count
                    neval_ += nodes.at(ruleLow).n;
                    neval_ += nodes.at(ruleHigh).n;
                    
                    // compute difference
                    AbsT absdelta = std::abs(estimateLow - estimateHigh);
                    AbsT abshigh = std::abs(estimateHigh);
                    
                    // do not further subdivide this cell if
                    // - the absolute tolerance has been reached
                    // - the relative tolerance has been reached
                    if (Dom.level >= minlevel_ and (absdelta < epsabs_ or absdelta < epsrel_ * abshigh))
                    {
                        //std::cout << "    converged" << std::endl;
                        
                        // store the estimate
                        domains.back().set = true;
                        domains.back().val = estimateHigh;
                        
                        // loop next
                        continue;
                    }
                    // - the maximal subdivision level has been reached
                    if (Dom.level == maxlevel_)
                    {
                        // no further subdivision allowed; copy the rough estimate to parent
                        domains[Dom.parent].set = true;
                        domains[Dom.parent].val += estimateHigh;
                        
                        // drop this subdomain
                        domains.pop_back();
                        
                        // loop next
                        continue;
                    }
                    
                    // get index and level of this domain
                    int idx = domains.size() - 1;
                    int lvl = Dom.level;
                    
                    // subdivide this domain
                    for (dCube c : Dom.cube.subdivide())
                    {
                        // add subdomain
                        domains.push_back({ lvl + 1, c, idx, false, T(0) });
                        
                        // update subdivision count
                        ncells_++;
                    }
                }
            }
        }
        
        template <int dim, class Functor> bool integrate_adapt_v2
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
            
            // VTK debug data
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
                
                // debug info
                if (verbose_)
                {
                    std::cout << std::endl;
                    std::cout << prefix_ << "Level " << level << std::endl;
                    std::cout << prefix_ << "  initial estimate : " << globalEstimate << std::endl;
                    std::cout << prefix_ << "  cells to process : " << domains.size() << " (";
                    
                    double frac = domains.size() * 100. / special::pow2(level*dim);
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
                
                // check convergence for all not yet converged sub-domains
                if (level >= minlevel_) for (tDomain & dom : domains)
                {
                    // compute the fine estimate for this domain
                    T fineEstimate = integrate_fixed(F, dom.cube, ruleHigh);
                    neval_ += nodes.at(ruleHigh).n;
                    
                    // check absolute local convergence
                    if (std::abs(dom.val - fineEstimate) < epsabs_)
                    {
                        dom.set = true;
                        abslocal++;
                    }
                    
                    // check relative local convergence
                    if (std::abs(dom.val - fineEstimate) < epsrel_ * std::abs(fineEstimate))
                    {
                        dom.set = true;
                        rellocal++;
                    }
                    
                    // check absolute global convergence
                    if (std::abs(fineEstimate) < global_epsabs_)
                    {
                        dom.set = true;
                        absglobal++;
                    }
                    
                    // check relative global convergence
                    if (std::abs(fineEstimate) < global_epsrel_ * std::abs(globalEstimate))
                    {
                        dom.set = true;
                        relglobal++;
                    }
                    
                    // store better estimate
                    dom.val = fineEstimate;
                }
                
                // debug info
                if (verbose_)
                {
                    std::cout << prefix_ << "  converged cells due to" << std::endl;
                    if (epsabs_ > 0)
                        std::cout << prefix_ << "    absolute local change between rules : " << abslocal << std::endl;
                    if (epsrel_ > 0)
                        std::cout << prefix_ << "    relative local change between rules : " << rellocal << std::endl;
                    if (global_epsabs_ > 0)
                        std::cout << prefix_ << "    absolute global threshold : " << absglobal << std::endl;
                    if (global_epsrel_ > 0)
                        std::cout << prefix_ << "    relative global threshold : " << relglobal << std::endl;
                    
                    std::size_t nsub = std::count_if
                    (
                        domains.begin(),
                        domains.end(),
                        [](tDomain const & dom) -> bool { return not dom.set; }
                    );
                    double frac = nsub * 100. / domains.size();
                    
                    std::cout << prefix_ << "  cells to bisect : " << nsub << " (";
                    if (0 < frac and frac < 1)
                        std::cout << "< 1";
                    else
                        std::cout << "~ " << (int)std::ceil(frac);
                    std::cout << "% of processed cells)" << std::endl;
                }
                
                // stop if we reached the max level
                if (level == maxlevel_)
                {
                    if (verbose_)
                        std::cout << prefix_ << "  maximal level " << maxlevel_ << " reached" << std::endl;
                    
                    break;
                }
                
                // further subdivide all not yet converged sub-domains
                std::vector<tDomain> oldDomains = std::move(domains);
                for (tDomain const & oldDom : oldDomains)
                {
                    if (oldDom.set)
                    {
                        // this domain has converged -- store its value
                        finalEstimate += oldDom.val;
                        
                        // save data for potential VTK debug output
                        if (write_vtk_ and dim <= 3)
                            vtk[oldDom.cube.history()] = F(dim,oldDom.cube.centre().data());
                    }
                    else
                    {
                        // this domain needs subdivision
                        for (dCube const & c : oldDom.cube.subdivide())
                        {
                            domains.push_back({ c, false, integrate_fixed(F, c, ruleLow) });
                            neval_ += nodes.at(ruleLow).n;
                            ncells_++;
                        }
                    }
                }
                
                if (verbose_)
                    std::cout << prefix_ << "  function evaluations up to now : " << neval_ << std::endl;
                
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
                    std::cout << prefix_ << "  integral did not converge in some cells" << std::endl;
                std::cout << prefix_ << "  final estimate : " << result_ << std::endl;
            }
            
            if (write_vtk_ and dim <= 3)
            {
                if (verbose_)
                {
                    std::cout << std::endl;
                    std::cout << prefix_ << "Writing VTK to file \"" << vtkfilename_ << "\"." << std::endl;
                }
                
                std::map<std::vector<std::bitset<dim>>,std::size_t,order_vector_of_bitsets<dim>> point_histories_map;
                std::vector<std::pair<std::array<std::size_t,special::pow2(dim)>,T>> cell_data;
                
                // for every cell in grid: extrude its origin to the whole cell
                for (auto data : vtk)
                {
                    // start with the root corner of the cube (use its subdivision history)
                    std::vector<std::vector<std::bitset<dim>>> pthistory = extrude_point_to_regular_VTK_cell<dim>(data.first);
                    
                    // point indices for this cell
                    std::array<std::size_t,special::pow2(dim)> pts;
                    
                    // check whether these points exist in database (add if they don't, add with a new consecutive id)
                    for (std::size_t i = 0; i < special::pow2(dim); i++)
                    {
                        if (point_histories_map.find(pthistory[i]) == point_histories_map.end())
                        {
                            // these two lines need to be split ('=' is not a sequence point!)
                            std::size_t npoints = point_histories_map.size();
                            point_histories_map[pthistory[i]] = npoints;
                        }
                        pts[i] = point_histories_map[pthistory[i]];
                    }
                    
                    // add cell point data to database
                    cell_data.push_back(std::make_pair(pts, data.second));
                }
                
                // reorder points so that they match their indices (flip the data structure)
                std::vector<std::vector<std::bitset<dim>>> point_histories(point_histories_map.size());
                for (auto data : point_histories_map)
                    point_histories[data.second] = data.first;
                
                // write data file header
                std::ofstream vtkf(vtkfilename_.c_str());
                vtkf << "# vtk DataFile Version 2.0" << std::endl;
                vtkf << "Sparse grid dump." << std::endl;
                vtkf << "ASCII" << std::endl;
                vtkf << "DATASET UNSTRUCTURED_GRID" << std::endl;
                vtkf << "POINTS " << point_histories.size() << " float" << std::endl;
                
                // decode subdivision histories into final coordinates
                for (auto ptsubs : point_histories)
                {
                    // compute coordinates of this point given its subdivision history
                    std::vector<double> coords = root.origin();
                    for (unsigned sub = 0; sub < ptsubs.size(); sub++)
                    for (int axis = 0; axis < dim; axis++)
                    if (ptsubs[sub][axis])
                        coords[axis] += root.edge() * std::pow(0.5,sub);
                    
                    // write exactly three coordinates
                    for (unsigned axis = 0; axis < dim; axis++)
                        vtkf << coords[axis] << " ";
                    for (int axis = dim; axis < 3; axis++)
                        vtkf << 0 << " ";
                    
                    vtkf << std::endl;
                }
                
                vtkf << "CELLS " << cell_data.size() << " " << cell_data.size() * (special::pow2(dim) + 1) << std::endl;
                
                // write cell vertex list in the correct order
                for (auto cell : cell_data)
                {
                    vtkf << special::pow2(dim) << " ";
                    for (unsigned ptid : cell.first)
                        vtkf << ptid << " ";
                    vtkf << std::endl;
                }
                
                vtkf << "CELL_TYPES " << cell_data.size() << std::endl;
                for (std::size_t i = 0; i < cell_data.size(); i++)
                    vtkf << vtk_cell[dim] << std::endl;
                vtkf << "CELL_DATA " << cell_data.size() << std::endl;
                vtkf << "FIELD EvaluatedIntegrand " << typeinfo<T>::ncmpt << std::endl;
                
                // for all components of the data type 'T'
                for (std::size_t i = 0; i < typeinfo<T>::ncmpt; i++)
                {
                    vtkf << "Component" << i << " 1 " << cell_data.size() << " float" << std::endl;
                    for (auto celldata : cell_data)
                        vtkf << typeinfo<T>::cmpt(i,celldata.second) << std::endl;
                }
                vtkf.close();
            }
            
            // terminate
            return ok_ = true;
        }

}; // class SparseGrid

} // namespace spgrid

#endif // HEX_SPARSE_GRID
