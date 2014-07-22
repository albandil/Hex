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

#include <cassert>
#include <map>
#include <vector>

#include "ndcube.h"

namespace spgrid
{

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
        
typedef struct
{
    int dim;
    int n;
    std::vector<double> x;
    std::vector<double> w;
}
SparseGridNodes;

extern const std::map<SparseGridId,SparseGridNodes> nodes;

template <class T> class SparseGrid
{
    private:
        
        T result_;
        T error_;
        std::size_t neval_;
        std::size_t ncells_;
        double epsabs_;
        double epsrel_;
        int minlevel_;
        int maxlevel_;
        bool ok_;
        bool verbose_;
    
    public:
        
        SparseGrid ()
            : epsabs_(1e-8), epsrel_(1e-8), minlevel_(0), maxlevel_(5), ok_(true), verbose_(false)
        {}
        
        T result () const { return result_; }
        T error () const { return error_; }
        std::size_t evalcount () const { return neval_; }
        std::size_t cellcount () const { return ncells_; }
        int minlevel () const { return minlevel_; }
        int maxlevel () const { return maxlevel_; }
        bool ok () const { return ok_; }
        
        void setEpsAbs (double epsabs) { epsabs_ = epsabs; }
        void setEpsRel (double epsrel) { epsrel_ = epsrel; }
        void setMinLevel (int minlevel) { minlevel_ = minlevel; }
        void setMaxLevel (int maxlevel) { maxlevel_ = maxlevel; }
        
        template <class Functor, int dim> T integrate_fixed (Functor F, ndCube<dim> const & domain, SparseGridId i) const
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

        template <int dim, class Functor> bool integrate_adapt
        (
            Functor F, ndCube<dim> const & root,
            SparseGridId ruleLow, SparseGridId ruleHigh
        )
        {
            typedef ndCube<dim> dCube;
            typedef decltype(std::abs(T(0))) AbsT;
            
            typedef struct
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
                
                /// DEBUG
                //std::cout << "Domain " << domains.size() - 1 << " (parent " << Dom.parent << ")" << std::endl;
                //std::cout << "    vertices : " << Dom.cube << std::endl;
                //std::cout << "    volume = " << Dom.cube.volume() << std::endl;
                
                // check if we already have an estimate in this subdomain
                if (Dom.set)
                {
                    //std::cout << "    final estimate = " << Dom.val << std::endl;
                    
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
                    
                    //std::cout << "    low-order estimate : " << estimateLow << std::endl;
                    //std::cout << "    high-order estimate : " << estimateHigh << std::endl;
                    //std::cout << "    absdelta : " << absdelta << std::endl;
                    //std::cout << "    reldelta : " << absdelta / abshigh << std::endl;
                    
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
                        
                        // print info
                        //std::cout << "    subdivision inhibited (parent " << Dom.parent << ")" << std::endl;
                        
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

}; // class SparseGrid

} // namespace spgrid

#endif // HEX_SPARSE_GRID
