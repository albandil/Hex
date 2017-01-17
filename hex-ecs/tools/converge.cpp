//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2016, Jakub Benda, Charles University in Prague                    //
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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <sys/stat.h>

double tolerance = 5e-3;
double cs_threshold = 1e-10;

std::vector<std::string> hex_ecs_args;

template <class T> T read_param
(
    std::map<std::string,std::vector<std::string>> const & data,
    std::string keyword
)
{
    if (data.find(keyword) == data.end())
    {
        return T(0);
    }
    
    if (data.at(keyword).empty())
    {
        return T(0);
    }
    
    if (data.at(keyword).size() > 1)
    {
        std::cout << "Warning: Keyword \"" << keyword << "\" should define just one value (given " << data.at(keyword).size() << ")." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    T result;
    
    std::istringstream iss (data.at(keyword).front());
    iss >> result;
    
    return result;
}

template <class T> std::vector<T> read_param_vector
(
    std::map<std::string,std::vector<std::string>> const & data,
    std::string keyword
)
{
    if (data.find(keyword) == data.end())
    {
        return std::vector<T>();
    }
    
    if (data.at(keyword).empty())
    {
        return std::vector<T>();
    }
    
    std::vector<T> result;
    T x;
    
    for (std::string const & s : data.at(keyword))
    {
        std::istringstream iss (s);
        iss >> x;
        result.push_back(x);
    }
    
    return result;
}

typedef struct
{
    double E;
    std::vector<double> pcs;
    
    double Ra, R0, Rmax;
    int L, Pi, nL, limit;
    int ni, li;
    double h;
    
    bool adjust_Ra, adjust_R0, adjust_Rmax, adjust_nL, adjust_limit, exchange;
}
calcdata;

std::vector<double> calculate (calcdata & c)
{
    //  clear previous cross section data
    c.pcs.clear();
    
    // get lengths of the individual components of the radial grid
    double igrid = c.Ra;
    double ogrid = c.R0 - c.Ra;
    double cgrid = c.Rmax - c.R0;
    
    // compose the directory name and create the directory for the calculation
    std::ostringstream oss;
    oss << "E" << c.E << "-" << "L" << c.L << "-Pi" << c.Pi << "-G" << igrid << "+" << ogrid << "+" << cgrid << "-nL" << c.nL;
    if (c.adjust_limit)
        oss << "-limit" << c.limit;
    mkdir(oss.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    std::cout << "\tcalculation " << oss.str() << std::endl;
    
    // compose names of the cross section files
    std::ostringstream oss0, oss1;
    oss0 << oss.str() << "/ics-L" << c.L << "-S0-Pi" << c.Pi << ".dat";
    oss1 << oss.str() << "/ics-L" << c.L << "-S1-Pi" << c.Pi << ".dat";
    
    // maximal principal quantum number
    int max_n = (c.E > 0 ? 5 : 1.0 / std::sqrt(-c.E));
    
    std::ifstream singlet (oss0.str()), triplet (oss1.str());
    if (not singlet.is_open() or not triplet.is_open())
    {
        singlet.close();
        triplet.close();
        
        // compose the input file
        std::ofstream ecsinp (oss.str() + "/ecs.inp");
        ecsinp << "# B-spline order.\n";
        ecsinp << "   4\n";
        ecsinp << "\n";
        ecsinp << "# ECS rotation angle in radians.\n";
        ecsinp << "   0.63\n";
        ecsinp << "\n";
        ecsinp << "# B-spline knots.\n";
        ecsinp << "# a) Real knots of the basis that is common to atomic and projectile electron.\n";
        ecsinp << "  L  0.0  0.0   4\n";
        ecsinp << "  G  0.1 10.0  0.1  1.02\n";
        ecsinp << "  L   " << 10 + c.h << " " << (unsigned long long)igrid << "  " << (unsigned long long)((igrid - 10)/c.h) << "\n";
        ecsinp << " -1\n";
        ecsinp << "# b) Real knots that are exclusive to the projectile, if any. (Start from zero.)\n";
        if (ogrid != 0)
        ecsinp << "  L    0  " << (unsigned long long)ogrid << "  " << (unsigned long long)((ogrid + 1)/c.h) << "\n";
        ecsinp << " -1\n";
        ecsinp << "# c) Complex region knots. (Start from zero.)\n";
        ecsinp << "  G    0  " << cgrid << "   " << c.h << "  1.02\n";
        ecsinp << " -1\n";
        ecsinp << "\n";
        ecsinp << "# Initial atomic states (ni, li, mi).\n";
        ecsinp << "  " << c.ni << " -1\n";
        ecsinp << "  " << c.li << "\n";
        ecsinp << "  " << "*" << "\n";
        ecsinp << "\n";
        ecsinp << "# Final atomic states (nf, lf).\n";
        for (int n = 1; n <= max_n; n++)
            ecsinp << "  " << n;
        ecsinp << "  -1\n";
        for (int n = 1; n <= max_n; n++)
            ecsinp << "  *";
        ecsinp << "\n";
        ecsinp << "\n";
        ecsinp << "# Angular momenta.\n";
        ecsinp << "# L  S  Pi nL limit exchange\n";
        ecsinp << "  " << c.L << "  " << (c.exchange ? "*" : "0") << "  " << c.Pi << "  " << c.nL << " " << c.limit << " " << c.exchange << "\n";
        ecsinp << "\n";
        ecsinp << "# Projectile charge.\n";
        ecsinp << "  -1\n";
        ecsinp << "\n";
        ecsinp << "# Atom + projectile total energies in Rydbergs.\n";
        ecsinp << "  E   " << c.E << "   -1\n";
        ecsinp << " -1\n";
        ecsinp << "\n";
        ecsinp << "# Weak magnetic field in atomic units.\n";
        ecsinp << " 0\n";
        ecsinp.close();
        
        // execute the solver
        std::ostringstream cmd;
        cmd << "( cd " << oss.str() << " ; ";
        for (std::string const & s : hex_ecs_args)
            cmd << s << " ";
//         if (igrid <= 200)
//             cmd << "--lu umfpack "; // override
        cmd << "2>&1 > ecs.log )";
        if (std::system(cmd.str().c_str()) != 0)
        {
            std::cout << "Execution failed." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        // clean the calculation directory
        std::ostringstream cl;
        cl << "( cd " << oss.str() << "; rm -f *ooc* *.bin rad-* kpa-* )";
        std::system(cl.str().c_str());
        
        // reopen the cross section files
        singlet.open(oss0.str());
        triplet.open(oss1.str());
    }
    
    // get last line from the cross section files
    std::string pcs_singlet, pcs_triplet;
    std::string s;
    while (std::getline(singlet, s)) pcs_singlet = s;
    while (std::getline(triplet, s)) pcs_triplet = s;
    
    double x;
    
    // read singlet cross sections
    std::istringstream iss_singlet (pcs_singlet);
    iss_singlet >> x; // drop energy
    while (iss_singlet >> x) c.pcs.push_back(x);
    
    // read triplet cross sections
    std::istringstream iss_triplet (pcs_triplet);
    iss_triplet >> x; // drop energy
    while (iss_triplet >> x) c.pcs.push_back(x);
    
    std::ostringstream cmd;
    cmd << "cp " << oss.str() << "/tmat-L" << c.L << "-Pi" << c.Pi << ".sql E" << c.E << "-L" << c.L << "-Pi" << c.Pi << ".sql";
    std::system(cmd.str().c_str());
    
    return c.pcs;
}

double cs_difference (std::vector<double> const & A, std::vector<double> const & B)
{
    std::cout << "\t\told cs:";
    for (std::size_t i = 0; i < A.size(); i++)
    {
        std::cout << "\t" << A[i];
    }
    std::cout << std::endl;
    
    std::cout << "\t\tnew cs:";
    for (std::size_t i = 0; i < B.size(); i++)
    {
        std::cout << "\t" << B[i];
    }
    std::cout << std::endl;
    
    double max_rel_diff = 0;
    
    double suma = std::accumulate(B.begin(), B.end(), 0.0);
    
    std::cout << "\t\tdelta:";
    for (std::size_t i = 0; i < A.size(); i++)
    {
        double rel_diff = std::abs(B[i] - A[i]) / B[i];
        std::cout << "\t" << rel_diff;
        if (std::abs(B[i]) > cs_threshold * suma)
            max_rel_diff = std::max(max_rel_diff, rel_diff);
    }
    std::cout << std::endl;
    
    return max_rel_diff;
}

void converge_energy (calcdata & c)
{
    // previous-step partial cross sections
    std::vector<double> pcs;
    
    std::cout << "E = " << c.E << std::endl;
    
    // run initial calculation
    pcs = calculate(c);
    
    if (c.adjust_Rmax)
    {
        do
        {
            // double the length of the complex grid
            pcs = c.pcs;
            c.Rmax = c.R0 + 2 * (c.Rmax - c.R0);
            calculate(c);
        }
        while (cs_difference(pcs, c.pcs) > tolerance);
    }
    
    std::cout << std::endl;
    
    if (c.adjust_Ra)
    {
        double step = c.Ra;
        
        do
        {
            // increase the length of the inner grid
            pcs = c.pcs;
            c.Ra += step;
	    if (c.R0 < c.Ra)
	    {
	        c.Rmax = c.Ra + c.Rmax - c.R0;
		c.R0 = c.Ra;
	    }
            calculate(c);
        }
        while (cs_difference(pcs, c.pcs) > tolerance);
    }
    
    std::cout << std::endl;
    
    if (c.adjust_R0)
    {
        do
        {
            // double the length of the real grid keeping the inner and complex grids as before
            pcs = c.pcs;
            c.Rmax += c.R0;
            c.R0 += c.R0;
            calculate(c);
        }
        while (cs_difference(pcs, c.pcs) > tolerance);
    }
    
    std::cout << std::endl;
    
    if (c.adjust_limit)
    {
        do
        {
            // raise the angular limit by one
            pcs = c.pcs;
            c.limit++;
            calculate(c);
        }
        while (cs_difference(pcs, c.pcs) > tolerance and c.limit <= c.L + c.Pi + c.nL);
    }
    
    std::cout << std::endl;
    
    if (c.adjust_nL)
    {
        do
        {
            // add another angular state
            pcs = c.pcs;
            c.nL++;
            if (c.limit >= 0) c.limit++;
            calculate(c);
        }
        while (cs_difference(pcs, c.pcs) > tolerance);
    }
}

double clamp (double x, double xmin, double xmax)
{
    return std::max(xmin, std::min(xmax, x));
}

double calc_angle (double x1, double y1, double x2, double y2, double x3, double y3)
{
    double ax = x1 - x2, ay = y1 - y2, amag = std::hypot(ax, ay);
    double bx = x3 - x2, by = y3 - y2, bmag = std::hypot(bx, by);
    double cos_theta = clamp((ax * bx + ay * by) / (amag * bmag), -1, +1);
    return std::acos(cos_theta);
}

int main (int argc, char* argv[])
{
    // get name of the setup file
    std::string settings = "converge.inp";
    if (argc > 1)
        settings = argv[1];
    
    // open the setup file
    std::ifstream setupfile (settings);
    if (not setupfile.is_open())
    {
        std::cout << "File \"" << settings << "\" not found." << std::endl;
        return EXIT_FAILURE;
    }
    
    // read all data
    std::string line;
    std::map<std::string, std::vector<std::string>> data;
    while (std::getline(setupfile, line))
    {
        // get rid of all leading white spaces
        while (not line.empty() and std::isspace(line.front()))
            line.erase(line.begin());
        
        // skip empty lines and lines introduced by hash symbol
        if (line.empty() or line.front() == '#')
            continue;
        
        // get first token
        std::istringstream iss (line);
        std::string s;
        iss >> s;
        
        // get other tokens
        std::string t;
        std::vector<std::string> v;
        while (iss >> t)
            v.push_back(t);
        
        // save the data
        data[s] = v;
    }
    
    // get all parameters
    
        calcdata c;
        std::vector<calcdata> alldata;
        
        std::cout << "Input:" << std::endl;
        
        //- knot spacing
        c.h = read_param<double>(data, "h");
        
        //- inner region radius
        c.Ra = read_param<double>(data, "Ra");
        c.adjust_Ra = read_param<bool>(data, "adjust_Ra");
        std::cout << "\tRa = " << c.Ra << " (adjust: " << c.adjust_Ra << ")" << std::endl;
        
        //- outer region radius
        c.R0 = read_param<double>(data, "R0");
        c.adjust_R0 = read_param<bool>(data, "adjust_R0");
        std::cout << "\tR0 = " << c.R0 << " (adjust: " << c.adjust_R0 << ")" << std::endl;
        
        //- complex region radius
        c.Rmax = read_param<double>(data, "Rmax");
        c.adjust_Rmax = read_param<bool>(data, "adjust_Rmax");
        std::cout << "\tRmax = " << c.Rmax << " (adjust: " << c.adjust_Rmax << ")" << std::endl;
        
        //- total angular momentum
        c.L = read_param<int>(data, "L");
        c.Pi = read_param<int>(data, "Pi");
        std::cout << "\tL = " << c.L << std::endl;
        std::cout << "\tS = *" << std::endl;
        std::cout << "\tPi = " << c.Pi << std::endl;
        
        //- angular limit
        c.nL = read_param<int>(data, "nL");
        c.adjust_nL = read_param<int>(data, "adjust_nL");
        std::cout << "\tnL = " << c.nL << " (adjust: " << c.adjust_nL << ")" << std::endl;
        c.limit = read_param<int>(data, "limit");
        c.adjust_limit = read_param<int>(data, "adjust_limit");
        std::cout << "\tlimit = " << c.limit << " (adjust: " << c.adjust_limit << ")" << std::endl;
        c.exchange = read_param<int>(data, "exchange");
        std::cout << "\texchange = " << c.exchange << std::endl;
        
//         //- initial and final state
        c.ni = read_param<int>(data, "ni");
        c.li = read_param<int>(data, "li");
//         c.mi = read_param<int>(data, "mi");
//         c.nf = read_param<int>(data, "nf");
//         c.lf = read_param<int>(data, "lf");
//         std::cout << "\tistate = (" << c.ni << "," << c.li << "," << c.mi << ")" << std::endl;
//         std::cout << "\tfstate = (" << c.nf << "," << c.lf << ",*)" << std::endl;
        
        //- cross sections (uninitialized)
        c.pcs = std::vector<double>{ -1 };
        
        //- convergence tolerance
        tolerance = read_param<double>(data, "tolerance");
        std::cout << "\ttolerance = " << tolerance << std::endl;
        
        //- threshold for cross sections considered in the convergence check
        cs_threshold = read_param<double>(data, "cs_threshold");
        std::cout << "\tcs_threshold = " << cs_threshold << std::endl;
        
        //- total energies
        for (double E : read_param_vector<double>(data, "E"))
        {
            c.E = E;
            alldata.push_back(c);
        }
        int nenergies = read_param<int>(data, "nenergies");
        double emin = read_param<double>(data, "emin");
        double emax = read_param<double>(data, "emax");
        for (int ie = 0; ie < nenergies; ie++)
        {
            c.E = (nenergies == 1 ? emin : emin + (emax - emin) * ie / (nenergies - 1));
            alldata.push_back(c);
        }
        
        //- total energies without separate convergence process
        for (double E : read_param_vector<double>(data, "E_add"))
        {
            c.E = E;
            c.adjust_Ra = c.adjust_R0 = c.adjust_Rmax = c.adjust_nL = 0;
            alldata.push_back(c);
        }
        
        //- sort total energies in ascending order
        std::sort(alldata.begin(), alldata.end(), [](calcdata const & a, calcdata const & b){ return a.E < b.E; });
        
        //- whether to automaticall refine the energy mesh
        bool refine_E = read_param<bool>(data, "refine_E");
        
        std::cout << "\tE = ";
        for (calcdata const & e : alldata)
            std::cout << e.E << ",";
        std::cout << std::endl;
        std::cout << "\trefine_E = " << refine_E << std::endl;
        
        //- energy convergence parameters
        int maxenergies = read_param<int>(data, "maxenergies");
        double minespacing = read_param<double>(data, "minespacing");
        std::cout << "\tminespacing = " << minespacing << std::endl;
        std::cout << "\tmaxenergies = " << maxenergies << std::endl;
        
        //- solver arguments
        hex_ecs_args = read_param_vector<std::string>(data, "args");
        std::cout << "\tadditional hex-ecs arguments:";
        for (std::string const & s : hex_ecs_args)
            std::cout << " " << s;
        std::cout << std::endl;
        std::cout << std::endl;
    
    // energy refinement loop
    do
    {
        // loop over all energies, but calculate only those, that have some uncalculated partial cross sections
        for (std::size_t ie = 0; ie < alldata.size(); ie++) if (std::any_of(alldata[ie].pcs.begin(), alldata[ie].pcs.end(), [](double x){ return x < 0; }))
        {
            converge_energy(alldata[ie]);
            
            // summary
            std::ofstream cs ("cs.dat");
            cs << "# E\tsigma" << std::endl;
            for (calcdata cd : alldata)
            {
                if (cd.pcs.front() >= 0)
                {
                    cs << cd.E;
                    for (double pcs : cd.pcs)
                        cs << "\t" << pcs;
                    cs << std::endl;
                }
            }        
        }
        
        // All energies are converged now, so we need to analyze where to place a new energy point
        if (refine_E and alldata.size() > 1 and alldata.size() < (unsigned)maxenergies)
        {
            std::vector<double> angles (alldata.size(), M_PI), dist_left (alldata.size(), 1e+30), dist_righ (alldata.size(), 1e+30);
            
            // for all transitions
            for (unsigned tr = 0; tr < alldata.front().pcs.size(); tr++)
            {
                // get bounds
                double cs_min = std::min_element(alldata.begin(), alldata.end(), [&](calcdata const & a, calcdata const & b){ return a.pcs[tr] < b.pcs[tr]; })->pcs[tr];
                double cs_max = std::max_element(alldata.begin(), alldata.end(), [&](calcdata const & a, calcdata const & b){ return a.pcs[tr] < b.pcs[tr]; })->pcs[tr];
                double E_min = std::min_element(alldata.begin(), alldata.end(), [](calcdata const & a, calcdata const & b){ return a.E < b.E; })->E;
                double E_max = std::max_element(alldata.begin(), alldata.end(), [](calcdata const & a, calcdata const & b){ return a.E < b.E; })->E;
                
                // calculate polyline angles at its vertices (except for the first and last, which will be pi)
                for (std::size_t i = 1; i < alldata.size() - 1; i++)
                {
                    // skip irrelevant data
                    if (alldata[i-1].pcs[tr] < cs_threshold or
                        alldata[i  ].pcs[tr] < cs_threshold or
                        alldata[i+1].pcs[tr] < cs_threshold)
                        continue;
                    
                    // calculate vertex angle
                    double new_angle = calc_angle
                    (
                        (alldata[i-1].E - E_min) / (E_max - E_min), (alldata[i-1].pcs[tr] - cs_min) / (cs_max - cs_min),
                        (alldata[i  ].E - E_min) / (E_max - E_min), (alldata[i  ].pcs[tr] - cs_min) / (cs_max - cs_min),
                        (alldata[i+1].E - E_min) / (E_max - E_min), (alldata[i+1].pcs[tr] - cs_min) / (cs_max - cs_min)
                    );
                    
                    // update sharpest angle and distances to the nearest vertices
                    if (new_angle < angles[i])
                    {
                        angles[i] = new_angle;
                        dist_left[i] = std::hypot(alldata[i-1].E - alldata[i].E, alldata[i-1].pcs[tr] - alldata[i].pcs[tr]);
                        dist_righ[i] = std::hypot(alldata[i+1].E - alldata[i].E, alldata[i+1].pcs[tr] - alldata[i].pcs[tr]);
                    }
                }
            }
            
            // sort angles by magnitude (ascending)
            std::vector<int> indices (angles.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](int i, int j){ return angles[i] < angles[j]; });
            
            // vertex with the largest angle
            for (int idx : indices)
            {
                // subdivide the longer edge
                if (alldata[idx].E - alldata[idx-1].E >= minespacing and dist_left[idx] >= dist_righ[idx])
                {
                    // no more convergence checking
                    c.adjust_Ra = c.adjust_R0 = c.adjust_Rmax = c.adjust_nL = 0;
                    
                    // get optimal convergence parameters from the neighbours
                    c.Ra = std::max(alldata[idx-1].Ra, alldata[idx].Ra);
                    c.R0 = c.Ra + std::max(alldata[idx-1].R0 - alldata[idx-1].Ra, alldata[idx].R0 - alldata[idx].Ra);
                    c.Rmax = c.R0 + std::max(alldata[idx-1].Rmax - alldata[idx-1].R0, alldata[idx].Rmax - alldata[idx].R0);
                    c.nL = std::max(alldata[idx-1].nL, alldata[idx].nL);
                    
                    // interpolate energy
                    c.E = 0.5 * (alldata[idx-1].E + alldata[idx].E);
                    alldata.insert(alldata.begin() + idx, c);
                    
                    std::cout << "Added energy " << alldata.size() << " of " << maxenergies << std::endl;
                    break;
                }
                else if (alldata[idx+1].E - alldata[idx].E >= minespacing)
                {
                    // no more convergence checking
                    c.adjust_Ra = c.adjust_R0 = c.adjust_Rmax = c.adjust_nL = 0;
                    
                    // get optimal convergence parameters from the neighbours
                    c.Ra = std::max(alldata[idx+1].Ra, alldata[idx].Ra);
                    c.R0 = c.Ra + std::max(alldata[idx+1].R0 - alldata[idx+1].Ra, alldata[idx].R0 - alldata[idx].Ra);
                    c.Rmax = c.R0 + std::max(alldata[idx+1].Rmax - alldata[idx+1].R0, alldata[idx].Rmax - alldata[idx].R0);
                    c.nL = std::max(alldata[idx+1].nL, alldata[idx].nL);
                    
                    // interpolate energy
                    c.E = 0.5 * (alldata[idx+1].E + alldata[idx].E);
                    alldata.insert(alldata.begin() + idx + 1, c);
                    
                    std::cout << "Added energy " << alldata.size() << " of " << maxenergies << std::endl;
                    break;
                }
            }
        }
    }
    
    // loop while any partial cross sections is less then zero (= uncalculated)
    while
    (
        std::any_of
        (
            alldata.begin(), alldata.end(),
            [](calcdata const & c)
            {
                return std::any_of(c.pcs.begin(), c.pcs.end(), [](double x){ return x < 0; });
            }
        )
    );
    
    // final summary
    std::cout << std::endl;
    std::cout << "Final summary" << std::endl;
    std::cout << "E\tsigma" << std::endl;
    for (calcdata cd : alldata)
    {
        std::cout << cd.E;
        for (double cs : cd.pcs)
            std::cout << "\t" << cs;
        std::cout << std::endl;
    }
    
    return EXIT_SUCCESS;
}
