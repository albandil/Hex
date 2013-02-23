#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <limits>

const double inf = std::numeric_limits<double>::infinity();

int main (int argc, char *argv[])
{
    std::ifstream in(argv[1]);
    if (not in.good())
    {
        std::cerr << "Cannot open file \"" << argv[1] << "\"" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::vector<double> stacked_min, stacked_max;
    std::string line;
    
    while (std::getline(in, line), not in.eof())
    {
        if (line.empty())
            continue;
        
        std::istringstream ss(line);
        double min = inf, max = -inf;
        
        while (not ss.eof())
        {
            double x;
            ss >> x;
            
            min = std::min(min, x);
            max = std::max(max, x);
        }
        
        stacked_min.push_back(min);
        stacked_max.push_back(max);
    }
    
    for (size_t i = 0; i < stacked_max.size(); i++)
        std::cout << stacked_min[i] << "\t" << stacked_max[i] << std::endl;
    
    return EXIT_SUCCESS;
}

