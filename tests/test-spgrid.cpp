#include <iostream>

#include "../src/spgrid.cpp"

// template <class T> using integrand_t = void (*) (int n, int dim, double const * origin, double range, const double * scale, T * eval);

// Result according to Maxima 5.35.1 : 2
void F0 (int n, int dim, double const * origin, double range, const double * scale, double * eval)
{
    for (int i = 0; i < n; i++)
    {
        eval[i] = 0;
        for (int k = 0; k < dim; k++)
        {
            double x = origin[k] + range * scale[i * dim + k];
            eval[i] += x*x;
        }
    }
}

// Result according to Maxima 5.35.1 : 31π⁶/30240 ≈ 0.9855510912974348
void F1 (int n, int dim, double const * origin, double range, const double * scale, double * eval)
{
    for (int i = 0; i < n; i++)
    {
        eval[i] = 1;
        for (int k = 0; k < dim; k++)
        {
            double x = origin[k] + range * scale[i * dim + k];
            eval[i] *= x;
        }
        eval[i] = 1 / (1 + eval[i]);
    }
}

// Result according to Maxima 5.35.1 : (log(-1)*log(2) + li₂(2) - π²/4)²/4 ≈ 0.1691130052673656
void F2 (int n, int dim, double const * origin, double range, const double * scale, double * eval)
{
    for (int i = 0; i < n; i++)
    {
        double x[dim];
        for (int k = 0; k < dim; k++)
            x[k] = origin[k] + range * scale[i * dim + k];
        
        eval[i] = x[0] / (1 + x[1] * x[2]) * x[3] / (1 + x[4] * x[5]);
    }
}

std::pair < spgrid::SparseGrid<double>::integrand_t, double > tests[3] = {
    { &F0, 2.0000000000000000 },
    { &F1, 0.9855510912974348 },
    { &F2, 0.1691130052673656 }
};

int main (void)
{
    spgrid::SparseGrid<double> G;
    
    for (int i = 0; i < 3; i++)
    {
        G.integrate_adapt<6>(tests[i].first, spgrid::d6l4n257, spgrid::d6l5n737);
        
        std::cout << "Test " << i + 1 << std::endl;
        std::cout << "\tresult = " << G.result() << " (exact: " << tests[i].second << ")" << std::endl;
        std::cout << "\tnevals = " << G.evalcount() << std::endl;
        std::cout << "\tncells = " << G.cellcount() << std::endl;
    }
    
    return 0;
}
