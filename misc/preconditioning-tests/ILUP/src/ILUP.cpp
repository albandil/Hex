#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

// Use the following to run the animation in shell:
//
//     i=1; while [ -f lvl-$i.txt ] ; do usleep 500000 ; clear ; cat lvl-$i.txt ; i=$(( i+1 )) ; done
//

void putlevel (std::ostream & os, int lvl, char l = ' ', char r = ' ')
{
    if (lvl == -1)
        os << l << '.' << r;
    else
        os << l << lvl << r;
}

void plot (int Anz, std::vector<std::vector<int>> const & L, std::vector<std::vector<int>> const & A, int k, int i = -1, int j = -1, bool done = false)
{
    static int id = 0;
    
    id++;
    
    const char* STARTRED = "\e[1;31m";
    const char* STARTGREEN = "\e[1;32m";
    const char* STARTYELLOW = "\e[1;33m";
    const char* STARTWHITE = "\e[1;35m";
    const char* STARTORANGE = "\e[1;36m";
    const char* STOPCOLOR = "\e[0m";
    
    std::ofstream os ("lvl-" + std::to_string(id) + ".txt");
    
    for (unsigned a = 0; a < A.size(); a++)
    {
        for (unsigned b = 0; b < A[a].size(); b++)
        {
            if (b == k and a > b and a < i)
            {
                os << STARTORANGE; 
                putlevel(os, L[a][b]);
                os << STOPCOLOR;
            }
            else if (b == k and a == i)
            {
                os << STARTYELLOW;
                putlevel(os, L[a][b]);
                os << STOPCOLOR;
            }
            else
            {
                putlevel(os, L[a][b]);
            }
        }
        os << "             ";
        
        
        for (unsigned b = 0; b < A[a].size(); b++)
        {
            if (a < k or b < k) os << STARTGREEN;
            else if (a == k) os << STARTORANGE;
            else if (a == i /* or b == j */) os << STARTYELLOW;
            else os << STARTRED;
            
            if (i == -1)
            {
                // new focus
                if (a == k and b == k) putlevel(os, A[a][b], '[', ']');
                else putlevel(os, A[a][b]);
            }
            else if (not done)
            {
                // arrows
                if (a == k and b == j and A[a][b] != -1) putlevel(os, A[a][b], 'v', 'v');
                else if (a == i and b == j) putlevel(os, A[a][b], '|', '|');
                else putlevel(os,A[a][b]);
            }
            else
            {
                // update
                if (a == i and b == j) putlevel(os, A[a][b], '[', ']');
                else putlevel(os, A[a][b]);
            }
            
            os << STOPCOLOR;
        }
        
        os << std::endl << std::endl;
    }
    
    // calculate number of non-zeros in LU
    int maxlvl = 0;
    for (int i = 0; i < A.size(); i++)
    for (int j = 0; j < i; j++)
    if (L[i][j] > maxlvl)
        maxlvl = L[i][j];
    
    // get non-zeros per level
    std::vector<int> Lnz (maxlvl + 1);
    Lnz[0] = A.size(); // diagonal
    for (int i = 0; i < A.size(); i++)
    for (int j = 0; j < i; j++)
    if (L[i][j] != -1)
        Lnz[L[i][j]] += 2;
    
    // print stats
    os << std::endl;
    os << "A nonzeros: " << Anz << std::endl;
    os << "LU nonzeros: " << std::accumulate(Lnz.begin(), Lnz.end(), 0) << std::endl;
    for (unsigned lvl = 0; lvl <= maxlvl; lvl++)
        os << " - level " << lvl << ": " << Lnz[lvl] << std::endl;
}

#define TWODIM

int main (void)
{
    int n = 5;
    int N = n*n;
    
    std::vector<std::vector<int>> A (N, std::vector<int>(N,-1));
    std::vector<std::vector<int>> L (N, std::vector<int>(N,-1));
    
    srand(time(0));
    
    // initial non-zero pattern
#if defined ( TRIDIAG )
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int j = 0; j < N; j++)
    {
        switch (std::abs(i - j))
        {
            case 0: case 1: A[i][j] = 0;
        }
        
    }
#elif defined ( PENTADIAG )
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    for (int j = 0; j < N; j++)
    {
        switch (std::abs(i - j))
        {
            case 0: case 2: case 4: A[i][j] = 0;
        }
        
    }
#elif defined ( MIXED )
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    {
        switch (std::abs(i - j))
        {
            case 0: case 1: case 4: case 5: A[i][j] = 0;
        }
        
    }
#elif defined ( TWODIM )
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++)
    for (int l = 0; l < n; l++)
    {
        if (std::abs(i - k) <= 1 and std::abs(j - l) <= 1)
        {
            A[i*n+j][k*n+l] = 0;
        }
        
    }
#elif defined ( LEADFULL )
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    {
        if (i == 0 or j == 0 or std::abs(i - j) <= 1)
            A[i][j] = 0;
    }
#elif defined ( TAILFULL )
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    {
        if (i == N - 1 or j == N - 1 or std::abs(i - j) <= 1)
            A[i][j] = 0;
    }
#elif defined ( RANDOM )
    for (int i = 0; i < N; i++)
    for (int j = i; j < N; j++)
    if (i == j or rand() > 0.75 * RAND_MAX)
    {
        A[i][j] = A[j][i] = 0;
    }
#else
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    {
        A[i][j] = 0;
    }
#endif
    
    // calculate number of non-zeros in A
    int Anz = 0;
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    if (A[i][j] != -1)
        Anz++;
    
    // initial non-zero L factor
    for (int i = 0; i < N; i++)
        L[i][i] = 0;
    
    // Gauss elimination
    for (int k = 0; k < N; k++) // along diagonal
    {
        plot(Anz, L, A, k);
        
        // "subtract this row from below rows"
        for (int i = k + 1; i < N; i++) if (A[i][k] != -1)
        {
            L[i][k] = A[i][k];
            
            for (int j = k; j < N; j++) if (A[k][j] != -1)
            {
                plot(Anz, L, A, k, i, j, false);
                if (j == k)
                {
                    // erase leading element of the row
                    //A[i][j] = -1;
                }
                else if (A[i][j] == -1)
                {
                    // add next fill-in
                    A[i][j] = std::max(A[k][j],A[i][k]) + 1;
                }
                else
                {
                    // update level
                    A[i][j] = std::min(A[i][j], (A[k][j],A[i][k]) + 1);
                }
                plot(Anz, L, A, k, i, j, true);
            }
        }
    }
    
    plot(Anz, L, A, N);
    
    return 0;
}
