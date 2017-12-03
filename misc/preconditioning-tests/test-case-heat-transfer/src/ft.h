#ifndef FT_H
#define FT_H

/**
 * Transform all M segments of the source vector @c x to @c y, assuming that
 * all segments have the same length N.
 */
void DFT (int sign, const cArrayView x, cArrayView y, int M, int N)
{
    assert(x.size() == M*N);
    assert(y.size() == M*N);
    
    // precompute the weights
    cArray weights (N);
    for (int i = 0; i < N; i++)
    {
        weights[i].real(std::cos(sign * i * 6.283185307179586 / N));
        weights[i].imag(std::sin(sign * i * 6.283185307179586 / N));
    }
    
    // transform
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            y[i * N + j] = 0;
            
            for (int k = 0; k < N; k++)
                y[i * N + j] += weights[(j * k) % N] * x[i * N + k];
            
            y[i * N + j] /= std::sqrt(N);
        }
    }
}

#endif
