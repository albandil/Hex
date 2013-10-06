/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEX_PARALLEL
#define HEX_PARALLEL

#ifndef NO_MPI
#include <mpi.h>
#endif

class Parallel
{
    public:
        
        // constructor
        Parallel(bool active) : active_(active), iproc_(0), Nproc_(1)
        {
#ifndef NO_MPI
            MPI_Init(nullptr, nullptr);
            MPI_Comm_size(MPI_COMM_WORLD, &Nproc_);
            MPI_Comm_rank(MPI_COMM_WORLD, &iproc_);
#else
            active_ = false;
#endif
        }
        
        //
        // getters
        //
        
        inline bool IamMaster() const { return iproc_ == 0; }
        
        inline bool active() const { return active_; }
        inline int iproc() const { return iproc_; }
        inline int Nproc() const { return Nproc_; }
        
    private:
        
        // whether the MPI is on
        bool active_;
        
        // communicator rank
        int iproc_;
        
        // communicator size
        int Nproc_;
};

#endif
