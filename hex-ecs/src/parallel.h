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

/**
 * @brief MPI info.
 * 
 * The class Parallel holds some useful MPI data, like the rank and size of the
 * communicator.
 */
class Parallel
{
    public:
        
        /**
         * @brief Constructor.
         * 
         * Initializes MPI (using the MPI_Init) function. Also determines
         * the rank and size of the communicator.
         * @param active Whether to turn the MPI on or keep with the computation
         *               sequential.
         */
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
        
        /// Returns true if this process is the master process (id = 0).
        inline bool IamMaster() const { return iproc_ == 0; }
        
        /// Returns true if the chunk is assigned to this process.
        inline bool isMyWork(int i) const { return i % Nproc_ == iproc_; }
        
        /// Returns true if the MPI is active.
        inline bool active() const { return active_; }
        
        /// Returns the rank of the communicator (process is).
        inline int iproc() const { return iproc_; }
        
        /// Returns the size of the communicator (process count).
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
