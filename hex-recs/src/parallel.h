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

#ifndef HEX_PARALLEL
#define HEX_PARALLEL

#ifndef NO_MPI
    #include <mpi.h>
#endif

#include "arrays.h"

/**
 * @brief MPI info.
 * 
 * The class Parallel holds some useful MPI data, like the rank and size of the
 * communicator. This class should be used for all parallel communication using
 * MPI.
 */
class Parallel
{
    public:
        
        /**
         * @brief Constructor.
         * 
         * Initializes MPI (using the MPI_Init) function. Also determines
         * the rank and size of the communicator.
         * @param argc Argc as received by main().
         * @param argv Argv as received by main().
         * @param active Whether to turn the MPI on or keep with the computation
         *               sequential.
         */
        Parallel (int* argc, char*** argv, bool active)
            : active_(active), iproc_(0), Nproc_(1)
        {
#ifndef NO_MPI
            MPI_Init (argc, argv);
            MPI_Comm_size (MPI_COMM_WORLD, &Nproc_);
            MPI_Comm_rank (MPI_COMM_WORLD, &iproc_);
#else
            active_ = false;
#endif
        }
        
        ~Parallel ()
        {
#ifndef NO_MPI
            if (active_)
                MPI_Finalize();
#endif
        }
        
        //
        // getters
        //
        
        /**
         * @brief Returns true if this process is the master process.
         * 
         * This function returns true if this process is the master process,
         * i.e. it has ID zero.
         */
        inline bool IamMaster () const
        {
            return iproc_ == 0;
        }
        
        /**
         * @brief Returns true if the work item is assigned to this process.
         * 
         * If there are several work items, i-th one is assigned to the
         * (i % Nproc_)-th process. This function returns true whenever given
         * wirk item (indexed by "i") is to be processed by the current process.
         * 
         * @param i Work item ID.
         */
        inline bool isMyWork (int i) const
        {
            return i % Nproc_ == iproc_;
        }
        
        /**
         * @brief Returns true if the MPI is active.
         */
        inline bool active () const
        {
            return active_;
        }
        
        /**
         * @brief Returns the rank of the communicator (process is).
         */
        inline int iproc () const
        {
            return iproc_;
        }
        
        /**
         * @brief Returns the size of the communicator (process count).
         */
        inline int Nproc () const
        {
            return Nproc_;
        }
        
        /**
         * @brief Synchronize across processes.
         * 
         * Synchronize array across processes. It is assumed that i-th chunk of the array is
         * present on the (i % Nproc_)-th process. That process will be used as the broadcast root.
         * This behaviour is compatible with the member function @ref isMyWork.
         * 
         * @param array Array to synchronize,
         * @param chunksize Size of the per-process segment.
         * @param Nchunk Total number of chunks in the array. Altogether chunksize*Nchunk elements
         *               will be synchronized. If there are some elements more, they will be left
         *               untouched (and un-broadcast).
         */
        template <class T> void sync (ArrayView<T> array, std::size_t chunksize, std::size_t Nchunk) const
        {
#ifndef NO_MPI
            if (active_)
            {
                for (unsigned ichunk = 0; ichunk < Nchunk; ichunk++)
                {
                    // relevant process will broadcast this chunk's data
                    MPI_Bcast
                    (
                        array.data() + ichunk * chunksize,
                        chunksize,
                        MPI_DOUBLE_COMPLEX,
                        ichunk % Nproc_,
                        MPI_COMM_WORLD
                    );
                }
            }
#endif
        }
        
        /**
         * @brief Wait for completition of all running tasks.
         * 
         * Inserts a MPI BARRIER.
         */
        void wait () const
        {
#ifndef NO_MPI
            if (active_)
                MPI_Barrier(MPI_COMM_WORLD);
#endif
        }
        
    private:
        
        // whether the MPI is on
        bool active_;
        
        // communicator rank
        int iproc_;
        
        // communicator size
        int Nproc_;
};

#endif
