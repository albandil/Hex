//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2015, Jakub Benda, Charles University in Prague                    //
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
            if (active_)
            {
    #ifndef _OPENMP
                // initialize MPI
                MPI_Init(argc, argv);
    #else
                // initialize MPI compatible with OpenMP
                int req_flag = MPI_THREAD_FUNNELED, prov_flag;
                MPI_Init_thread(argc, argv, req_flag, &prov_flag);
                
                // check thread support
                if (prov_flag == MPI_THREAD_SINGLE)
                {
                    std::cout << "Warning: The MPI implementation doesn't support MPI_THREAD_FUNNELED. ";
                    std::cout << "Every MPI process may thus run only on a single core." << std::endl;
                }
    #endif
                // get number of processes and ID of this process
                MPI_Comm_size(MPI_COMM_WORLD, &Nproc_);
                MPI_Comm_rank(MPI_COMM_WORLD, &iproc_);
            }
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
         * @brief Broadcast array from owner to everyone.
         */
        template <class T> void bcast (int owner, NumberArray<T> & data) const
        {
            // owner : broadcast size
            int size = data.size();
            MPI_Bcast(&size, 1, MPI_INT, owner, MPI_COMM_WORLD);
            
            // all : resize data array
            data.resize(size);
            
            // owner : broadcast data
            MPI_Bcast
            (
                data.data(),
                size * typeinfo<T>::ncmpt,
                typeinfo<T>::mpicmpttype(),
                owner,
                MPI_COMM_WORLD
            );
        }
        
        /**
         * @brief Synchronize across processes by composition.
         * 
         * Synchronize array across processes. It is assumed that i-th chunk of the array is
         * present on the (i % Nproc_)-th process. That process will be used as the broadcast root.
         * This behaviour is compatible with the member function @ref isMyWork.
         * 
         * @param array Pointer to data array to synchronize.
         * @param chunksize Size of the per-process segment.
         * @param Nchunk Total number of chunks in the array. Altogether chunksize*Nchunk elements
         *               will be synchronized. If there are some elements more, they will be left
         *               untouched (and un-broadcast).
         * 
         * It is expected that the array has length equal or greater than chunksize * Nchunk.
         */
        template <class T> void sync (T * array, std::size_t chunksize, std::size_t Nchunk) const
        {
#ifndef NO_MPI
            if (active_)
            {
                for (unsigned ichunk = 0; ichunk < Nchunk; ichunk++)
                {
                    MPI_Bcast
                    (
                        array + ichunk * chunksize,
                        chunksize * typeinfo<T>::ncmpt,
                        typeinfo<T>::mpicmpttype(),
                        ichunk % Nproc_,
                        MPI_COMM_WORLD
                    );
                }
            }
#endif
        }
        
        /**
         * @brief Sum arrays to node.
         * 
         */
        template <class T> void sum (T* array, std::size_t N, int owner = 0) const
        {
#ifndef NO_MPI
            if (active_)
            {
                NumberArray<T> tmp(N);
                MPI_Status status;
                
                // slaves will send their data to owner
                if (iproc_ != owner)
                    MPI_Send(array, typeinfo<T>::ncmpt * N, MPI_DOUBLE, owner, iproc_, MPI_COMM_WORLD);
                
                // master node will receive and sum the data
                if (iproc_ == owner)
                {
                    // for all non-master processes
                    for (int i = 0; i < Nproc_; i++)
                    {
                        // skip self
                        if (i == iproc_)
                            continue;
                        
                        // receive data from a particular node into 'tmp'
                        MPI_Recv(&tmp[0], typeinfo<T>::ncmpt * N, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
                        
                        // sum the arrays
                        ArrayView<T>(N, array) += tmp;
                    }
                }
            }
#endif
        }
        
        /**
         * @brief Synchronize across processes by summing.
         * 
         * Synchronize array across processes by summing. There are no assumptions on what
         * segment is on which node. The whole arrays are broadcasted and summed on every node.
         * 
         * @param array Pointer to data array to synchronize.
         * @param Nchunk Total number of chunks in the array. Altogether chunksize*Nchunk elements
         *               will be synchronized. If there are some elements more, they will be left
         *               untouched (and un-broadcast).
         * 
         * It is expected that the array has length equal or greater than chunksize * Nchunk.
         */
        template <class T> void syncsum (T* array, std::size_t N) const
        {
#ifndef NO_MPI
            if (active_)
            {
                // master node will receive and sum the data ...
                sum(array, N, 0);
                
                // ... and broadcast the result back to all the slaves
                MPI_Bcast(array, typeinfo<T>::ncmpt * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
