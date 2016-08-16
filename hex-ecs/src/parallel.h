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

#ifndef HEX_PARALLEL
#define HEX_PARALLEL

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "hex-arrays.h"

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
        Parallel (int* argc, char*** argv, bool active, int groupsize)
            : active_(active), iproc_(0), Nproc_(1), igroup_(0), Ngroup_(1), groupsize_(groupsize)
        {
#ifdef WITH_MPI
            if (active_)
            {
    #ifndef _OPENMP
                // initialize MPI
                MPI_Init(argc, argv);
    #else
                // initialize MPI compatible with OpenMP
                int req_flag = MPI_THREAD_SERIALIZED, prov_flag;
                MPI_Init_thread(argc, argv, req_flag, &prov_flag);
                
                // check thread support
                if (prov_flag != MPI_THREAD_SERIALIZED)
                {
                    std::cout << "Warning: The MPI implementation doesn't support MPI_THREAD_SERIALIZED. ";
                    std::cout << "Every MPI process may thus run only on a single core." << std::endl;
                }
    #endif
                // get number of processes and ID of this process
                MPI_Comm_size(MPI_COMM_WORLD, &Nproc_);
                MPI_Comm_rank(MPI_COMM_WORLD, &iproc_);
                
                // check parameter compatibility
                if (Nproc_ % groupsize_ != 0)
                    HexException("Number of processes (currently %d) must be integer mutiple of groupsize (currently %d).", Nproc_, groupsize_);
                
                // calculate number of groups
                Ngroup_ = Nproc_ / groupsize_;
                
                // calculate group rank
                igroup_ = iproc_ / groupsize_;
                
                // prepare group communicators
                for (int i = 0; i < Ngroup_; i++)
                {
                    // split a new communicator for this group
                    MPI_Comm newcomm;
                    MPI_Comm_split
                    (
                        MPI_COMM_WORLD,
                        igroup_ == i ? 0 : MPI_UNDEFINED,
                        igroupproc(),
                        &newcomm
                    );
                    
                    // store this communicator
                    if (igroup_ == i)
                        groupcomm_ = newcomm;
                }
                
                // construct also set of groups' master processes
                MPI_Comm_split
                (
                    MPI_COMM_WORLD,
                    IamGroupMaster() ? 0 : MPI_UNDEFINED,
                    igroup_,
                    &mastergroup_
                );
            }
#else
            active_ = false;
#endif // WITH_MPI
        }
        
        ~Parallel ()
        {
#ifdef WITH_MPI
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
        bool IamMaster () const
        {
            return iproc_ == 0;
        }
        
        /**
         * @brief Returns true if this process is the master process within a group.
         * 
         * This function returns true if this process is the master process within its,
         * i.e. it has local ID zero.
         */
        bool IamGroupMaster () const
        {
            return iproc_ % groupsize_ == 0;
        }
        
        /**
         * @brief Returns true if the work item is assigned to this process.
         * 
         * If there are several work items, i-th one is assigned to the
         * (i % Nproc_)-th process. This function returns true whenever given
         * work item (indexed by "i") is to be processed by the current process.
         * 
         * @param i Work item ID.
         */
        bool isMyWork (int i) const
        {
            return i % Nproc_ == iproc_;
        }
        
        /**
         * @brief Returns true if the work item is assigned to this process' group.
         * 
         * If there are several work items, i-th one is assigned to the
         * (i % Ngroup)-th group. This function returns true whenever given
         * work item (indexed by "i") is to be processed by the current process' group.
         * 
         * @param i Work item ID.
         */
        bool isMyGroupWork (int i) const
        {
            return i % Ngroup_ == igroup_;
        }
        
        /// Returns true if the MPI is active.
        bool active () const { return active_; }
        
        /// Returns the rank of the communicator (process is).
        int iproc () const { return iproc_; }
        
        /// Returns the rank of process group.
        int igroup () const { return igroup_; }
        
        /// Returns index of rank within the group.
        int igroupproc () const { return iproc_ % groupsize_; }
        
        /// Returns the size of the communicator (process count).
        int Nproc () const { return Nproc_; }
        
        /// Returns the number of groups.
        int Ngroup () const { return Ngroup_; }
        
        /// Return group size.
        int groupsize () const { return groupsize_; }
        
        /**
         * @brief Broadcast array from owner to everyone.
         */
        template <class T> void bcast (int owner, NumberArray<T> & data) const
        {
#ifdef WITH_MPI
            if (active_)
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
#endif
        }
        
        /**
         * @brief Broadcast array from owner to everyone in the group.
         */
        template <class T> void bcast_g (int igroup, int groupowner, T * data, std::size_t N) const
        {
#ifdef WITH_MPI
            if (active_ and igroup == igroup_ and groupsize_ > 1)
            {
                // groupowner : broadcast data
                MPI_Bcast
                (
                    data,
                    N * typeinfo<T>::ncmpt,
                    typeinfo<T>::mpicmpttype(),
                    groupowner,
                    groupcomm_
                );
            }
#endif
        }
        
        /**
         * @brief Send data to a process.
         */
        template <class T> void send (T const * array, std::size_t size, int origin, int destination) const
        {
#ifdef WITH_MPI
            if (active_ and iproc_ == origin)
            {
                MPI_Send
                (
                    const_cast<T*>(array),
                    typeinfo<T>::ncmpt * size,
                    typeinfo<T>::mpicmpttype(),
                    destination,
                    origin,
                    MPI_COMM_WORLD
                );
            }
#endif
        }
        
        /**
         * @brief Receive data from a process.
         */
        template <class T> void recv (T * array, std::size_t size, int origin, int destination) const
        {
#ifdef WITH_MPI
            if (active_ and iproc_ == destination)
            {
                MPI_Status status;
                
                MPI_Recv
                (
                    array,
                    typeinfo<T>::ncmpt * size,
                    typeinfo<T>::mpicmpttype(),
                    origin,
                    origin,
                    MPI_COMM_WORLD,
                    &status
                );
            }
#endif
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
#ifdef WITH_MPI
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
        
        template <class T> void sync_m (T * array, std::size_t chunksize, std::size_t Nchunk) const
        {
#ifdef WITH_MPI
            if (active_ and IamGroupMaster())
            {
                for (unsigned ichunk = 0; ichunk < Nchunk; ichunk++)
                {
                    MPI_Bcast
                    (
                        array + ichunk * chunksize,
                        chunksize * typeinfo<T>::ncmpt,
                        typeinfo<T>::mpicmpttype(),
                        ichunk % Ngroup_,
                        mastergroup_
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
#ifdef WITH_MPI
            if (active_)
            {
                MPI_Reduce
                (
                    (iproc_ == owner ? MPI_IN_PLACE : array),
                    (iproc_ == owner ? array : nullptr),
                    typeinfo<T>::ncmpt * N,
                    typeinfo<T>::mpicmpttype(),
                    MPI_SUM,
                    owner,
                    MPI_COMM_WORLD
                );
            }
#endif
        }
        
        /**
         * @brief Synchronize across processes by summing.
         * 
         * Synchronize array across processes by summing.
         * 
         * @param array Pointer to data array to synchronize.
         * @param Nchunk Number of elements in the array to sum-synchronize. If there are some
         *               elements more, they will be left untouched (and un-broadcast).
         */
        template <class T> void syncsum (T* array, std::size_t N) const
        {
#ifdef WITH_MPI
            if (active_)
            {
                MPI_Allreduce
                (
                    MPI_IN_PLACE,
                    array,
                    typeinfo<T>::ncmpt * N,
                    typeinfo<T>::mpicmpttype(),
                    MPI_SUM,
                    MPI_COMM_WORLD
                );
            }
#endif
        }
        
        /**
         * @brief Synchronize across group's processes by summing.
         * 
         * Synchronize array across group's processes by summing.
         * 
         * @param array Pointer to data array to synchronize.
         * @param Nchunk Number of elements in the array to sum-synchronize. If there are some
         *               elements more, they will be left untouched (and un-broadcast).
         */
        template <class T> void syncsum_g (T* array, std::size_t N) const
        {
#ifdef WITH_MPI
            if (active_)
            {
                MPI_Allreduce
                (
                    MPI_IN_PLACE,
                    array,
                    typeinfo<T>::ncmpt * N,
                    typeinfo<T>::mpicmpttype(),
                    MPI_SUM,
                    groupcomm_
                );
            }
#endif
        }
        
        /**
         * @brief Sum array within group to a given process.
         */
        template <class T> void sum_g (T* array, std::size_t N, int destination) const
        {
#ifdef WITH_MPI
            if (active_ and groupsize_ > 1)
            {
                MPI_Reduce
                (
                    (igroupproc() == destination ? MPI_IN_PLACE : array),
                    (igroupproc() == destination ? array : nullptr),
                    typeinfo<T>::ncmpt * N,
                    typeinfo<T>::mpicmpttype(),
                    MPI_SUM,
                    destination,
                    groupcomm_
                );
            }
#endif
        }
        
        /**
         * @brief Sum array from all group masters to one of the group masters.
         */
        template <class T> void mastersum (T* array, std::size_t N, int destgroup) const
        {
#ifdef WITH_MPI
            if (active_ and IamGroupMaster() and Ngroup_ > 1)
            {
                MPI_Reduce
                (
                    (igroup_ == destgroup ? MPI_IN_PLACE : array),
                    (igroup_ == destgroup ? array : nullptr),
                    typeinfo<T>::ncmpt * N,
                    typeinfo<T>::mpicmpttype(),
                    MPI_SUM,
                    destgroup,
                    mastergroup_
                );
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
#ifdef WITH_MPI
            if (active_)
                MPI_Barrier(MPI_COMM_WORLD);
#endif
        }
        
        void wait_g () const
        {
#ifdef WITH_MPI
            if (active_)
                MPI_Barrier(groupcomm_);
#endif
        }
        
    private:
        
        // whether the MPI is on
        bool active_;
        
        // global communicator rank
        int iproc_;
        
        // global communicator size
        int Nproc_;
        
        // group index
        int igroup_;
        
        // group count
        int Ngroup_;
#ifdef WITH_MPI
        // MPI group communicator
        MPI_Comm groupcomm_;
        
        // MPI master group communicator
        MPI_Comm mastergroup_;
#endif
        // local communicator size
        int groupsize_;
};

#endif
