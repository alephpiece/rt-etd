#ifndef PARTITION_H_
#define PARTITION_H_

#include <stdio.h>


/// \brief   Returns the first index of a block.
/// \details first = int( r*n/p )
/// \param   ntasks Number of tasks or data points to be distributed.
/// \param   rank Index of a process (or thread).
/// \param   nprocs number of processes (or threads).
/// \return  Index of the first element of a block.
size_t getBlockFirst(size_t ntasks, int rank, int nprocs) {
    size_t quo = ntasks / nprocs;
    size_t res = ntasks % nprocs;
    size_t first = (quo * rank) + (res * rank / nprocs);   // int( r*n/p )
    return first;
}


/// \brief   Returns the last index of a block.
/// \details first = int( (r+1)*n/p ) - 1
/// \param   ntasks Number of tasks or data points to be distributed.
/// \param   rank Index of a process (or thread).
/// \param   nprocs number of processes (or threads).
/// \return  Index of the last element of a block.
size_t getBlockLast(size_t ntasks, int rank, int nprocs) {
    size_t quo = ntasks / nprocs;
    size_t res = ntasks % nprocs;
    size_t last = (quo * (rank + 1)) + (res * (rank + 1) / nprocs) - 1;   // int( (r+1)*n/p ) - 1
    return last;
}


/// \brief  Returns the block size by block distribution formula.
/// \param  ntasks Number of tasks or data points to be distributed.
/// \param  rank Index of a process (or thread).
/// \param  nprocs number of processes (or threads).
/// \return Size of the block for a rank.
size_t getBlockSize(size_t ntasks, int rank, int nprocs) {
    size_t first = getBlockFirst(ntasks, rank, nprocs);
    size_t last = getBlockLast(ntasks, rank, nprocs);
    return last - first + 1;
}


/// \brief Returns the first and last block indices for each thread.
/// \param  starts Buffer to store first element indices. [OUT]
/// \param  ends Buffer to store last element indices. [OUT]
/// \param  ntasks Number of tasks or data points to be distributed.
/// \param  nthreads number of threads.
void getDeviceBlockIndices(size_t *starts, size_t *ends, size_t ntasks, int nthreads) {

    int tid;
    for (tid = 0; tid < nthreads; tid++) {
        starts[tid] = getBlockFirst(ntasks, tid, nthreads);
        ends[tid] = getBlockLast(ntasks, tid, nthreads);
        #ifdef DEBUG
        printf("[tid %d] first element = %lu, last element = %lu\n", tid, starts[tid], ends[tid]);
        #endif
    }
}


#endif  // PARTITION_H_
