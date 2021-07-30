#include "device_structs.h"

#include <simd.h>
#include <slave.h>
#include <dma.h>

#include <stddef.h>
#include <stdio.h>

#ifndef DEBUGTID
#define DEBUGTID 0
#endif

#define TASKSIZE 1980           // number of data points processed at one time (must be divisible by 4)

__thread_local_fix int rank;    // MPI rank
__thread_local_fix size_t block_size;    // iend - istart + 1

//__thread_local double k;
__thread_local double *fx;    // valid length = TASKSIZE
__thread_local double *x;     // valid length = TASKSIZE+2
__thread_local double *beta;  // valid length = TASKSIZE+1
__thread_local double *alpha; // valid length = TASKSIZE+1

__thread_local volatile unsigned long num_gets, num_puts;  // counters for DMA replies


extern DeviceArgs args_dev;

/// \brief Initialize and hold some constants on device.
void initializeConstantsOnDevice(DeviceConstants *constants_dev) {
    int tid = athread_get_id(-1);

    rank = constants_dev->rank;
    size_t dev_n = constants_dev->dev_n;
    block_size = dev_n / 64;

    #ifdef DEBUG
    if (tid == DEBUGTID) {
        printf("[Rank %d][tid %d] block size = %lu, range = [%lu, %lu]\n",
            rank, tid, block_size, tid * block_size, (tid + 1) * block_size - 1);
        printf("[Rank %d][tid %d] iterations = %lu, res = %lu\n",
            rank, tid, block_size / TASKSIZE, block_size % TASKSIZE);
    }
    #endif
}


/// \brief   Compute F() on device.
/// \details To compute F(), we need to allocate 4 vectors:
///             fx: length n
///             x: length n+3
///             alpha: length n+1
///             beta: length n+1
///          so that the number of data points is bounded by
///             (4n+4)*8 <= 64KB = 64*1024
///          => n <= 2047
void computeFOnDevice() {

    // Initialize variables
    double x0 = args_dev.x0;
    double k = args_dev.k;

    // Partioning
    int tid = athread_get_id(-1);
    size_t batch_size = TASKSIZE;
    size_t num_iterations = block_size / TASKSIZE;
    size_t res_points = block_size % TASKSIZE;

    // Allocation
    fx = (double *)ldm_malloc((TASKSIZE+3)*sizeof(double));
    x = (double *)ldm_malloc((TASKSIZE+3)*sizeof(double));
    alpha = (double *)ldm_malloc((TASKSIZE+3)*sizeof(double));
    beta = (double *)ldm_malloc((TASKSIZE+3)*sizeof(double));

    #ifdef DEBUG
    if (tid == DEBUGTID) {
        printf("[Rank %d][tid %d] x0 = %e, k = %f\n", rank, tid, x0, k);
    }
    #endif

    #ifdef USE_DEVICE_SIMD
    // SIMD variables
    doublev4 vfx, vx, vx1, vx2, vbeta, vbeta1, valpha1;
    doublev4 vx0 = x0;
    doublev4 vk = k;
    #endif

    size_t j;
    for (j = 0; j <= num_iterations; j++) {
        size_t istart = tid * block_size + j * TASKSIZE;  // the start index of the block of this task
        if (j == num_iterations)
            batch_size = res_points;

        #ifdef DEBUG
        if (tid == DEBUGTID) {
            printf("[Rank %d][tid %d] iteration %lu, batch size = %lu, istart = %lu\n", rank, tid, j, batch_size, istart);
        }
        #endif

        // Fetch data from the host
        num_gets = 0;
        athread_get(PE_MODE, args_dev.x + istart, &x[0], (batch_size+2)*8, (void *)&num_gets, 0, 0, 0);
        athread_get(PE_MODE, args_dev.alpha + istart, &alpha[0], (batch_size+1)*8, (void *)&num_gets, 0, 0, 0);
        athread_get(PE_MODE, args_dev.beta + istart, &beta[0], (batch_size+1)*8, (void *)&num_gets, 0, 0, 0);
        while (num_gets != 3);

        size_t i;
        #ifdef USE_DEVICE_SIMD
        // Computing (simd)
        for (i = 0; i < (batch_size/4)*4; i+=4) {
            simd_load(vx, &(x[i]));
            simd_load(vbeta, &(beta[i]));
            simd_loadu(vx1, &(x[i+1]));
            simd_loadu(vbeta1, &(beta[i+1]));
            simd_loadu(vx2, &(x[i+2]));
            simd_loadu(valpha1, &(alpha[i+1]));

            vfx = simd_vmad(valpha1, vx2, simd_vmss(vx0, vbeta * vx, vx1 * simd_vmss(vx0, vbeta1, vk)));

            simd_store(vfx, &(fx[i]));
        }

        for (i = (batch_size/4)*4; i < batch_size; i++) {
            fx[i] = x0 * (beta[i] * x[i]) - x[i+1] * (x0 * beta[i+1] - k) + alpha[i+1] * x[i+2];
        }
        #else
        // Computing
        for (i = 0; i < batch_size; i++) {
            fx[i] = x0 * (beta[i] * x[i]) - x[i+1] * (x0 * beta[i+1] - k) + alpha[i+1] * x[i+2];
        }
        #endif

        // Send data back to the host
        num_puts = 0;
        athread_put(PE_MODE, &fx[0], args_dev.fx + istart, batch_size*8, (void *)&num_puts, 0, 0);
        while (num_puts != 1);
    }

    // Deallocation
    ldm_free(fx, (TASKSIZE+3)*sizeof(double));
    ldm_free(x, (TASKSIZE+3)*sizeof(double));
    ldm_free(alpha, (TASKSIZE+3)*sizeof(double));
    ldm_free(beta, (TASKSIZE+3)*sizeof(double));
}
