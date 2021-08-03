#include "device_structs.h"

#include <simd.h>
#include <slave.h>
#include <dma.h>

#include <stddef.h>
#include <stdio.h>

#ifndef DEBUGTID
#define DEBUGTID 0
#endif

#define THREAD_BATCH_SIZE 256  // number of data points processed at one time (must be divisible by 4)

__thread_local_fix int rank;    // MPI rank
__thread_local_fix size_t block_size;    // iend - istart + 1

//__thread_local double k;
__thread_local double *(fx[2]);    // valid length = THREAD_BATCH_SIZE
__thread_local double *(x[2]);     // valid length = THREAD_BATCH_SIZE+2
__thread_local double *(beta[2]);  // valid length = THREAD_BATCH_SIZE+1
__thread_local double *(alpha[2]); // valid length = THREAD_BATCH_SIZE+1

__thread_local volatile unsigned long num_gets[2], num_puts[2];  // counters for DMA replies


extern DeviceArgs args_dev;

/// \brief Initialize and hold some constants on device.
void initializeConstantsOnDevice(DeviceConstants *constants_dev) {
    int tid = athread_get_id(-1);

    rank = constants_dev->rank;
    size_t dev_n = constants_dev->dev_n;
    block_size = dev_n / 64;

    if (tid == 0) {
        printf("[Rank %d][tid %d] block size = %lu, range = [%lu, %lu]\n",
            rank, tid, block_size, tid * block_size, (tid + 1) * block_size - 1);
        printf("[Rank %d][tid %d] iterations = %lu, res = %lu\n",
            rank, tid, block_size / THREAD_BATCH_SIZE + 1, block_size % THREAD_BATCH_SIZE);
    }
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
    size_t num_iterations = block_size / THREAD_BATCH_SIZE;
    size_t res_points = block_size % THREAD_BATCH_SIZE;

    // Allocation
    size_t batch_size = THREAD_BATCH_SIZE;
    size_t batch_size_next;
    size_t istart = tid * block_size;
    int buf_id;
    int buf_id_next;

    for (buf_id = 0; buf_id < 2; buf_id++) {
        fx[buf_id] = (double *)ldm_malloc((THREAD_BATCH_SIZE+3)*sizeof(double));
        x[buf_id] = (double *)ldm_malloc((THREAD_BATCH_SIZE+3)*sizeof(double));
        alpha[buf_id] = (double *)ldm_malloc((THREAD_BATCH_SIZE+3)*sizeof(double));
        beta[buf_id] = (double *)ldm_malloc((THREAD_BATCH_SIZE+3)*sizeof(double));
    }

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

    // Fetch data from the host for the first buffer
    batch_size_next = THREAD_BATCH_SIZE;

    if (num_iterations == 0) {
        batch_size_next = res_points;
    }

    num_gets[0] = 0;
    athread_get(PE_MODE, args_dev.x + istart, &x[0][0], (batch_size_next+2)*8, (void *)&num_gets[0], 0, 0, 0);
    athread_get(PE_MODE, args_dev.alpha + istart, &alpha[0][0], (batch_size_next+1)*8, (void *)&num_gets[0], 0, 0, 0);
    athread_get(PE_MODE, args_dev.beta + istart, &beta[0][0], (batch_size_next+1)*8, (void *)&num_gets[0], 0, 0, 0);

    // Skip the first send
    num_puts[1] = 1;

    size_t j;
    for (j = 0; j <= num_iterations; j++) {
        istart = tid * block_size + j * THREAD_BATCH_SIZE;  // the start index of the block of this task
        buf_id = j % 2;
        buf_id_next = (j+1) % 2;

        if (j == num_iterations) {
            batch_size = res_points;
            batch_size_next = 0;
        }
        else { // need to fetch the next batch of data
            if (j == num_iterations - 1) {
                batch_size_next = res_points;
            }

            num_gets[buf_id_next] = 0;
            athread_get(PE_MODE, args_dev.x + istart, &x[buf_id_next][0], (batch_size_next+2)*8, (void *)&num_gets[buf_id_next], 0, 0, 0);
            athread_get(PE_MODE, args_dev.alpha + istart, &alpha[buf_id_next][0], (batch_size_next+1)*8, (void *)&num_gets[buf_id_next], 0, 0, 0);
            athread_get(PE_MODE, args_dev.beta + istart, &beta[buf_id_next][0], (batch_size_next+1)*8, (void *)&num_gets[buf_id_next], 0, 0, 0);
        }

        #ifdef DEBUG
        if (tid == DEBUGTID) {
            printf("[Rank %d][tid %d] iteration %lu, istart = %lu\n", rank, tid, j, istart);
            printf("[Rank %d][tid %d] buffer id = %d (current), %d (next)\n", rank, tid, buf_id, buf_id_next);
            printf("[Rank %d][tid %d] batch size = %lu (current), %lu (next)\n", rank, tid, batch_size, batch_size_next);
        }
        #endif

        #ifdef DEBUG
        if (tid == DEBUGTID) {
            printf("[Rank %d][tid %d] waiting for buffer %d (current)\n", rank, tid, buf_id);
        }
        #endif
        // Wait for the data fetched during the previous loop (stored in the current buffer)
        while (num_gets[buf_id] != 3);

        size_t i;
        #ifdef USE_DEVICE_SIMD
        // Computing (simd)
        for (i = 0; i < (batch_size/4)*4; i+=4) {
            simd_load(vx, &(x[buf_id][i]));
            simd_load(vbeta, &(beta[buf_id][i]));
            simd_loadu(vx1, &(x[buf_id][i+1]));
            simd_loadu(vbeta1, &(beta[buf_id][i+1]));
            simd_loadu(vx2, &(x[buf_id][i+2]));
            simd_loadu(valpha1, &(alpha[buf_id][i+1]));

            vfx = simd_vmad(valpha1, vx2, simd_vmss(vx0, vbeta * vx, vx1 * simd_vmss(vx0, vbeta1, vk)));

            simd_store(vfx, &(fx[buf_id][i]));
        }

        for (i = (batch_size/4)*4; i < batch_size; i++) {
            fx[buf_id][i] = x0 * (beta[buf_id][i] * x[buf_id][i])
                            - x[buf_id][i+1] * (x0 * beta[buf_id][i+1] - k)
                            + alpha[buf_id][i+1] * x[buf_id][i+2];
        }
        #else
        // Computing
        for (i = 0; i < batch_size; i++) {
            fx[buf_id][i] = x0 * (beta[buf_id][i] * x[buf_id][i])
                            - x[buf_id][i+1] * (x0 * beta[buf_id][i+1] - k)
                            + alpha[buf_id][i+1] * x[buf_id][i+2];
        }
        #endif

        #ifdef DEBUG
        if (tid == DEBUGTID) {
            printf("[Rank %d][tid %d] sending back buffer %d (current)\n", rank, tid, buf_id);
        }
        #endif
        // Send data back to the host
        num_puts[buf_id] = 0;
        athread_put(PE_MODE, &fx[buf_id][0], args_dev.fx + istart, batch_size*8, (void *)&num_puts[buf_id], 0, 0);

        #ifdef DEBUG
        if (tid == DEBUGTID) {
            printf("[Rank %d][tid %d] waiting for sending back buffer %d (previous)\n", rank, tid, buf_id_next);
        }
        #endif
        // Wait for the data sent during the previous loop (stored in another buffer)
        while (num_puts[buf_id_next] != 1);
    }

    // Wait for the last batch of data to be sent
    while (num_puts[buf_id] != 1);

    // Deallocation
    for (buf_id = 0; buf_id < 2; buf_id++) {
        ldm_free(fx[buf_id], (THREAD_BATCH_SIZE+3)*sizeof(double));
        ldm_free(x[buf_id], (THREAD_BATCH_SIZE+3)*sizeof(double));
        ldm_free(alpha[buf_id], (THREAD_BATCH_SIZE+3)*sizeof(double));
        ldm_free(beta[buf_id], (THREAD_BATCH_SIZE+3)*sizeof(double));
    }
}
