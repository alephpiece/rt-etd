#ifndef DEVICE_STRUCTS_H_
#define DEVICE_STRUCTS_H_

#include <stddef.h>


/// \struct DeviceArgs
/// \brief A struct of time-dependent values.
#pragma pack (4)
typedef struct {

    double *fx;     // F()
    double *x;      // C, A, B, or D
    double *alpha;  // \alpha
    double *beta;   // \beta
    double x0;      // C0, A0, B0, or D0
    double k;   // k

} DeviceArgs;


/// \struct DeviceArgs
/// \brief A struct of constant values.
#pragma pack (4)
typedef struct {

    int rank;    // MPI rank
    size_t dev_n;   // Total number of data points for device

} DeviceConstants;


#endif  // DEVICE_STRUCTS_H_
