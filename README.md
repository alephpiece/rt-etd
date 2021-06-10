RT ETD
======

# Serial versions

```bash
serial:origin  # the origin without MPI
serial:v1      # from origin, remove unnecessary branches and variables, rearrange formula evaluation
serial:v1a     # from v1, modify the computation of F(0)
serial:v2      # from v1, pad vector F with 1 leading space
serial:v2a     # from v2, pad vector F with 2 leading space
serial:v3      # from v2, support simd reduction
serial:v4      # from v1, combine adjacent for loops together
```

# MPI versions

```bash
mpi:v1      # from serial:v3
```
