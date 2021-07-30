#ifndef UTILS_H_
#define UTILS_H_

#include <stdlib.h>
#include <string.h>


static inline unsigned long rpcc() {
    unsigned long time;
    asm("rtc %0": "=r" (time) : );
    return time;
}


/// \brief Allocates aligned memory.
void *alloc_aligned(size_t alignment, size_t size) {
    void *ptr;
    posix_memalign(&ptr, alignment, size);
    memset(ptr, 0, size);
    return ptr;
}


#endif  // UTILS_H_
