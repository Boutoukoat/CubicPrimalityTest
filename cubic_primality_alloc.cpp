// -----------------------------------------------------------------------
// Cubic primality test
//
// allocation functions with alignment and false-sharing care
// -----------------------------------------------------------------------

#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "cubic_primality_alloc.h"

void *cubic_allocate_function(size_t alloc_size)
{
    void *ptr = aligned_alloc(64, alloc_size + 64);
    if (!ptr)
    {
        // catastrophic failure
        perror("aligned_alloc");
        abort();
    }
    return ptr;
}

void *cubic_reallocate_function(void *ptr, size_t old_size, size_t new_size)
{
    ptr = realloc(ptr, new_size + 64);
    if (!((uintptr_t)ptr & 63))
    {
        return ptr;
    }
    void *newptr = aligned_alloc(64, new_size + 64);
    memcpy(newptr, ptr, new_size);
    return newptr;
}

void cubic_free_function(void *ptr, size_t size)
{
    free(ptr);
}
