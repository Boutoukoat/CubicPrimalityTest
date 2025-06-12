#pragma once

// -----------------------------------------------------------------------
// Cubic primality test
//
// interface to a multithreaded version for heavy computations
// -----------------------------------------------------------------------

#include "gmp.h"
#include <stdbool.h>
#include <stdint.h>

struct mod_precompute_t
{
    mpz_t a;
    mpz_t b;
    mpz_t m;
    uint64_t n;
    uint64_t n32;
    uint64_t n2;
    bool special_case;
    bool montg;
    bool proth;
    bool power2me;
    bool power2pe;
    uint64_t e;
};

void mpz_inner_multithread_exponentiate(mpz_t s, mpz_t t, mpz_t u, mpz_t e, uint64_t a, mod_precompute_t *p);
void mpz_mod_fast_reduce(mpz_t r, struct mod_precompute_t *p);

