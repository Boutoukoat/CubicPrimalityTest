#pragma once

// -----------------------------------------------------------------------
// Cubic primality test
//
// mpz_cubic_primality():
//    true: composite for sure
//    false: might be prime
//
// cubic_primality_self_test()
//    simplified unit tests to detect a possible compiler/platform issue.
//    assert when fail (this should not happen).
// -----------------------------------------------------------------------

#include "gmp.h"
#include <stdbool.h>

bool mpz_cubic_primality(mpz_t v, bool verbose = false, bool multithread = false);
void cubic_primality_self_test(void);
