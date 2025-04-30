#pragma once

// -----------------------------------------------------------------------
// Cubic primality test
//
// true: composite for sure
// false: might be prime
//
// -----------------------------------------------------------------------

#include <stdbool.h>
#include "gmp.h"

bool mpz_cubic_primality(mpz_t v);
void cubic_primality_self_test(void);
