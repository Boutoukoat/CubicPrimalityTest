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

bool cubic_primality(mpz_t v);
