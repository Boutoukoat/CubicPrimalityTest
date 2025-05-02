// Requires file cubic-sq-disc-B.dat containing start and end like:
// 11 100000000
//
// Compile with: g++ -03 -o cubic-sq-disc-B cubic-sq-disc-B.cpp -lprimesieve
//
// Author : Paul Underwood http://www.worldofprimes.co.uk/corner
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "primesieve/soe/PrimeSieve.h" // ver 4.4

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t n) {
	uint64_t res = a;
	__asm__("mulq %1        \n\t"
	      "divq %2        \n\t"
	      "movq %%rdx, %0 \n\t"
		: "+a"(res)
		: "g"(b), "g"(n)
		: "%rdx", "cc");
	return res;
}

uint64_t gcd_test( uint64_t P, uint64_t n ) { 
  uint64_t r;
  if ( P == 0 ) return ( 0 );
  r = n % P;
  if ( r == 0 ) return ( P );
  return ( gcd_test( r, P ) );
}

int fermat_with_mask ( uint64_t a, uint64_t n, uint64_t bit_start ) {
  uint64_t result = a;
  for( uint64_t bit_mask = bit_start; bit_mask; bit_mask = bit_mask >> 1 ) {
    result *= result;
    result %= n;
    if( bit_mask & n ) {
      result *= a ;
      result %= n;
    }
  }
  return ( result );
}

uint64_t n, n_end, base_bit_mask;

void callback ( unsigned long prime ) {
  FILE *log_file;
  uint64_t k, a, k_limit, diff;
  uint64_t s, t, u;
  uint64_t S, T, U;
  uint64_t s2, t2, u2;
  uint64_t st, tu, us;
  uint64_t bit_start, bit_mask, nm1;
  for ( ; n < prime; n += 2ULL ) {
    bit_start = base_bit_mask;
    while ( bit_start <= n ) bit_start = bit_start << 1;
    base_bit_mask = bit_start;
    bit_start = bit_start >> 2;
    k_limit = ( n + 1ULL ) >> 1;
    nm1 = n - 1ULL;
    diff = 0ULL;
    a = 7ULL;
    for ( k = 1ULL; k <= k_limit; k++ ) {
      a += diff;
      diff += 2ULL;
      if ( a >= n ) a -= n;
      if ( gcd_test ( a * ( 2ULL * k - 1ULL ) % n, n ) != 1ULL ) continue;
      if ( fermat_with_mask ( a, n, bit_start ) != a ) continue;
      s = 0ULL;
      t = 1ULL;
      u = 0ULL;
      for ( bit_mask = bit_start; bit_mask; bit_mask = bit_mask >> 1 ) {
        s2 = mulmod ( s, s, n );
	s2 = a * s2;
	t2 = t * t;
	u2 = u * u;
        st = mulmod ( s, t, n );
	st = a * st;
	tu = t * u;
	us = u * s;
	st = st << 1;
	tu = tu << 1;
	us = us << 1;
	if ( bit_mask & nm1 ) {
	  u = a * ( ( s2 + us + t2 ) %n );
	  s = s2 + st + tu;
	  t = u + u2 + st;
	} else {
	  s = s2 + us + t2;
	  t = s2 + st + tu;
	  u = st + u2;
	}
	s %= n;
	t %= n;
	u %= n;
      }
      if ( t == 0ULL && s == 0ULL && u == 1ULL ) continue;
      S = s;
      T = t;
      U = u;
      s2 = mulmod ( s, s, n );
      s2 = a * s2;
      t2 = t * t;
      u2 = u * u;
      st = mulmod ( s, t, n );
      st = a * st;
      tu = t * u;
      us = u * s;
      st = st << 1;
      tu = tu << 1;
      us = us << 1;
      s = s2 + us + t2 + S;
      t = s2 + st + tu + T;
      u = st + u2 + U + 1;
      s %= n;
      t %= n;
      u %= n;
      if ( s == nm1 && t == 1ULL && u == a ) {
	printf ( "[%llu, %llu]\n", n, k );
	log_file = fopen ( "cubic-sq-disc-B.log", "a" );
	fprintf ( log_file, "[%llu, %llu]\n", n, k );
	fclose ( log_file );
      }
    }
  }
  n += 2ULL;
}

int main ( int argc, char *argv[] ) {
  FILE *dat_file;
  PrimeSieve ps;
  dat_file = fopen ( "cubic-sq-disc-B.dat", "r" );
  fscanf ( dat_file, "%llu %llu", &n, &n_end );
  fclose( dat_file );
  if ( n % 2ULL == 0 ) n++;
  base_bit_mask = 1ULL;
  while ( base_bit_mask <= n ) base_bit_mask = base_bit_mask << 1;
  while ( n < n_end ) {
    ps.generatePrimes ( n, n + 1000, callback );
    dat_file = fopen ( "cubic-sq-disc-B.dat", "w" );
    fprintf ( dat_file, "%llu %llu\n", n, n_end );
    fclose ( dat_file );
  }
}

