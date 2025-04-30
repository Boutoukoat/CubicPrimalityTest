// -----------------------------------------------------------------------
// Cubic primality test
// -----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <string.h>

#include "cubic_primality.h"
#include "bison.gmp_expr.h"

static void cubic_primality_file(char *name)
{
	long prime_count = 0;
	long composite_count = 0;
	FILE *f = fopen(name, "rt");
	if (f) {
		const int buff_len = 1000000;   // max line length
		char buff[buff_len];
		mpz_t v;
		mpz_init(v);

		// get a line
		while (fgets(buff, buff_len, f)) {
			// trim trailing spaces
			int len = strlen(buff);
			while (len > 0 && isspace(buff[len - 1])) {
				len--;
			}
			buff[len] = 0;
			// trim leading spaces
			char *pt = buff;
			while (isspace(*pt)) {
				pt++;
			}
			if (*pt && *pt != '#')	// discard empty lines or comments
			{
				mpz_expression_parse(v, pt);
				bool is_prime = mpz_cubic_primality(v);
				prime_count += (is_prime == true);
				composite_count += (is_prime == false);
			}
		}
		mpz_clear(v);
	}
	printf("File %s done, %ld primes, %ld composites\n", name, prime_count, composite_count);
}

int main(int argc, char **argv)
{
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-st")) {
			// internal sanity self tests
			cubic_primality_self_test();
			printf("Self tests completed\n");
			exit(0);
		} else if (!strcmp(argv[i], "-f")) {
			cubic_primality_file(argv[++i]);
			exit(0);
		} else {
			// command line argument must be a number, or an expression
			struct timespec ts1, ts2;
			mpz_t n;
			mpz_init(n);

			// Read an expression from the command line
			// Supported operators are +/-*^() with usual precedence.
			mpz_expression_parse(n, argv[i]);
			// gmp_printf("Test n=%Zd.\n",n);

			clock_gettime(CLOCK_REALTIME, &ts1);

			// Get the number of steps. Assuming the Collatz conjecture, it is not an infinite task. 
			bool is_prime = mpz_cubic_primality(n);
			clock_gettime(CLOCK_REALTIME, &ts2);

			// Display the input number, the number of steps, and the time it took to run it, in milliseconds.
			double diff = ts2.tv_sec - ts1.tv_sec;
			diff *= 1e9;
			diff += ts2.tv_nsec;
			diff -= ts1.tv_nsec;
			diff /= 1e6;
			printf("%s %s, time=%12.3f msecs.\n", argv[i], (is_prime ? "might be prime" : "is composite for sure"),
			       diff);
			fflush(stdout);
			mpz_clear(n);
		}
	}

	return 0;
}
