
# -----------------------------------------
# single thread version
# -----------------------------------------
all:
	g++ -O3 -o gcd_ext gcd_ext_check.cpp

# -----------------------------------------
#  parallel version
# -----------------------------------------
omp:
	g++ -O3 -fopenmp -o gcd_ext gcd_ext_check.cpp -l gomp

test: all
	./gcd_ext -p 1000000000 -l 10000

