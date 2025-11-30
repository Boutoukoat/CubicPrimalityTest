
# -----------------------------------------
# single thread version
# -----------------------------------------
all:
	g++ -O3 -o all_a all_a_check.cpp

# -----------------------------------------
#  parallel version
# -----------------------------------------
omp:
	g++ -O3 -fopenmp -o all_a all_a_check.cpp -l gomp

test: all
	./all_a -p 10000 -Q 2000

