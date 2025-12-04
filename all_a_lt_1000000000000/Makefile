
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

# -----------------------------------------
#  parallel version
# -----------------------------------------
test: all
	./all_a -max 200000

unit_test: math_tests.cpp math_utils.cpp
	g++ -O3 -o unit_tests math_tests.cpp
	./unit_tests



