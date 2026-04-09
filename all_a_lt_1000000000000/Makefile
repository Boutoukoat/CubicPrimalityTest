
GGG=clang++ -O3 -std=c++2b -Wno-deprecated-volatile

# -----------------------------------------
# single thread version
# -----------------------------------------
all: ./all_a ./all_ab

./all_a: ./all_a_check.cpp
	$(GGG) -o ./all_a ./all_a_check.cpp

./all_ab: ./all_ab_check.cpp
	$(GGG) -o ./all_ab ./all_ab_check.cpp -lm

# -----------------------------------------
#  parallel version
# -----------------------------------------
omp:
	g++ -O3 -fopenmp -o ./all_a ./all_a_check.cpp -l gomp

# -----------------------------------------
#  parallel version
# -----------------------------------------
test: ./all_a
	./all_a -max 200000

unit_test: math_tests.cpp math_utils.cpp
	g++ -O3 -o unit_tests math_tests.cpp
	./unit_tests



