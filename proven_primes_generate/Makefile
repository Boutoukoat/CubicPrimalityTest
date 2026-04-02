


all: ./generate example

./generate: generate.cpp
	g++ -O3 -o ./generate generate.cpp -lgmp

example: ./generate
	# 3 threads, 128 primes of size 256 bits, starting from well known prime 103
	@echo
	@echo "generate 1000 primes"
	@echo
	./generate -t 3 -c 1000 -b 256 -s 103 | tee example_256.txt
	@echo
	@echo "Check primality"
	@echo
	../cubic -f example_256.txt
	@echo
	@echo "done"
	@echo

