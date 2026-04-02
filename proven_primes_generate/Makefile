
# GGG=clang -O3
GGG=g++ -O3

all: ./generate

./generate: generate.cpp bison.gmp_expr.h expression_parser.a
	$(GGG) -o ./generate generate.cpp expression_parser.a -lgmp

example: ./generate
	# 3 threads, 1000 primes of size 256 bits, starting from well known prime 103
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

gigantic: ./generate
	# 32 threads, compute 80 primes of size up to 33500 bits (approx 10100 digits), 
	# starting from a 10000 digit proven prime
	# (This takes hours)
	./generate -t 32 -c 80 -verbose -b 33500 -s " (10^10000 - 1)/3 - 10^6333" 

expression_parser.a : bison.gmp_expr.o lex.gmp_expr.o bison.gmp_expr.tab.h
	ar vr expression_parser.a bison.gmp_expr.o lex.gmp_expr.o

bison.gmp_expr.o : bison.gmp_expr.tab.c bison.gmp_expr.h
	$(GGG) -c -o bison.gmp_expr.o bison.gmp_expr.tab.c

bison.gmp_expr.tab.c bison.gmp_expr.tab.h : parser.y
	bison -d parser.y

lex.gmp_expr.o : lex.gmp_expr.c
	$(GGG) -Wno-unused-function -c -o lex.gmp_expr.o lex.gmp_expr.c

lex.gmp_expr.c : parser.l bison.gmp_expr.tab.h
	flex parser.l


