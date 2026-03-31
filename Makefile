

GGG = g++ -O3 -Wall -march=native -fomit-frame-pointer -fexpensive-optimizations
# GGG = clang++ -O3 -Wall -march=native -fomit-frame-pointer

OBJ = cubic_primality_main.o \
      cubic_primality.o \
      cubic_primality_alloc.o \
      cubic_primality_precompute.o \
      expression_parser.a

cubic: $(OBJ)
	$(GGG) -static -o cubic $(OBJ) -lgmp -lpthread -lm

cubic_primality_main.o: cubic_primality_main.cpp cubic_primality.h cubic_primality_alloc.h bison.gmp_expr.tab.h
	$(GGG) -c -o cubic_primality_main.o cubic_primality_main.cpp

cubic_primality_alloc.o: cubic_primality_alloc.cpp cubic_primality_alloc.h
	$(GGG) -c -o cubic_primality_alloc.o cubic_primality_alloc.cpp

cubic_primality.o: cubic_primality.cpp cubic_primality.h
	$(GGG) -c -o cubic_primality.o cubic_primality.cpp

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

check: cubic
	./cubic -st

clean:
	rm -f ./cubic $(OBJ) bison.gmp_expr.o bison.gmp_expr.tab.c bison.gmp_expr.tab.h lex.gmp_expr.o lex.gmp_expr.c


