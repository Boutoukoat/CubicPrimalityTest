

GGG = g++ -O3 -march=native -fomit-frame-pointer -fexpensive-optimizations

OBJ = cubic_primality_main.o cubic_primality.o bison.gmp_expr.o lex.gmp_expr.o


cubic: $(OBJ)
	$(GGG) -static -o cubic $(OBJ) -lgmp

cubic_primality_main.o: cubic_primality_main.cpp cubic_primality.h
	$(GGG) -c -o cubic_primality_main.o cubic_primality_main.cpp

cubic_primality.o: cubic_primality.cpp cubic_primality.h
	$(GGG) -c -o cubic_primality.o cubic_primality.cpp

bison.gmp_expr.o : bison.gmp_expr.tab.c bison.gmp_expr.h
	$(GGG) -c -o bison.gmp_expr.o bison.gmp_expr.tab.c

bison.gmp_expr.tab.c bison.gmp_expr.tab.h : parser.y
	bison -d parser.y

lex.gmp_expr.o : lex.gmp_expr.c
	$(GGG) -c -o lex.gmp_expr.o lex.gmp_expr.c

lex.gmp_expr.c : parser.l bison.gmp_expr.tab.h
	flex parser.l

clean:
	rm -f ./cubic $(OBJ)


