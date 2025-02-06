TRACE=NDEBUG
#SCIPDIR=/usr
#LDFLAGS=-L $(SCIPDIR)
LDFLAGS=
CFLAGS=-g -std=c11 -Wall -D$(TRACE) -D SCIP_VERSION_MAJOR

bin/mochila: bin/cmain.o bin/probdata_mochila.o bin/problem.o bin/heur_problem.o bin/heur_myrounding.o  bin/heur_aleatoria.o bin/heur_grasp.o
	gcc -o bin/mochila-$(TRACE) bin/cmain.o bin/probdata_mochila.o bin/problem.o bin/heur_problem.o bin/heur_myrounding.o bin/heur_aleatoria.o bin/heur_grasp.o -lscip -lm

bin/cmain.o: src/cmain.c
	gcc $(LDFLAGS) $(CFLAGS) -c -o bin/cmain.o src/cmain.c

bin/problem.o: src/problem.c
	gcc $(LDFLAGS) $(CFLAGS) -c -o bin/problem.o src/problem.c

bin/heur_problem.o: src/heur_problem.c src/heur_problem.h
	gcc $(LDFLAGS) $(CFLAGS) -c -o bin/heur_problem.o src/heur_problem.c

bin/probdata_mochila.o: src/probdata_mochila.c src/probdata_mochila.h
	gcc $(LDFLAGS) $(CFLAGS) -c -o bin/probdata_mochila.o src/probdata_mochila.c

bin/heur_myrounding.o: src/heur_myrounding.c src/heur_myrounding.h
	gcc $(LDFLAGS) $(CFLAGS) -c -o bin/heur_myrounding.o src/heur_myrounding.c

bin/heur_aleatoria.o: src/heur_aleatoria.c src/heur_aleatoria.h
	gcc $(LDFLAGS) $(CFLAGS) -c -o bin/heur_aleatoria.o src/heur_aleatoria.c

bin/heur_grasp.o: src/heur_grasp.c src/heur_grasp.h
	gcc $(LDFLAGS) $(CFLAGS) -c -o bin/heur_grasp.o src/heur_grasp.c

.PHONY: clean

clean:
	rm -f bin/*.o bin/mochila

