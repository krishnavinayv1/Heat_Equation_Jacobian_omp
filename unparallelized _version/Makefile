CC     =  gcc 
TAREADOR_CFLAGS = -Wall -std=c99 -O0 -g

OPT3 	= -O3
OPTG0 	= -g -O0
CFLAGS  = -Wall -std=c99
LFLAGS  = -lm 

ALL	= heat heat-tareador 
all: $(ALL)

misc.o: misc.c
	$(CC) -c  $(OPT3) $(CFLAGS) $< -o $@

heat: heat.c solver.c misc.o
	$(CC) $(OPT3) $(CFLAGS) $+ $(LFLAGS) -fPIC -o $@


clean:
	rm -fr $(ALL) .tareador.* *.o *.xdot *.ppm *.out *.err *.prv *.pcf *.row 

