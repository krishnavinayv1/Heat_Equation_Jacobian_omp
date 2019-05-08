
#ifndef _EXTRAE_
   #define _EXTRAE_ 0
#else 
   #define _EXTRAE_ 1
#endif

#include <stdio.h>
#include <stdlib.h>
#include "heat.h"
#if _EXTRAE_
#include "extrae_user_events.h"
#endif

void usage( char *s )
{
	fprintf(stderr, 
		"Usage: %s <input file> [result file]\n\n", s);
}

inline void swap(double** a, double** b) {
	double *tmp = *a;
	*a = *b;
	*b = tmp;
}

int main( int argc, char *argv[] )
{
	unsigned iter;
	FILE *infile, *resfile;
	char *resfilename;

    //omp_set_max_threads(8);


	int max_threads = omp_get_max_threads();


	algoparam_t param;
	int np;

	double runtime, flop;
	double residual=0.0;


	if( argc < 2 ) {
		usage( argv[0] );
		return 1;
	}

	// check input file
	if(!(infile=fopen(argv[1], "r"))) {
		fprintf(stderr, "\nError: Cannot open \"%s\" for reading.\n\n", argv[1]);
		usage(argv[0]);
		return 1;
	}

	// check result file
	resfilename = (argc>=3) ? argv[2]:"heat.ppm";

	if(!(resfile=fopen(resfilename, "w"))) {
		fprintf(stderr, "\nError: Cannot open \"%s\" for writing.\n\n", resfilename);
		usage(argv[0]);
		return 1;
	}

	// check input
	if(!read_input(infile, &param)) {
		fprintf(stderr, "\nError: Error parsing input file.\n\n");
		usage(argv[0]);
		return 1;
	}

	// print parameters
	printf("Max Threads       : %d\n", max_threads);
	print_params(&param);

	// initialize solver
	if(!initialize(&param)) {
		fprintf(stderr, "Error in Solver initialization.\n\n");
		usage(argv[0]);
		return 1;
	}


	np = param.resolution + 2;
	
#if _EXTRAE_
	Extrae_init();
#endif

	// starting time
	runtime = wtime();

	iter = 0;
	while(1) {
	switch( param.algorithm ) {
		case 0: 
			residual = relax_jacobi(param.u, param.uhelp, np, np);
			swap(&param.u, &param.uhelp);
			break;
		}
                if(iter<2)
                {
                    for(int i=0;i<256;i++)
                    {
                        for(int j=0;j<256;j++)
                            printf("%.3f ",param.u[i*10+j]);
                        printf("\n");
                    }
                printf("%.5f\n",residual);
                }
		iter++;
		// solution good enough ?
		if (residual < 0.00005) 
                {
                    break;
                }

		// max. iteration reached ?
		if (param.maxiter>0 && iter>=param.maxiter) break;
	}

	// Flop count after iter iterations
	flop = iter * 11.0 * param.resolution * param.resolution;
	// stopping time
	runtime = wtime() - runtime;
#if _EXTRAE_
	Extrae_fini();
#endif


	fprintf(stdout, "Time: %04.3f ", runtime);
	fprintf(stdout, "Convergence to residual=%f: %d iterations\n", residual, iter);

	// for plot...
	coarsen(param.u, np, np, param.uvis, param.visres+2, param.visres+2);
  
	write_image(resfile, param.uvis, param.visres+2, param.visres+2 );

	finalize(&param);

	return 0;
}
