#include "heat.h"

#define NB 100

#define min(a,b) ( ((a) < (b)) ? (a) : (b) )

double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
	double diff, sum=0.0;
	int nbx, bx, nby, by;

	nbx = omp_get_max_threads();//NB;
	bx = sizex/nbx + ((sizex%nbx) ? 1 : 0);//sizex/nbx;
	nby = 1;//NB;
	by = sizey/nby;

	#pragma omp parallel for reduction(+:sum) private(diff)
	for (int ii=0; ii<nbx; ii++) {
		for (int jj=0; jj<nby; jj++)  {
			for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) {
				for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
					utmp[i*sizey+j]= 0.25 * (u[ i*sizey	    + (j-1) ]+  // left
					                         u[ i*sizey     + (j+1) ]+  // right
					                         u[ (i-1)*sizey + j     ]+  // top
					                         u[ (i+1)*sizey + j     ]); // bottom
					diff = utmp[i*sizey+j] - u[i*sizey + j];
					sum += diff * diff; 
				}
			}
		}
	}

	return sum;
}

