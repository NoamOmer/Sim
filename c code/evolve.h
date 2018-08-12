/*
	Header for magnetization evolution C functions.
	Contains infrastructure API.

	Version 1.0 16-Jul-2007 Noam Ben-Eliezer
	Version 2.0 22-Oct-2014 Noam Ben-Eliezer
*/
# define MAC

#include "mex.h"
#include "matrix.h"

#ifndef MAC
	#include "malloc.h"
#endif

#include <math.h>
#include <time.h>

typedef double ** Matrix;

double *unit_matrix(double *R_mat)
{
	R_mat[0] = 1;
	R_mat[1] = 0;
	R_mat[2] = 0;
	R_mat[3] = 0;
	R_mat[4] = 1;
	R_mat[5] = 0;
	R_mat[6] = 0;
	R_mat[7] = 0;
	R_mat[8] = 1;

	return R_mat;
}

/* ================================================
                  Debugging
   ================================================ */
void print_2D_M_mat(Matrix Mx, Matrix My, Matrix Mz, int Nx, int Ny, int n_jumps, char *context)
{
	mxArray *input_mxArr[7];
	int cols_idx, rows_idx;
	double *M_outx, *M_outy, *M_outz;

	input_mxArr[0] = mxCreateDoubleMatrix(1,Nx*Ny, mxREAL);
	input_mxArr[1] = mxCreateDoubleMatrix(1,Nx*Ny, mxREAL);
	input_mxArr[2] = mxCreateDoubleMatrix(1,Nx*Ny, mxREAL);
	input_mxArr[3] = mxCreateDoubleScalar(Nx);
	input_mxArr[4] = mxCreateDoubleScalar(Ny);
	input_mxArr[5] = mxCreateDoubleScalar(n_jumps);
	input_mxArr[6] = mxCreateString(context);

	M_outx = mxGetPr(input_mxArr[0]);
	M_outy = mxGetPr(input_mxArr[1]);
	M_outz = mxGetPr(input_mxArr[2]);

	// Reshape the output magnetization 2D matrices into 1D vectors
	for (cols_idx = 0; cols_idx < Ny; cols_idx++) {
		for (rows_idx = 0; rows_idx < Nx; rows_idx++) {
			M_outx[cols_idx*Nx + rows_idx]  = Mx[rows_idx][cols_idx];
			M_outy[cols_idx*Nx + rows_idx]  = My[rows_idx][cols_idx];
			M_outz[cols_idx*Nx + rows_idx]  = Mz[rows_idx][cols_idx];
		}
	}

	mexCallMATLAB(0, NULL, 7, input_mxArr, "plot_2D_M_mat_from_C_ifs");

	mxDestroyArray(input_mxArr[0]);
	mxDestroyArray(input_mxArr[1]);
	mxDestroyArray(input_mxArr[2]);
	mxDestroyArray(input_mxArr[3]);
	mxDestroyArray(input_mxArr[4]);
	
	return;
}

/* ================================================
                  Rotation matrices
   ================================================ */
// This function should have had an integer return value to indicate status instead of (double *)
double *rot_RHR(double x, double y, double z, double phi, double *R_mat)
{
	double t,s,c;
//	double tx,ty;
	/* Skip check
	if ( (abs(sqrt(sum(unit_vec.^2)) - 1) > 1E-13) || isnan(x) || isnan(y) || isnan(z) )
		error('rot: x,y,z do not represent a unit vector (x=%d, y=%d, z=%d). Exiting',x,y,z);
	end;
	*/

	t = 1 - cos(phi);
	s = sin(phi);
	c = cos(phi);

	R_mat[0] = t*x*x + c;
	R_mat[1] = t*x*y - s*z;
	R_mat[2] = t*x*z + s*y;
	R_mat[3] = t*x*y + s*z;
	R_mat[4] = t*y*y + c;
	R_mat[5] = t*y*z - s*x;
	R_mat[6] = t*x*z - s*y;
	R_mat[7] = t*y*z + s*x;
	R_mat[8] = t*z*z + c;

/*
	tx = t*x,
	ty = t*y;
	R_mat[0] = tx*x  + c;
	R_mat[1] = tx*y  - s*z;
	R_mat[2] = tx*z  + s*y;
	R_mat[3] = tx*y  + s*z;
	R_mat[4] = ty*y  + c;
	R_mat[5] = t*y*z - s*x;
	R_mat[6] = tx*z  - s*y;
	R_mat[7] = ty*z  + s*x;
	R_mat[8] = t*z*z + c;
*/
	return R_mat;
}

double *rot_LHR(double x, double y, double z, double phi, double *R_mat)
{
	double t,s,c;

	phi = -phi;
	t = 1 - cos(phi);
	s = sin(phi);
	c = cos(phi);

	R_mat[0] = t*x*x + c;
	R_mat[1] = t*x*y - s*z;
	R_mat[2] = t*x*z + s*y;
	R_mat[3] = t*x*y + s*z;
	R_mat[4] = t*y*y + c;
	R_mat[5] = t*y*z - s*x;
	R_mat[6] = t*x*z - s*y;
	R_mat[7] = t*y*z + s*x;
	R_mat[8] = t*z*z + c;

	return R_mat;
}

void print_R_mat(double *R)
{
	mexPrintf("R = %6.4f  %6.4f  %6.4f\n    %6.4f  %6.4f  %6.4f\n    %6.4f  %6.4f  %6.4f\n\n",
		       R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7],R[8]);
	return;
}

# ifndef MAC
void wait(int seconds)
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
}
#endif

/* ================================================
                    Relaxation API
   ================================================ */
void T1_relaxation(double *Mx, double *My, double *Mz, double M0z, double t, double T1)
{
	double newMz = *Mz;

	newMz = newMz*exp(-t/T1) + M0z*(1 - exp(-t/T1));
//mexPrintf("*Mz = %3.3f, M0z = %3.3f, newMz = %3.3f, t=%3.3f, T1=%3.3f\n",*Mz,M0z,newMz,t,T1);
	*Mz = newMz;

	return;
}

void T2_relaxation(double *Mx, double *My, double *Mz, double t, double T2)
{
	double T2_relax_factor = exp(-t/T2),
	       newMx = *Mx,
		   newMy = *My;

	newMx = newMx * T2_relax_factor;
	newMy = newMy * T2_relax_factor;
	*Mx = newMx;
	*My = newMy;

	return;
}

/* ================================================
           Matrix  manipulations  Interface        
   ================================================ */

Matrix NewMatrix(int rows, int cols)
{
    int i;
    Matrix newM;
	
	newM = (double **) mxCalloc(rows, sizeof(double *));
	if (!newM)
		mexErrMsgTxt("Error allocating new matrix");

	for(i = 0; i < rows; i++) {
		newM[i] = (double *) mxCalloc(cols, sizeof(double));
		if (!newM[i])
			mexErrMsgTxt("Error allocating new matrix row");
	}

return newM;
}

void FreeMatrix(Matrix mat, int rows)
{
    int i;
    for(i = 0; i < rows; i++)
		mxFree(mat[i]);
    mxFree(mat);
}

void TransposeMatrix(Matrix inM, Matrix outM, int cols, int rows)
{
    int tempI, tempJ;
    for(tempI=0; tempI < rows; tempI++)
	for(tempJ=0; tempJ < cols; tempJ++)
	    outM[tempI][tempJ] = inM[tempJ][tempI];    
}

void MultMatrix(Matrix firstM, Matrix secondM, Matrix outM, int firstrows, int cols, int secondcols)
{
    int i,j,k;
    double sum;

    for(i=0; i < secondcols; i++)
	for(j=0; j < firstrows; j++)
	{
	    sum = 0.0;
	    for(k=0; k < cols; k++)
		sum += firstM[j][k] * secondM[k][i];
	    outM[j][i] = sum;
	}
}

double InvertMatrix(Matrix mat, int actual_size)
// Matrix mat;			/* Holds the original and inverse */
// int actual_size;	    /* Actual size of matrix in use, (high_subscript+1)*/
{
    int i,j,k;
					/* Locations of pivot elements */
    int *pvt_i, *pvt_j;
    double pvt_val;                     /* Value of current pivot element */
    double hold;                        /* Temporary storage */
    double determ;                      /* Determinant */

    determ = 1.0;

    pvt_i = (int *) mxCalloc(actual_size, sizeof(int));
    pvt_j = (int *) mxCalloc(actual_size, sizeof(int));

    for (k = 0; k < actual_size; k++)
    {
        /* Locate k'th pivot element */
        pvt_val = mat[k][k];            /* Initialize for search */
        pvt_i[k] = k;
        pvt_j[k] = k;
        for (i = k; i < actual_size; i++)
          for (j = k; j < actual_size; j++)
            if (fabs(mat[i][j]) > fabs(pvt_val))
            {
                pvt_i[k] = i;
                pvt_j[k] = j;
                pvt_val = mat[i][j];
            }
        /* Product of pivots, gives determinant when finished */
        determ *= pvt_val;
        if (determ == 0.0) {    
         /* Matrix is singular (zero determinant). */
	    mxFree(pvt_i);
	    mxFree(pvt_j);
            return (0.0);              
	}

        /* "Interchange" rows (with sign change stuff) */
        i = pvt_i[k];
        if (i != k)                     /* If rows are different */
          for (j = 0; j < actual_size; j++)
          {
            hold = -mat[k][j];
            mat[k][j] = mat[i][j];
            mat[i][j] = hold;
          }

        /* "Interchange" columns */
        j = pvt_j[k];
        if (j != k)                     /* If columns are different */
          for (i = 0; i < actual_size; i++)
          {
            hold = -mat[i][k];
            mat[i][k] = mat[i][j];
            mat[i][j] = hold;
          }
        /* Divide column by minus pivot value */
        for (i = 0; i < actual_size; i++)
          if (i != k)                   /* Don't touch the pivot entry */
            mat[i][k] /= ( -pvt_val) ;  /* (Tricky C syntax for division) */

        /* Reduce the matrix */
        for (i = 0; i < actual_size; i++)
        {
            hold = mat[i][k];
            for (j = 0; j < actual_size; j++)
              if ( i != k && j != k )   /* Don't touch pivot. */
                mat[i][j] += hold * mat[k][j];
        }

        /* Divide row by pivot */
        for (j = 0; j < actual_size; j++)
          if (j != k)                   /* Don't touch the pivot! */
            mat[k][j] /= pvt_val;

        /* Replace pivot by reciprocal (at last we can touch it). */
        mat[k][k] = 1.0/pvt_val;
    }

    /* That was most of the work, one final pass of row/column interchange */
    /* to finish */
    for (k = actual_size-2; k >= 0; k--)  /* Don't need to work with 1 by 1 */
                                        /* corner */
    {
        i = pvt_j[k];		 /* Rows to swap correspond to pivot COLUMN */
        if (i != k)                     /* If rows are different */
          for(j = 0; j < actual_size; j++)
          {
            hold = mat[k][j];
            mat[k][j] = -mat[i][j];
            mat[i][j] = hold;
          }

        j = pvt_i[k];           /* Columns to swap correspond to pivot ROW */
        if (j != k)                     /* If columns are different */
          for (i = 0; i < actual_size; i++)
          {
            hold = mat[i][k];
            mat[i][k] = -mat[i][j];
            mat[i][j] = hold;
          }
    }

    mxFree(pvt_i);
    mxFree(pvt_j);
    return(determ);
}
      

