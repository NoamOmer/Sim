/*
  2D (XY) Right-handed magnetization evolution in the presence of a 2D gradient
  Version 1.0 24-February-2008 Noam Ben-Eliezer

  Features:
   (1) 2D sample
   (2) Time independent 2D gradient
   (3) Time independent 3D RF Field
   (4) Spatial-dependent 2D B0 field inhomogeneity
   (5) No signal acqusition

Syntax:
 [Mx,My,Mz] = evolve_M_n_acq_2D_tIndependent_C(ParamsSt);

Input parameters
 M_initx      :        : 1D vector : x-coordinate of spins as a function of (x,y) location in sample
 M_inity      :        : 1D vector : y-coordinate of spins as a function of (x,y) location in sample
 M_initz      :        : 1D vector : z-coordinate of spins as a function of (x,y) location in sample
 B_eff_rot    : [T]    : 1D vector : [Bx,By,Bz]
 dB0_XY       : [T]    : 1D vector : Magnetic field inhomogeneity as a function of (x,y) location in sample
 zGrad        : [G/cm] : 1D vector : [Gx,Gy,Gz] magnetic field gradient
 x_axis       : [cm]   : 1D vector : X-axis
 y_axis       : [cm]   : 1D vector : Y-axis
 Tevolution   : [sec]  : Scalar    : Evolution duration
 inhomo_f     : [none] : Scalar    : Enable / disable B0 inhomogeneity
 relax_flag   : [none] : Scalar    : Enable / disable T1 & T2 Relaxation
 plot_flag    : [none] : Scalar    : Enable / disable plotting
 T1           : [sec]  : Scalar    : T1 relaxation constant
 T2           : [sec]  : Scalar    : T2 relaxation constant
 context      :        : String    : Context name

Output parameters
 Mx           :        : 1D vector : Post evolution x magnetization(z)
 My           :        : 1D vector : Post evolution y magnetization(z)
 Mz           :        : 1D vector : Post evolution z magnetization(z)
*/

#include "../C code/evolve.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double  gammaHz = 4.2574e+3,          // [Hz/G]
	        pi      = 3.1415926535897,
	        gamma_T = gammaHz*2*pi*1e+4,  // [rad/(sec*Tesla)] = [rad*amper*sec/kg]  ('H' protons)
	        Tevolution,
	        T1,T2,
	        B_zGrad,
	        B_zGrad_y,
	        dx,dy,
	        B_vec_sz,
	        Bx, Bnorm_x,
	        By, Bnorm_y,
	        Bz, Bnorm_z,
	        Mx,My,Mz,
	        phi_rot;
	double *M_initx    = NULL,            // x-component of Input  Magnetization(z)
	       *M_inity    = NULL,            // y-component of Input  Magnetization(z)
	       *M_initz    = NULL,            // z-component of Input  Magnetization(z)
	       *M0z        = NULL,
	       *M_outx     = NULL,            // x-component of Output Magnetization(z)
	       *M_outy     = NULL,            // y-component of Output Magnetization(z)
	       *M_outz     = NULL,            // z-component of Output Magnetization(z)
	       *B_eff_rot  = NULL,
	       *dB0_XY     = NULL,
	       *zGrad      = NULL,
	       *R_mat      = NULL,            // Rotation matrix. Rows are serially spreaded [(1,1)(1,2)(1,3)(2,1)...(3,2)(3,3)]
	       *x_axis     = NULL,
	       *y_axis     = NULL;
	Matrix M_initx_mat = NULL,
	       M_inity_mat = NULL,
	       M_initz_mat = NULL,
	       M_outx_mat  = NULL,
	       M_outy_mat  = NULL,
	       M_outz_mat  = NULL,
	       dB0_XY_mat  = NULL;
	int p_idx = 0,
	    t_idx,
	    x_idx,
	    y_idx,
	    cols_idx,
	    rows_idx,
	    Msz,
	    Nx,Ny,
	    RH_flag = 1,                      // Use only RH rotation
	    inhomo_flag,
	    relax_flag,
	    plot_flag,
	    n_struct,                         // number of field in structure
	    n_fields;                         // number of structures

	int  ctxlen = 100;
	char context[100];

fprintf(stdout,"[nbe] 2014_06_08 \n");
fprintf(stdout,"Before running this function check the z-grad an inhomogeneity parameters\n");
fprintf(stdout,"(and possibly other parameters) to see if the MATLAB values are overwritten\n");
fprintf(stdout,"by the C code. Merge improved code from: evolve_M_n_acq1Dgrad_t_CPP.c\n");



	/* -------------------------
	   Retrieve input parameters
	   ------------------------- */
	if (nrhs != 1) {
		mexErrMsgTxt("One input required");                  // See API header at the top of this page
	}
	else if (nlhs != 3) {
		mexErrMsgTxt("Five output arguments required");      // See API header at the top of this page
	}

	n_fields    =                  mxGetNumberOfFields(prhs[0]        );
	n_struct    =                mxGetNumberOfElements(prhs[0]        );
	// Magnetization
	plhs[0]     =  mxDuplicateArray(mxGetFieldByNumber(prhs[0],0,p_idx));   // Initialize M_outx with M_inx
	M_initx     =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	plhs[1]     =  mxDuplicateArray(mxGetFieldByNumber(prhs[0],0,p_idx));   // Initialize M_outy with M_iny
	M_inity     =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	plhs[2]     =  mxDuplicateArray(mxGetFieldByNumber(prhs[0],0,p_idx));   // Initialize M_outz with M_inz
	M_initz     =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   // don't advance idx. Retrieve vector size.
	Msz         =           mxGetN (mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	M0z         =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;

	// RF Field
	B_eff_rot   =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	// Inhomogeneity
	dB0_XY      =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	// z-Gradient
	zGrad       =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	// Spatial axis
	x_axis      =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   // don't advance idx. Retrieve vector size.
	Nx          =           mxGetN (mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	y_axis      =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   // don't advance idx. Retrieve vector size.
	Ny          =           mxGetN (mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	// Temporal resolution
	Tevolution  =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	// Flags
	inhomo_flag =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	relax_flag  =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	plot_flag   =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	T1          =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	T2          =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	mxGetString(mxGetFieldByNumber(prhs[0],0,p_idx),context,ctxlen);        p_idx++;

	// Basic validations (partial)
	if (!M_initx || !M_inity || !M_initz || !B_eff_rot || !dB0_XY || !zGrad || !x_axis || !y_axis) {
		mexPrintf("Error parsing params: %p,%p,%p,%p,%p,%p,%p,%p,%p,%p",M_initx,M_inity,M_initz,B_eff_rot,dB0_XY,zGrad,x_axis,y_axis);
		mexErrMsgTxt("Exiting...");
	}
	if (Msz != (Nx*Ny)) {
		mexPrintf("Error parsing params Wrong M_init size (%d). Should be %d*%d\n",Msz,Nx,Ny);
		mexErrMsgTxt("Exiting...");
	}
	if (p_idx != n_fields) {
		mexPrintf("Error parsing params: Wrong number of struct fields (%d). Should be %d\n",p_idx,n_fields);
		mexErrMsgTxt("Exiting...");
	}
	dx = x_axis[1] - x_axis[0];
	dy = y_axis[1] - y_axis[0];

//	mexPrintf("Input params:\n n_fields=%d\n n_struct=%d\n Msz=%d\n Tevolution=%f\n Bx=%f\n By=%f\n Bz=%f\n Gx=%f\n Gy=%f\n Gz=%f\n dx=%f\n dy=%f\n Nx=%d\n Ny=%d\n inhomo_flag=%d\n relax_flag=%d\n plot_flag=%d\n T1=%f\n T2=%f\n",
//	           n_fields,n_struct,Msz,Tevolution,B_eff_rot[0],B_eff_rot[1],B_eff_rot[2],zGrad[0],zGrad[1],zGrad[2],dx,dy,Nx,Ny,inhomo_flag,relax_flag,plot_flag,T1,T2);
//	mexPrintf("Input pointers:\n %p\n %p\n %p\n %p\n %p\n %p\n %p\n %p\n",M_initx,M_inity,M_initz,B_eff_rot,dB0_XY,zGrad,x_axis,y_axis);

	/* ------------------------
	   Create output parameters 
	   ------------------------ */
	M_outx   = mxGetPr(plhs[0]);
	M_outy   = mxGetPr(plhs[1]);
	M_outz   = mxGetPr(plhs[2]);

	/* ---------------
	   Initializations
	   --------------- */
	// Allocate space for rotation matrix
	R_mat = (double *)mxCalloc(9,sizeof(double));

	// Verify inhomogeneity flag value
	inhomo_flag = inhomo_flag? 1:0;                                 // Flag must be either 0 or 1

	// Allocate space for input & output magnetization matrices and for inhomogeneity matrix
	M_initx_mat = NewMatrix(Nx, Ny);
	M_inity_mat = NewMatrix(Nx, Ny);
	M_initz_mat = NewMatrix(Nx, Ny);
	M_outx_mat  = NewMatrix(Nx, Ny);
	M_outy_mat  = NewMatrix(Nx, Ny);
	M_outz_mat  = NewMatrix(Nx, Ny);
	dB0_XY_mat  = NewMatrix(Nx, Ny);

	// Reshape input magnetization, output magnetization & inhomogeneity vectors
	for (cols_idx = 0; cols_idx < Ny; cols_idx++) {
		for (rows_idx = 0; rows_idx < Nx; rows_idx++) {
			M_initx_mat[rows_idx][cols_idx] = M_initx[cols_idx*Nx + rows_idx];
			M_inity_mat[rows_idx][cols_idx] = M_inity[cols_idx*Nx + rows_idx];
			M_initz_mat[rows_idx][cols_idx] = M_initz[cols_idx*Nx + rows_idx];
			M_outx_mat [rows_idx][cols_idx] = M_outx [cols_idx*Nx + rows_idx];
			M_outy_mat [rows_idx][cols_idx] = M_outy [cols_idx*Nx + rows_idx];
			M_outz_mat [rows_idx][cols_idx] = M_outz [cols_idx*Nx + rows_idx];
			dB0_XY_mat [rows_idx][cols_idx] = dB0_XY [cols_idx*Nx + rows_idx];
		}
	}
	/* -----------------------
	   Magnetization evolution
	   ----------------------- */
	if(plot_flag)
		print_2D_M_mat(M_outx_mat,M_outy_mat,M_outz_mat,Nx,Ny,0,context); // once before the evolution

	Bx = B_eff_rot[0];
	By = B_eff_rot[1];

	for (y_idx = 0; y_idx < Ny; y_idx++) {                  // X-Spatial forloop
//		mexPrintf("    Spatial forloop #1: x_idx=%d\n",x_idx);

		B_zGrad_y = (zGrad[1]) * 1E-4 * (y_axis[y_idx]);                     // [T] = Gy[G/cm] * 1E-4 * [cm]
		for (x_idx = 0; x_idx < Nx; x_idx++) {              // Y-Spatial forloop
//			mexPrintf("        Spatial forloop #2: y_idx=%d\n",y_idx);
			B_zGrad  = (zGrad[0]) * 1E-4 * (x_axis[x_idx]) + B_zGrad_y;      // [T] = Gx[G/cm] * 1E-4 * [cm]
//			B_zGrad += (zGrad[1]) * 1E-4 * (y_axis[y_idx]);                  Moved this line to outer forloop
			
			Bz  = B_eff_rot[2] + B_zGrad + dB0_XY_mat[x_idx][y_idx] * inhomo_flag;
			B_vec_sz  = sqrt(Bx*Bx + By*By + Bz*Bz);

			if (B_vec_sz == 0)
				R_mat = unit_matrix(R_mat);                        // No magnetic field --> no precession
			else {
				Bnorm_x = Bx / B_vec_sz;
				Bnorm_y = By / B_vec_sz;
				Bnorm_z = Bz / B_vec_sz;
				phi_rot = gamma_T*B_vec_sz*Tevolution;             // [rad] = [rad/(sec*T)] * [T] * [sec]
				R_mat = rot_RHR(Bnorm_x, Bnorm_y, Bnorm_z, phi_rot, R_mat);
			}

			Mx = M_outx_mat[x_idx][y_idx];
			My = M_outy_mat[x_idx][y_idx];
			Mz = M_outz_mat[x_idx][y_idx];
			M_outx_mat[x_idx][y_idx] = R_mat[0]*Mx + R_mat[1]*My + R_mat[2]*Mz;
			M_outy_mat[x_idx][y_idx] = R_mat[3]*Mx + R_mat[4]*My + R_mat[5]*Mz;
			M_outz_mat[x_idx][y_idx] = R_mat[6]*Mx + R_mat[7]*My + R_mat[8]*Mz;

			if (relax_flag) {
				T1_relaxation(&(M_outx_mat[x_idx][y_idx]),
					          &(M_outy_mat[x_idx][y_idx]),
					          &(M_outz_mat[x_idx][y_idx]),
					          M_initz_mat[x_idx][y_idx]  , Tevolution, T1);  Fix this to use M0z instead of Minitz_mat... which is a BUG!
				T2_relaxation(&(M_outx_mat[x_idx][y_idx]),
					          &(M_outy_mat[x_idx][y_idx]),
					          &(M_outz_mat[x_idx][y_idx]), Tevolution, T2);
	 		}

		}
	}		

	// Reshape the output magnetization 2D matrices into 1D vectors
	for (cols_idx = 0; cols_idx < Ny; cols_idx++) {
		for (rows_idx = 0; rows_idx < Nx; rows_idx++) {
			M_outx[cols_idx*Nx + rows_idx] = M_outx_mat[rows_idx][cols_idx];
			M_outy[cols_idx*Nx + rows_idx] = M_outy_mat[rows_idx][cols_idx];
			M_outz[cols_idx*Nx + rows_idx] = M_outz_mat[rows_idx][cols_idx];
		}
	}
//	mexPrintf("END Reshape\n\n");

	// Free allocated memory
	FreeMatrix(M_initx_mat, Nx);
	FreeMatrix(M_inity_mat, Nx);
	FreeMatrix(M_initz_mat, Nx);
	FreeMatrix(M_outx_mat , Nx);
	FreeMatrix(M_outy_mat , Nx);
	FreeMatrix(M_outz_mat , Nx);
	FreeMatrix(dB0_XY_mat , Nx);
	mxFree(R_mat);
//	mexPrintf("END Free\n\n");
	
	return;
}

