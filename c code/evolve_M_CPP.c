/*
  Magnetization evolution C function
  Version 1.0 16-July-2007 Noam Ben-Eliezer

  Features:
   (1) 1D sample
   (2) Time-independant 1D gradient
   (3) Time-independant 3D RF Field
   (4) Spatial-dependant 1D B0 field inhomogeneity

Syntax:
  [Mx,My,Mz] = evolve_M_CPP(stParams);

Input parameters
 M_initx      :  [none]   :   1D vector : x-coordinate of spins as a function of 'z' spatial location along sample
 M_inity      :  [none]   :   1D vector : y-coordinate of spins as a function of 'z' spatial location along sample
 M_initz      :  [none]   :   1D vector : z-coordinate of spins as a function of 'z' spatial location along sample
 B_eff_rotx   :  [T]      :   Scalar    : Bx
 B_eff_roty   :  [T]      :   Scalar    : By
 B_eff_rotz   :  [T]      :   Scalar    : Bz
 Tev          :  [sec]    :   Scalar    : Evolution duration
 dB0r         :  [T]      :   1D vector : Magnetic field inhomogeneity along the spatial axis
 rGrad        :  [G/cm]   :   Scalar    : Gz - Time independent magnetic field gradient
 r_axis       :  [cm]     :   1D vector : Spatial axis
 RH_flag      :  [none]   :   Scalar    : Right hand or left hand rotation
 inhomo_f     :  [none]   :   Scalar    : Enable / disable B0 inhomogeneity
 relax_flag   :  [none]   :   Scalar    : Enable / disable T1 & T2 Relaxation
 T1           :  [sec]    :   Scalar    : T1 relaxation constant
 T2           :  [sec]    :   Scalar    : T2 relaxation constant

Output parameters
 Mx           : [none]    :   1D vector : Post evolution x magnetization(z)
 My           : [none]    :   1D vector : Post evolution y magnetization(z)
 Mz           : [none]    :   1D vector : Post evolution z magnetization(z)
*/
#include "evolve.h"

double *unit_matrix(double *R_mat);
double *rot_RHR(double x, double y, double z, double phi, double *R_mat);
void print_R_mat(double *R);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double  gammaHz = 4.2574e+3,          // [Hz/G]
	        pi      = 3.1415926535897,
	        gamma_T = gammaHz*2*pi*1e+4,  // [rad/(sec*Tesla)] = [rad*amper*sec/kg]  ('H' protons)
	        Tev,
	        T1,T2,
	        B_grad_rt,
	        dz,
	        B_vec_sz,
	        B_eff_rotx,
	        B_eff_roty,
	        B_eff_rotz,
	        Bx, Bnorm_x,
	        By, Bnorm_y,
	        Bz, Bnorm_z,
	        Mx,My,Mz,
	        rGrad,
	        phi_rot;
	double *M_initx    = NULL,            // x-component of Input  Magnetization(z)
	       *M_inity    = NULL,            // y-component of Input  Magnetization(z)
	       *M_initz    = NULL,            // z-component of Input  Magnetization(z)
	       *M0z        = NULL,
	       *M_outx     = NULL,            // x-component of Output Magnetization(z)
	       *M_outy     = NULL,            // y-component of Output Magnetization(z)
	       *M_outz     = NULL,            // z-component of Output Magnetization(z)
	       *dB0r_in    = NULL,
	       *dB0r       = NULL,
	       *R_mat      = NULL,            // Rotation matrix. Rows are serially spreaded [(1,1)(1,2)(1,3)(2,1)...(3,2)(3,3)]
	       *r_axis     = NULL;
	int p_idx = 0,
	    t_idx,
	    r_idx,
	    n_r_steps,
	    Msz,
	    RH_flag,
	    inhomo_flag,
	    relax_flag,
	    n_struct,                             // number of field in structure
	    n_fields;                             // number of structures
	double Moutxsz,Moutysz,Moutzsz;
	int  ctxlen = 100;
	char context[100];

	/* -------------------------
	   Retrieve input parameters
	   ------------------------- */
	if (nrhs != 1) { 
		mexErrMsgTxt("One input required");
	}
	else if (nlhs != 3) {
		mexErrMsgTxt("Wrong number of output arguments");
	}

	n_fields    =         mxGetNumberOfFields(prhs[0]        );
	n_struct    =       mxGetNumberOfElements(prhs[0]        );
	// Magnetization
	M_initx     =  mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	M_inity     =  mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	M_initz     =  mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      // don't advance idx. Retrieve vector size.
	Msz         =  mxGetN (mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	M0z         =  mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;

	// RF Field
	B_eff_rotx  = *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	B_eff_roty  = *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	B_eff_rotz  = *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	// Evolution duration
	Tev         = *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	// Inhomogeneity
	dB0r_in     =  mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	// z-Gradient
	rGrad       = *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	// Spatial axis
	r_axis      =  mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      // don't advance idx. Retrieve vector size.
	n_r_steps   =  mxGetN (mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	// Flags
	RH_flag     = *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	inhomo_flag = *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	relax_flag  = *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	T1          = *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	T2          = *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));      p_idx++;
	mxGetString(mxGetFieldByNumber(prhs[0],0,p_idx),context,ctxlen);  p_idx++;

	// Basic validations (partial)
	if (!M_initx || !M_inity || !M_initz || !dB0r_in || !r_axis) {
		mexPrintf("Error parsing params: %p,%p,%p,%p,%p,%p,%p,%p,%p",M_initx,M_inity,M_initz,dB0r_in,r_axis);
		mexErrMsgTxt("Exiting...");
	}
	if (Msz != n_r_steps) {
		mexPrintf("Error parsing params Wrong M_init size (%d). Should be %d\n",Msz,n_r_steps);
		mexErrMsgTxt("Exiting...");
	}
	if (p_idx != n_fields) {
		mexPrintf("Error parsing params: Wrong number of struct fields (%d). Should be %d\n",p_idx,n_fields);
		mexErrMsgTxt("Exiting...");
	}
//	mexPrintf("Input params:\n n_fields=%d\n n_struct=%d\n B_eff_rotx=%f\n B_eff_roty=%f\n B_eff_rotz=%f\n rGrad=%f\n Tev=%f\n n_r_steps=%d\n RH_flag=%d\n inhomo_flag=%d\n relax_flag=%d\n T1=%f\n T2=%f\n",
//		      n_fields,n_struct,B_eff_rotx,B_eff_roty,B_eff_rotz,rGrad,Tev,n_r_steps,RH_flag,inhomo_flag,relax_flag,T1,T2);
//	mexPrintf("Input pointers:\n %p\n %p\n %p\n %p\n %p\n",M_initx,M_inity,M_initz,dB0r,r_axis);

	/* ------------------------
	   Create output parameters 
	   ------------------------ */
	plhs[0] = mxCreateDoubleMatrix(1,n_r_steps, mxREAL);   M_outx = mxGetPr(plhs[0]);  Moutxsz = mxGetN(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(1,n_r_steps, mxREAL);   M_outy = mxGetPr(plhs[1]);  Moutysz = mxGetN(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(1,n_r_steps, mxREAL);   M_outz = mxGetPr(plhs[2]);  Moutzsz = mxGetN(plhs[2]);
//	mexPrintf("Output:\n Moutxsz=%f\n Moutysz=%f\n Moutzsz=%f\n",Moutxsz,Moutysz,Moutzsz);
	
	/* ---------------
	   Initializations
	   --------------- */
	dz = r_axis[1] - r_axis[0];
	R_mat = (double *)calloc(9,sizeof(double));
	inhomo_flag = inhomo_flag? 1:0;                            // Flag must be either 0 or 1

	/* -----------------------------
	    For faster run time ...
	   -----------------------------*/
	dB0r = (double *)calloc(n_r_steps,sizeof(double));
	for (r_idx = 0; r_idx < n_r_steps; r_idx++) {
		dB0r[r_idx]  = dB0r_in[r_idx]*inhomo_flag;
	}

	/* -----------------------
	   Magnetization evolution
	   ----------------------- */
	Bx = B_eff_rotx;
	By = B_eff_roty;
//	mexPrintf(" gamma_T = %f\n Tev = %f\n",gamma_T,Tev);
	for (r_idx = 0; r_idx < n_r_steps; r_idx++) {              // Spatial forloop
		B_grad_rt = rGrad * 1E-4 * (r_axis[r_idx]);            // Gz converted to [T/cm]
		Bz = B_eff_rotz + B_grad_rt + dB0r[r_idx];
		B_vec_sz  = sqrt(Bx*Bx + By*By + Bz*Bz);
		if (B_vec_sz == 0)
			R_mat = unit_matrix(R_mat);                        // No magnetic field --> no precession
		else {
			Bnorm_x = Bx / B_vec_sz;
			Bnorm_y = By / B_vec_sz;
			Bnorm_z = Bz / B_vec_sz;
			phi_rot = gamma_T*B_vec_sz*Tev;                    // [rad] = [rad/(sec*T)] * [T] * [sec]

			if (RH_flag)
				R_mat = rot_RHR(Bnorm_x, Bnorm_y, Bnorm_z, phi_rot, R_mat);
            else
				R_mat = rot_LHR(Bnorm_x, Bnorm_y, Bnorm_z, phi_rot, R_mat);
		}
//		if ((r_idx % 50) == 0)
//			mexPrintf("Bx=%f  By=%f  Bz=%f  B_grad_rt=%f  B_vec_sz=%f Phi=%f\n",Bx,By,Bz,B_grad_rt,B_vec_sz,phi_rot);	
//		mexPrintf("B_grad_rt=%f  dB0r[r_idx]=%f  Bz=%f  B_vec_sz=%f Phi=%f\n",B_grad_rt,dB0r[r_idx],Bz,B_vec_sz,phi_rot);
		Mx = M_initx[r_idx];
		My = M_inity[r_idx];
		Mz = M_initz[r_idx];
		M_outx[r_idx] = R_mat[0]*Mx + R_mat[1]*My + R_mat[2]*Mz;
		M_outy[r_idx] = R_mat[3]*Mx + R_mat[4]*My + R_mat[5]*Mz;
		M_outz[r_idx] = R_mat[6]*Mx + R_mat[7]*My + R_mat[8]*Mz;

 		if (relax_flag) {
			T1_relaxation(&(M_outx[r_idx]), &(M_outy[r_idx]), &(M_outz[r_idx]), M0z[r_idx], Tev, T1);
			T2_relaxation(&(M_outx[r_idx]), &(M_outy[r_idx]), &(M_outz[r_idx]), Tev, T2);
 		}

	}

	free(R_mat);
	free(dB0r);
	
	return;
}

