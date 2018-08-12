/*
  Magnetization evolution C function
  Version 1.0 16-July-2007 Noam Ben-Eliezer

  Features:
   (1) 1D sample
   (2) Time-dependant 1D gradient
   (3) Time-dependant 3D RF Field
   (4) Spatial-dependant 1D B0 field inhomogeneity
   (5) Signal acqusition

Syntax:
 [Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(ParamsSt);
 
Input parameters
 M_initx      :  [none]   :   1D vector : x-coordinate of spins as a function of 'z' spatial location along sample
 M_inity      :  [none]   :   1D vector : y-coordinate of spins as a function of 'z' spatial location along sample
 M_initz      :  [none]   :   1D vector : z-coordinate of spins as a function of 'z' spatial location along sample
 B_eff_rotx   :  [T]      :   1D vector : Bx(t)
 B_eff_roty   :  [T]      :   1D vector : By(t)
 B_eff_rotz   :  [T]      :   1D vector : Bz(t)
 dB0z         :  [T]      :   1D vector : dB0z(z) Magnetic field inhomogeneity along the spatial axis
 zGrad        :  [G/cm]   :   1D vector : Gz(t)   Time dependent magnetic field gradient in z-direction.
 z_axis       :  [cm]     :   1D vector : Spatial axis
 RH_flag      :  [none]   :   Scalar    : Right hand or left hand rotation
 dt           :  [sec]    :   Scalar    : Time axis interval
 acquire_flag :  [none]   :   Scalar    : Enable / disable signal acquisition
 inhomo_f     :  [none]   :   Scalar    : Enable / disable B0 inhomogeneity
 relax_flag   :  [none]   :   Scalar    : Enable / disable T1 & T2 Relaxation
 T1           :  [sec]    :   Scalar    : T1 relaxation constant
 T2           :  [sec]    :   Scalar    : T2 relaxation constant

Output parameters
 Mx           : [none]    :   1D vector : Post evolution x magnetization(z)
 My           : [none]    :   1D vector : Post evolution y magnetization(z)
 Mz           : [none]    :   1D vector : Post evolution z magnetization(z)
 Sig_real     : []        :   1D vector : Real part of MRI Signal Vs. time
 Sig_imag     : []        :   1D vector : Imag part of MRI Signal Vs. time
*/
#include "evolve.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double  gammaHz = 4.2574e+3,          // [Hz/G]
	        pi      = 3.1415926535897,
	        gamma_T = gammaHz*2*pi*1e+4,  // [rad/(sec*Tesla)] = [rad*amper*sec/kg]  ('H' protons)
	        dt,
	        T1,T2,
	        B_grad_zt,
	        dz,
	        B_vec_sz,
	        Bx, Bnorm_x,
	        By, Bnorm_y,
	        Bz, Bnorm_z,
	        Mx,My,Mz,
	        phi_rot,
	        gamma_dt,
	        s_real,
	        s_imag;
	double *M_initx    = NULL,            // x-component of Input  Magnetization(z)
	       *M_inity    = NULL,            // y-component of Input  Magnetization(z)
	       *M_initz    = NULL,            // z-component of Input  Magnetization(z)
	       *M0z        = NULL,
	       *M_outx     = NULL,            // x-component of Output Magnetization(z)
	       *M_outy     = NULL,            // y-component of Output Magnetization(z)
	       *M_outz     = NULL,            // z-component of Output Magnetization(z)
	       *B_eff_rotx = NULL,
	       *B_eff_roty = NULL,
	       *B_eff_rotz = NULL,
	       *dB0z       = NULL,
	       *zGrad      = NULL,
	       *R_mat      = NULL,            // Rotation matrix. Rows are serially spreaded [(1,1)(1,2)(1,3)(2,1)...(3,2)(3,3)]
	       *R_mat_unit = NULL,
	       *R_mat_ptr  = NULL,
	       *z_axis     = NULL,
	       *Sig_real   = NULL,
	       *Sig_imag   = NULL,
	       *acc_s_real = NULL,
	       *acc_s_imag = NULL;
	int p_idx = 0,
	    t_idx,
	    n_t_steps,
	    r_idx,
	    n_r_steps,
	    Msz,
	    RH_flag,
	    acq_flag,
	    inhomo_flag,
	    relax_flag,
	    n_struct,                         // number of field in structure
	    n_fields,                         // number of structures
	    n_jumps = 100;                    // # of time epochs between debug points
	double Moutxsz,Moutysz,Moutzsz,Sigresz,Sigimsz;

	int  ctxlen = 100;
	char context[100];

#define acc_s_sz 10
	mxArray *acc_s_mxArr[acc_s_sz];       // Array for internal debugging information
	
	/* -------------------------
	   Retrieve input parameters
	   ------------------------- */
	if (nrhs != 1) {
		mexErrMsgTxt("One input required");
	}
	else if (nlhs != 5) {
		mexErrMsgTxt("Five output arguments required");
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
	B_eff_rotx  =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	B_eff_roty  =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	B_eff_rotz  =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   // don't advance idx. Retrieve vector size.
	// # of Temporal steps
	n_t_steps   =           mxGetN (mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	// Inhomogeneity
	dB0z        =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	// z-Gradient
	zGrad       =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	// Spatial axis
	z_axis      =           mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   // don't advance idx. Retrieve vector size.
	n_r_steps   =           mxGetN (mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	// Temporal resolution
	dt          =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	// Flags
	RH_flag     =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	acq_flag    =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	inhomo_flag =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	relax_flag  =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	T1          =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	T2          =          *mxGetPr(mxGetFieldByNumber(prhs[0],0,p_idx));   p_idx++;
	mxGetString(mxGetFieldByNumber(prhs[0],0,p_idx),context,ctxlen);        p_idx++;

	// Basic validations (partial)
	if (!M_initx || !M_inity || !M_initz || !B_eff_rotx || !B_eff_roty || !B_eff_rotz || !dB0z || !zGrad || !z_axis) {
		mexPrintf("Error parsing params: %p,%p,%p,%p,%p,%p,%p,%p,%p",M_initx,M_inity,M_initz,B_eff_rotx,B_eff_roty,B_eff_rotz,dB0z,zGrad,z_axis);
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
	dz = z_axis[1] - z_axis[0];
//	mexPrintf("Input params:\n n_fields=%d\n n_struct=%d\n n_t_steps=%d\n dt=%f\n n_r_steps=%d\n dz=%f\n RH_flag=%d\n acq_flag=%d\n inhomo_flag=%d\n relax_flag=%d\n T1=%f\n T2=%f\n",
//		      n_fields,n_struct,n_t_steps,dt,n_r_steps,dz,RH_flag,acq_flag,inhomo_flag,relax_flag,T1,T2);
//	mexPrintf("Input pointers:\n %p\n %p\n %p\n %p\n %p\n %p\n %p\n %p\n %p\n",M_initx,M_inity,M_initz,B_eff_rotx,B_eff_roty,B_eff_rotz,dB0z,zGrad,z_axis);
	
	/* ------------------------
	   Create output parameters 
	   ------------------------ */
	plhs[3]  = mxCreateDoubleMatrix(1,n_t_steps, mxREAL);
	plhs[4]  = mxCreateDoubleMatrix(1,n_t_steps, mxREAL);
	M_outx   = mxGetPr(plhs[0]);  Moutxsz = mxGetN(plhs[0]);
	M_outy   = mxGetPr(plhs[1]);  Moutysz = mxGetN(plhs[1]);
	M_outz   = mxGetPr(plhs[2]);  Moutzsz = mxGetN(plhs[2]);
	Sig_real = mxGetPr(plhs[3]);  Sigresz = mxGetN(plhs[3]);
	Sig_imag = mxGetPr(plhs[4]);  Sigimsz = mxGetN(plhs[4]);
//	mexPrintf("Output:\n Moutxsz=%f\n Moutysz=%f\n Moutzsz=%f\n Sigresz=%f\n Sigimsz=%f\n",Moutxsz,Moutysz,Moutzsz,Sigresz,Sigimsz);
	
	/* ---------------
	   Initializations
	   --------------- */
//	acc_s_mxArr[0] = mxCreateDoubleMatrix(1,n_r_steps+1, mxREAL);   // Real part of accumulated signal Vs. z (partially accumulated over part of the sample)
//	acc_s_mxArr[1] = mxCreateDoubleMatrix(1,n_r_steps+1, mxREAL);   // Imag part of accumulated signal Vs. z
//	acc_s_mxArr[2] = plhs[0];                                       // M_outx (in mid-evolution)
//	acc_s_mxArr[3] = plhs[1];                                       // M_outy (in mid-evolution)
//	acc_s_mxArr[4] = plhs[2];                                       // M_outz (in mid-evolution)
//	acc_s_mxArr[7] = mxCreateDoubleScalar(dt);                      // dt
//	acc_s_mxArr[8] = mxCreateDoubleScalar(n_jumps);                 // Number of debug points
//	acc_s_mxArr[9] = mxCreateString(context);                       // context
//	acc_s_real = mxGetPr(acc_s_mxArr[0]);    acc_s_real[0] = 0;
//	acc_s_imag = mxGetPr(acc_s_mxArr[1]);    acc_s_imag[0] = 0;
	
	R_mat      = (double *)calloc(9,sizeof(double));
	R_mat_unit = (double *)calloc(9,sizeof(double));
	R_mat_unit = unit_matrix(R_mat_unit);                        // No magnetic field --> no precession

	inhomo_flag = inhomo_flag? 1:0;                                 // Flag must be either 0 or 1

	/* -----------------------------
	    For faster run time ...
	   -----------------------------*/
	gamma_dt = gamma_T*dt;
	for (r_idx = 0; r_idx < n_r_steps; r_idx++) {
		dB0z[r_idx]  = dB0z[r_idx]*inhomo_flag;
	}
	for (t_idx = 0; t_idx < n_t_steps; t_idx++) {
		zGrad[t_idx] = zGrad[t_idx] * 1E-4;
	}

	/* -----------------------
	   Magnetization evolution
	   ----------------------- */
	for (t_idx = 0; t_idx < n_t_steps; t_idx++) {                   // Temporal forloop
		s_real = 0;
		s_imag = 0;

		Bx = B_eff_rotx[t_idx];
		By = B_eff_roty[t_idx];

		// Keep the current magnetization in order to compare it with its value after this time epoch
//		acc_s_mxArr[5] = mxDuplicateArray(acc_s_mxArr[2]);
//		acc_s_mxArr[6] = mxDuplicateArray(acc_s_mxArr[3]);
		for (r_idx = 0; r_idx < n_r_steps; r_idx++) {                       // Spatial forloop
//			B_grad_zt = (zGrad[t_idx]) * 1E-4 * (z_axis[r_idx]);            // Gz converted to [T]
			B_grad_zt = (zGrad[t_idx])        * (z_axis[r_idx]);            //                    ... done before forloop

//			Bz = B_eff_rotz[t_idx] + B_grad_zt + dB0z[r_idx]*inhomo_flag;
			Bz = B_eff_rotz[t_idx] + B_grad_zt + dB0z[r_idx];               // dB0z is multiplied by inhomo_flag before forloop

			B_vec_sz  = sqrt(Bx*Bx + By*By + Bz*Bz);

			if (B_vec_sz == 0)
			{ // (3) even faster: comment all ... there's really nothing to be done here -- M_outx,y,z remain unchanged
			
				// (2) faster: we know which elements of the R_mat_unit are =1 and which are =0
//				M_outx[r_idx] =               Mx;
//				M_outy[r_idx] =                                  My;
//				M_outz[r_idx] =                                                     Mz;

				// (1) original
//				R_mat = unit_matrix(R_mat);                        // No magnetic field --> no precession
//				Mx = M_outx[r_idx];
//				My = M_outy[r_idx];
//				Mz = M_outz[r_idx];
				
//				M_outx[r_idx] = R_mat_unit[0]*Mx + R_mat_unit[1]*My + R_mat_unit[2]*Mz;
//				M_outy[r_idx] = R_mat_unit[3]*Mx + R_mat_unit[4]*My + R_mat_unit[5]*Mz;
//				M_outz[r_idx] = R_mat_unit[6]*Mx + R_mat_unit[7]*My + R_mat_unit[8]*Mz;
			}
			else
			{
				Bnorm_x = Bx / B_vec_sz;
				Bnorm_y = By / B_vec_sz;
				Bnorm_z = Bz / B_vec_sz;
//				phi_rot = gamma_T*B_vec_sz*dt;                     // [rad] = [rad/(sec*T)] * [T] * [sec]
				phi_rot = gamma_dt*B_vec_sz;                       // [rad] = [rad/(sec*T)] * [T] * [sec]
				if (RH_flag)
					R_mat = rot_RHR(Bnorm_x, Bnorm_y, Bnorm_z, phi_rot, R_mat);
	            else
					R_mat = rot_LHR(Bnorm_x, Bnorm_y, Bnorm_z, phi_rot, R_mat);

//				mexPrintf(" B_grad_zt=%f  dB0z[r_idx]=%f  Bx=%f  By=%f  Bz=%f  B_vec_sz=%f\n",B_grad_zt,dB0z[r_idx],Bx,By,Bz,B_vec_sz);
//				mexPrintf("Phi = %f\n",phi_rot);
//				mexPrintf("rot_f=%f\n",phi_rot);
//				mexPrintf("dt = %f\n",dt);
				Mx = M_outx[r_idx];
				My = M_outy[r_idx];
				Mz = M_outz[r_idx];
				M_outx[r_idx] = R_mat[0]*Mx + R_mat[1]*My + R_mat[2]*Mz;
				M_outy[r_idx] = R_mat[3]*Mx + R_mat[4]*My + R_mat[5]*Mz;
				M_outz[r_idx] = R_mat[6]*Mx + R_mat[7]*My + R_mat[8]*Mz;
			}

			if (acq_flag) {
				s_real += M_outx[r_idx];
				s_imag += M_outy[r_idx];
			}			
//			acc_s_real[r_idx+1] = acc_s_real[r_idx] + M_outx[r_idx];
//			acc_s_imag[r_idx+1] = acc_s_imag[r_idx] + M_outy[r_idx];

			if (relax_flag) {
				T1_relaxation(&(M_outx[r_idx]), &(M_outy[r_idx]), &(M_outz[r_idx]), M0z[r_idx], dt, T1);
				T2_relaxation(&(M_outx[r_idx]), &(M_outy[r_idx]), &(M_outz[r_idx]), dt, T2);
	 		}
		}
		if (acq_flag) {
			Sig_real[t_idx] = s_real * dz;
			Sig_imag[t_idx] = s_imag * dz;

//			if ((t_idx % (int)(floor(n_t_steps / n_jumps))) == 0) {
//				mexPrintf(" t_idx=%f\n n_t_steps=%f\n n_jumps=%f\n",t_idx,n_t_steps,n_jumps);
//				mexCallMATLAB(0, NULL, acc_s_sz, acc_s_mxArr, "plot_accumulated_acq");
//			}
		}

		// Plot mid-pulse spatial profile
//		if ((t_idx % (int)(floor(n_t_steps / n_jumps))) == 0) {
//			mexPrintf(" t_idx=%f\n n_t_steps=%f\n n_jumps=%f\n",t_idx,n_t_steps,n_jumps);
//			mexCallMATLAB(0, NULL, acc_s_sz, acc_s_mxArr, "plot_accumulated_acq");
//		}

	}

	free(R_mat);
	free(R_mat_unit);
	return;
}

