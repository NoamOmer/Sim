#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: spinsout = ApplyPulse(spinsin, pulse)
// spins structure:
//      spins(i).r      3x1 double      mm
//      spins(i).M      3x1 double      a.u.
//      spins(i).cs     1x1 double      kHz
// pulse structure:
//      pulse.tp        1x1 double      mm
//      pulse.RFamp     1xNp double     kHz
//      pulse.RFphase   1xNp double     radians
//      pulse.Gx        1xNp double     kHz/mm
//      pulse.Gy        1xNp double     kHz/mm
//      pulse.Gz        1xNp double     kHz/mm
{
    int i,j,current_spin;                                  // counter variables
    double vx, vy, vz;                                     // temporary vector components
    double dt;                                             // Pulse time step in ms
    double *RFamp, *RFphase, *Gx, *Gy, *Gz;                // pointers to RF data
    double tp;                                             // RF pulse duration in ms
    int NSpins, NSteps;                                    // Number of spins and of pulse-steps
    double gm, cs;
    double *r;
    double *Mout;
    double *Bx, *By, *Bz;
    double norm, angle, ct, st, nx, ny, nz, dotproduct;    // Internal variables for bloch equations
    
    // Gyromagnetic ratio of Proton in kHz/Gauss/10
    gm = 4.257/10;
    
    // Retrieve Pulse data
    tp = *mxGetPr(mxGetFieldByNumber(prhs[1],0,0));
    RFamp = mxGetPr(mxGetFieldByNumber(prhs[1],0,1));
    RFphase = mxGetPr(mxGetFieldByNumber(prhs[1],0,2));
    Gx = mxGetPr(mxGetFieldByNumber(prhs[1],0,3));
    Gy = mxGetPr(mxGetFieldByNumber(prhs[1],0,4));
    Gz = mxGetPr(mxGetFieldByNumber(prhs[1],0,5));
    NSteps = mxGetN(mxGetFieldByNumber(prhs[1],0,1));
    
    // Time step, in ms
    dt = tp/NSteps;
    
    // Initialize arrays to contain effective B's x, y, z components (in kHz)
    Bx = (double*) malloc (sizeof(double)*NSteps);
    By = (double*) malloc (sizeof(double)*NSteps);
    Bz = (double*) malloc (sizeof(double)*NSteps);
    for (i=0; i<NSteps; i++)
    {
        Bx[i]=RFamp[i]*cos(RFphase[i]);
        By[i]=RFamp[i]*sin(RFphase[i]);
    }
    
    // Retrieve number of spins
    NSpins = mxGetNumberOfElements(prhs[0]);
    
    // Create output structure (identical to input structure)
    plhs[0] = mxDuplicateArray(prhs[0]);
    
    for (current_spin=0; current_spin<NSpins; current_spin++)          // Loop over spins
    {
        r = mxGetPr(mxGetFieldByNumber(prhs[0], current_spin, 0));     // Get pointer to vector of spin location
        Mout = mxGetPr(mxGetFieldByNumber(plhs[0], current_spin, 1));  // Get pointer to output magnetization of spin object
        vx = Mout[0]; vy = Mout[1]; vz = Mout[2];                      // Store (temporarily) x,y,z components of output magnetization
        cs = mxGetScalar(mxGetFieldByNumber(prhs[0],current_spin,2));  // Get chemical shift of current spin (in kHz)
        // Calculate Bz in kHz as a function of time
        for (i=0; i<NSteps; i++)
            Bz[i] = cs + Gx[i]*r[0] + Gy[i]*r[1] + Gz[i]*r[2];
        for (i=0; i<NSteps; i++)                                       // Loop over time
        {
            // calculate cosine and sine of rotation angles
            norm = sqrt(Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]* Bz[i])+ 0.000000000000001;
            angle = norm*dt;
            ct = cos(angle);
            st = sin(angle);
            nx = Bx[i]/norm;
            ny = By[i]/norm;
            nz = Bz[i]/norm;
            dotproduct = nx*vx + ny*vy + nz*vz;
            // Use Rodriguez Formula for rotation:
            Mout[0] = vx*ct + (ny*vz - nz*vy)*st + dotproduct*(1-ct)*nx;
            Mout[1] = vy*ct + (nz*vx - nx*vz)*st + dotproduct*(1-ct)*ny;
            Mout[2] = vz*ct + (nx*vy - ny*vx)*st + dotproduct*(1-ct)*nz;
            vx = Mout[0];
            vy = Mout[1];
            vz = Mout[2];
        }
    }
    
    // Free memory
    free(Bx);
    free(By);
    free(Bz);
    
}