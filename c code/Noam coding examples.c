
mexPrintf("MEX-file: Printing input structure\n");
// ----------------------------------------------------
// Get 'double' from field number 3 inside input struct
// ----------------------------------------------------
mxArray *p1;
double f1;
p3 = mxGetFieldByNumber(prhs[0],0,2);
f3 = *mxGetPr(mxGetFieldByNumber(prhs[0],0,2));
mexPrintf("   p1=%p\n",p1);
mexPrintf("   f1=%f\n",f1);

// ----------------------------------------------------------
// Get 'double' array from field number 1 inside input struct
// ----------------------------------------------------------
mxArray *p1;
double *f1;
p1 = mxGetFieldByNumber(prhs[0],0,0);
f1 = mxGetPr(mxGetFieldByNumber(prhs[0],0,0));
mexPrintf("   p1=%p\n",p1);
mexPrintf("   f1[0]=%f\n",f1[0]);
mexPrintf("   f1[1]=%f\n",f1[1]);
mexPrintf("   f1[2]=%f\n",f1[2]);

// --------------------------------------------------------
// Get string array from field number 5 inside input struct
// --------------------------------------------------------
mxArray *p      = NULL;
int      buflen = 0,
	     fnum   = 4,
         status = 0;
char    *buf    = NULL;

// Input must be a string
p = mxGetFieldByNumber(prhs[0],0,4);
if (mxIsChar(p) != 1)
	mexErrMsgTxt("Input must be a string.\n");
// Input must be a row vector.
if (mxGetM(p) != 1)
	mexErrMsgTxt("Input must be a row vector.");
// Get the length of the input string.
buflen = (mxGetM(p) * mxGetN(p) + 1);
// Allocate memory for input string
buf = mxCalloc(buflen, sizeof(char));

// Copy the string data from prhs[0] into a C string input_buf
status = mxGetString(p, buf, buflen);
if (status != 0) 
	mexWarnMsgTxt("Not enough space. String is truncated.");

mexPrintf("Field number %d is: %s\n",fnum,buf);


