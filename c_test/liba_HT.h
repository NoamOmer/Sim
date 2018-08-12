/*
 * MATLAB Compiler: 4.14 (R2010b)
 * Date: Tue Jul 26 11:53:07 2011
 * Arguments: "-B" "macro_default" "-W" "lib:liba_HT" "-T" "link:lib" "HT" 
 */

#ifndef __liba_HT_h
#define __liba_HT_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_liba_HT
#define PUBLIC_liba_HT_C_API __global
#else
#define PUBLIC_liba_HT_C_API /* No import statement needed. */
#endif

#define LIB_liba_HT_C_API PUBLIC_liba_HT_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_liba_HT
#define PUBLIC_liba_HT_C_API __declspec(dllexport)
#else
#define PUBLIC_liba_HT_C_API __declspec(dllimport)
#endif

#define LIB_liba_HT_C_API PUBLIC_liba_HT_C_API


#else

#define LIB_liba_HT_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_liba_HT_C_API 
#define LIB_liba_HT_C_API /* No special import/export declaration */
#endif

extern LIB_liba_HT_C_API 
bool MW_CALL_CONV liba_HTInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_liba_HT_C_API 
bool MW_CALL_CONV liba_HTInitialize(void);

extern LIB_liba_HT_C_API 
void MW_CALL_CONV liba_HTTerminate(void);



extern LIB_liba_HT_C_API 
void MW_CALL_CONV liba_HTPrintStackTrace(void);

extern LIB_liba_HT_C_API 
bool MW_CALL_CONV mlxHT(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

extern LIB_liba_HT_C_API 
long MW_CALL_CONV liba_HTGetMcrID();



extern LIB_liba_HT_C_API bool MW_CALL_CONV mlfHT(int nargout, mxArray** y, mxArray* x);

#ifdef __cplusplus
}
#endif
#endif
