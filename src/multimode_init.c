#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP BinDist2(SEXP, SEXP, SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(calcdist)(void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"BinDist2", (DL_FUNC) &BinDist2, 5},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"calcdist", (DL_FUNC) &F77_NAME(calcdist), 7},
    {NULL, NULL, 0}
};

void R_init_multimode(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
