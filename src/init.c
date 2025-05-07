#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP cpmebart(SEXP, SEXP, SEXP);
extern SEXP cmebart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mc_cores_openmp(void);

static const R_CallMethodDef CallEntries[] = {
    {"cpmebart", (DL_FUNC)&cpmebart, 3},
    {"cmebart", (DL_FUNC)&cmebart, 35},
    {"mc_cores_openmp", (DL_FUNC)&mc_cores_openmp, 0},
    {NULL, NULL, 0}};

void R_init_meBART(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
