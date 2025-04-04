#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP cpbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP clbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* extern SEXP cmbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); */
extern SEXP cpwbart(SEXP, SEXP, SEXP);
extern SEXP cwbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cmebart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cgbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/*extern SEXP cgbmm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);*/
extern SEXP cabart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/*extern SEXP cspbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);*/
extern SEXP mc_cores_openmp(void);
extern SEXP crtnorm(SEXP, SEXP, SEXP, SEXP);
extern SEXP crtgamma(SEXP, SEXP, SEXP, SEXP);
extern SEXP cdraw_lambda_i(SEXP, SEXP, SEXP, SEXP);
extern SEXP cdnorm(SEXP, SEXP, SEXP);



static const R_CallMethodDef CallEntries[] = {
    {"cpbart", (DL_FUNC)&cpbart, 27},
    {"clbart", (DL_FUNC)&clbart, 24},
    /*  {"cmbart",  (DL_FUNC) &cmbart,  29},*/
    {"cpwbart", (DL_FUNC)&cpwbart, 3},
    {"cwbart", (DL_FUNC)&cwbart, 31},
    {"cgbart", (DL_FUNC)&cgbart, 30},
    /*  {"cgbmm",   (DL_FUNC) &cgbmm,   34}, */
    {"cabart", (DL_FUNC)&cabart, 31},
    {"cmebart", (DL_FUNC)&cmebart, 32},
    {"mc_cores_openmp", (DL_FUNC)&mc_cores_openmp, 0},
    {"crtnorm", (DL_FUNC)&crtnorm, 4},
    {"crtgamma", (DL_FUNC)&crtgamma, 4},
    {"cdraw_lambda_i", (DL_FUNC)&cdraw_lambda_i, 4},
    {"cdnorm", (DL_FUNC)&cdnorm, 3},
    {NULL, NULL, 0}};

void R_init_meBART(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
