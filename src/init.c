/*
 * meBART: Bayesian Additive Regression Trees with Measurement Error
 * Copyright (C) 2025 Kevin McCoy, Zachary Wooten, and Christine Peterson
 *
 * This package is a modification of the BART package originally created by 
 * Robert McCulloch and Rodney Sparapani, who own the original copyright.
 * Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/GPL-2
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP cpmebart(SEXP, SEXP, SEXP);
extern SEXP cmebart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mc_cores_openmp(void);

static const R_CallMethodDef CallEntries[] = {
    {"cpmebart", (DL_FUNC)&cpmebart, 3},
    {"cmebart", (DL_FUNC)&cmebart, 35},
    {"cpbart", (DL_FUNC)&cpbart, 31},
    {"mc_cores_openmp", (DL_FUNC)&mc_cores_openmp, 0},
    {NULL, NULL, 0}};

void R_init_meBART(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
