#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" SEXP catmech_marginal_gibbs_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP catmech_marginal_hw_cpp(SEXP, SEXP, SEXP);
extern "C" SEXP catmech_marginal_predict_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP catmech_marginal_derived_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP catmech_marginal_ratio_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
  {"catmech_marginal_gibbs_cpp", (DL_FUNC) &catmech_marginal_gibbs_cpp, 5},
  {"catmech_marginal_hw_cpp", (DL_FUNC) &catmech_marginal_hw_cpp, 3},
  {"catmech_marginal_predict_cpp", (DL_FUNC) &catmech_marginal_predict_cpp, 6},
  {"catmech_marginal_derived_cpp", (DL_FUNC) &catmech_marginal_derived_cpp, 5},
  {"catmech_marginal_ratio_cpp", (DL_FUNC) &catmech_marginal_ratio_cpp, 9},
  {NULL, NULL, 0}
};

extern "C" void R_init_catmech(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
