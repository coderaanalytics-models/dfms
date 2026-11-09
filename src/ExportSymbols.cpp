#include <Rcpp.h>

RcppExport SEXP _dfms_SKF(SEXP XSEXP, SEXP ASEXP, SEXP CSEXP, SEXP QSEXP, SEXP RSEXP, SEXP F_0SEXP, SEXP P_0SEXP, SEXP retLLSEXP);
RcppExport SEXP _dfms_FIS(SEXP ASEXP, SEXP ZTfSEXP, SEXP ZTpSEXP, SEXP VTf_vSEXP, SEXP VTp_vSEXP, SEXP F_0SEXP, SEXP P_0SEXP);
RcppExport SEXP _dfms_SKFS(
    SEXP XSEXP,
    SEXP ASEXP,
    SEXP CSEXP,
    SEXP QSEXP,
    SEXP RSEXP,
    SEXP F_0SEXP,
    SEXP P_0SEXP,
    SEXP retLLSEXP,
    SEXP contributionsSEXP
);
RcppExport SEXP _dfms_Estep(SEXP XSEXP, SEXP ASEXP, SEXP CSEXP, SEXP QSEXP, SEXP RSEXP, SEXP F_0SEXP, SEXP P_0SEXP);
RcppExport SEXP _dfms_ainv(SEXP FsEXP);
RcppExport SEXP _dfms_apinv(SEXP FsEXP);
RcppExport SEXP _dfms_update_loadings_cpp(
    SEXP XSEXP,
    SEXP F_smoothSEXP,
    SEXP P_smoothSEXP,
    SEXP loading_mapsSEXP
);
RcppExport SEXP _dfms_update_factor_blocks_cpp(
	SEXP F_smoothSEXP,
    SEXP P_smoothSEXP,
    SEXP PPm_smoothSEXP,
    SEXP F_smooth_0SEXP,
    SEXP P_smooth_0SEXP,
    SEXP factor_layoutSEXP
);


static const R_CallMethodDef CallEntries[] = {
  {"Cpp_SKF", (DL_FUNC) &_dfms_SKF, 8},
  {"Cpp_FIS", (DL_FUNC) &_dfms_FIS, 7},
  {"Cpp_SKFS", (DL_FUNC) &_dfms_SKFS, 9},
  {"Cpp_Estep", (DL_FUNC) &_dfms_Estep, 7},
  {"Cpp_ainv", (DL_FUNC) &_dfms_ainv, 1},
  {"Cpp_apinv", (DL_FUNC) &_dfms_apinv, 1},
  {"_dfms_update_loadings_cpp", (DL_FUNC) &_dfms_update_loadings_cpp, 4},
  {"_dfms_update_factor_blocks_cpp", (DL_FUNC) &_dfms_update_factor_blocks_cpp, 6},
  {NULL, NULL, 0}
};

RcppExport void R_init_dfms(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
