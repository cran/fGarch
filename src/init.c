#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* in ./llhGarch.f :
 *       SUBROUTINE GARCHLLH(N, Y, Z, H, NF, X, DPARM, MDIST, MYPAR, F)
 */
extern void F77_NAME(garchllh)(int *N,
			       double *Y, double *Z, double *H,
			       int *NF, double *X,
			       double *DPARM,
			       int *MDIST, int *MYPAR,
			       double *F);

static const R_FortranMethodDef FortranEntries[] = {
    {"garchllh", (DL_FUNC) &F77_NAME(garchllh), 10},
    {NULL, NULL, 0}
};

void R_init_fGarch(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
