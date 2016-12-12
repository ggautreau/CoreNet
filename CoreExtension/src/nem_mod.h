/*\
    Prototypes of nem_noi.c exported functions

    1.05-a    12-JAN-1997  MissMode in ParaP*V*I
    1.05-b    17-JAN-1997  EmptyK_P and StatusET return in ParaP*V*I
    1.06-a    28-JUN-1998  GetDensityFunc <- nem_alg.c and del DensPkVkI
    1.06-b    28-JUN-1998  New EstimPara and del ParaP*V*I
\*/

#include "nem_typ.h"    /* NoiseParaT, ... */

#if 0
/*V1.06-a*/
int GetDensityFunc  /* STS_OK or STS_E_FUNCARG */
        (
            const ModelSpecT  *SpecP,           /* I */
            CompuDensFT**     CompuDensFP       /* O */
        ) ;

#endif

void ComputeFki
 (
  const DataT*      DataP,      /* I: observed records */
  const ModelSpecT* SpecP,      /* I: model specification */
  const ModelParaT* ParaP,      /* I: model parameters */
  double*           PkFki_NK,   /* O: f(x_i|para_k), i=0..N-1, k=0..K-1 */
  float*            LogPkFki_NK /* O: log f(x_i|para_k) */
 ) ;


StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimPara 
(
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* O : estimated parameters */
 ) ;


