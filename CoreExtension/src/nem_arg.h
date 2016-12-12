#include "nem_typ.h"    /* DataT, ... */


void NemArgs
        (
	  int *NInds,
	  int *Nvars,
	  int *Nclasses,
	  StatModelT    *StatModelP,    /* O */
          NemParaT      *NemParaP,      /* O */
	  DataT*        DataP,          /* O */
	  PtNeighsT*  PtsNeighsVP,         /* O */
	  TypeET*       SpatialTypeP    /* O */
	 );


extern const char *CritStrVC[ CRIT_NB ] ;
