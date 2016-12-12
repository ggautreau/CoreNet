#include "nem_typ.h"    /* DataT, ... */
#include "nem_ver.h"    /* PrintVersions */
#include "nem_hlp.h"    /* PrintHelpGeneral ... */
#include "lib_io.h"     /* ReadOpeningComments */ /*V1.06-a*/
#include "genmemo.h"    /* GenAlloc, ... */ /*V1.06-a*/
#include <stdio.h>      /* printf, ... */
#include <stdlib.h>     /* atof, ... */
#include <string.h>     /* strncpy, ... */
#include <time.h>       /* time() */          /*V1.04-e*/

#include "nem_arg.h"    /* prototypes */

#define DEFAULT_ALGO         ALGO_NEM
#define DEFAULT_BETA         1.0
#define DEFAULT_BTAMODE      BETA_FIX           /*V1.04-b*/
#define DEFAULT_BTAHEUSTEP   0.1                /*V1.04-b*/
#define DEFAULT_BTAHEUMAX    2.0                /*V1.04-b*/
#define DEFAULT_BTAHEUDDROP  0.8             	/*V1.04-b*/
#define DEFAULT_BTAHEUDLOSS  0.5             	/*V1.04-b*/
#define DEFAULT_BTAHEULLOSS  0.02            	/*V1.04-b*/
#define DEFAULT_BTAGRADNIT   1               	/*V1.06-f*/
#define DEFAULT_BTAGRADCVTH  0.001           	/*V1.06-f*/
#define DEFAULT_BTAGRADSTEP  0.0            	/*V1.06-f*/
#define DEFAULT_BTAGRADRAND  0            	/*V1.06-f*/
#define DEFAULT_CRIT            CRIT_M             /*V1.04-h*/
#define DEFAULT_CVTHRES     0.04               /*V1.04-d*/
#define DEFAULT_CVTEST       CVTEST_CLAS        /*V1.06-g*/
#define DEFAULT_FAMILY       FAMILY_BERNOULLI      /*V1.06-b*/
#define DEFAULT_DISPER        DISPER_KD          /*V1.06-b*/
#define DEFAULT_PROPOR      PROPOR_K         /*V1.06-b*/
#define DEFAULT_FORMAT      FORMAT_HARD
#define DEFAULT_INIT              INIT_RANDOM       
#define DEFAULT_MISSING      MISSING_REPLACE    /*V1.05-a*/
#define DEFAULT_SORTEDVAR    0
#define DEFAULT_NEIGHSPEC     NEIGH_FILE
#define DEFAULT_NBITERS         500
#define DEFAULT_NBEITERS       1
#define DEFAULT_NBRANDINITS  30                 /*V1.06-h*/
#define DEFAULT_ORDER        ORDER_RANDOM       /*V1.04-f*/
#define DEFAULT_UPDATE       UPDATE_SEQ         /*V1.06-d*/
#define DEFAULT_TIE          TIE_RANDOM
         /*V1.06-e*/


typedef enum
{

  HELP_GENERAL,
  HELP_OPTIONS,
  HELP_EXAMPLES,
  HELP_FILEIN,
  HELP_FILEOUT,
  HELP_VERSIONS,
  HELP_NB

} HelpET ;



const char *CritStrVC[ CRIT_NB ] = { "U" , "M", "D" , "L" } ;



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
        )
/* ------------------------------------------------------------------- */
{
  const char*  func = "NemArgs" ;
  int nk = *Nclasses;


  StatModelP->Spec.K = nk ;

  /*	ImageP->Nl    = *Ncol ; /*c'est bien les lignes sont les colonnes*/
  /*    ImageP->Nc    = *Nrow;*/
        DataP->NbPts  = *NInds ;
        DataP->NbVars = *Nvars;
        *SpatialTypeP = TYPE_SPATIAL ;

  StatModelP->Para.Prop_K    = GenAlloc( nk, sizeof(float),
					 1, func, "Prop_K" ) ;
  StatModelP->Para.Disp_KD   = GenAlloc( nk * DataP->NbVars, sizeof(float),
				       1, func, "Disp_KD" ) ;
  StatModelP->Para.Center_KD = GenAlloc( nk * DataP->NbVars, sizeof(float),
					 1, func, "Center_KD" ) ;

  StatModelP->Para.NbObs_K   = GenAlloc( nk, sizeof(float),
					 1, func, "NbObs_K" ) ;
  StatModelP->Para.NbObs_KD  = GenAlloc( nk * DataP->NbVars, sizeof(float),
					 1, func, "NbObs_KD" ) ;
  StatModelP->Para.Iner_KD   = GenAlloc( nk * DataP->NbVars, sizeof(float),
					 1, func, "NbObs_KD" ) ;

  StatModelP->Desc.DispSam_D = GenAlloc( DataP->NbVars, sizeof(float),
					  1, func, "DispSam_D" );
  StatModelP->Desc.MiniSam_D = GenAlloc( DataP->NbVars, sizeof(float),
					  1, func, "MiniSam_D" );
  StatModelP->Desc.MaxiSam_D = GenAlloc( DataP->NbVars, sizeof(float),
					  1, func, "MaxiSam_D" );

  /* Set default value of optional parameters */
  //  StatModelP->Spec.ClassFamily = DEFAULT_FAMILY ;
  StatModelP->Spec.ClassDisper = DEFAULT_DISPER ;
  StatModelP->Spec.ClassPropor = DEFAULT_PROPOR ;


  NemParaP->Algo          = DEFAULT_ALGO ;
  StatModelP->Para.Beta   = DEFAULT_BETA ;          /*V1.06-b*/
  StatModelP->Spec.BetaModel = DEFAULT_BTAMODE ;       /*V1.04-b*/
  NemParaP->BtaHeuStep    = DEFAULT_BTAHEUSTEP ;    /*V1.04-b*/
  NemParaP->BtaHeuMax     = DEFAULT_BTAHEUMAX ;
  NemParaP->BtaHeuDDrop   = DEFAULT_BTAHEUDDROP ;
  NemParaP->BtaHeuDLoss   = DEFAULT_BTAHEUDLOSS ;
  NemParaP->BtaHeuLLoss   = DEFAULT_BTAHEULLOSS ;
  NemParaP->BtaPsGrad.NbIter    = DEFAULT_BTAGRADNIT  ;/*V1.06-g*/
  NemParaP->BtaPsGrad.ConvThres = DEFAULT_BTAGRADCVTH ;
  NemParaP->BtaPsGrad.Step      = DEFAULT_BTAGRADSTEP ;
  NemParaP->BtaPsGrad.RandInit  = DEFAULT_BTAGRADRAND ;
  NemParaP->Crit          = DEFAULT_CRIT ;          /*V1.04-h*/
  NemParaP->CvThres       = DEFAULT_CVTHRES ;       /*V1.04-d*/
  NemParaP->CvTest        = CVTEST_CLAS ;           /*V1.06-g*/
  NemParaP->DoLog         = FALSE ;                 /*V1.03-a previously TRUE*/
  NemParaP->NbIters       = DEFAULT_NBITERS ;
  NemParaP->NbEIters      = DEFAULT_NBEITERS ;
  /*NemParaP->NbRandomInits = DEFAULT_NBRANDINITS*DataP->NbVars ;  /*V1.06-h*/
  NemParaP->Seed          = time( NULL ) ;          /*V1.04-e*/
  NemParaP->Format        = DEFAULT_FORMAT ;
  NemParaP->InitMode      = DEFAULT_INIT ;
  NemParaP->SortedVar     = DEFAULT_SORTEDVAR ;
  NemParaP->NeighSpec     = DEFAULT_NEIGHSPEC ;
  NemParaP->VisitOrder    = DEFAULT_ORDER ;         /*V1.04-f*/
  NemParaP->SiteUpdate    = DEFAULT_UPDATE ;        /*V1.06-d*/
  NemParaP->TieRule       = DEFAULT_TIE ;           /*V1.06-e*/
  NemParaP->Debug         = FALSE ;                 /*V1.04-g*/

} 

/* end of NemArgs() */

