/*\

    NEM_EXE.C

    Neighborhood EM project : 
    main program and input/output functions.

    June 96

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    This program clusters a given data set, using the 
    Neighborhood EM algorithm proposed by Ambroise 1996.
    All parameters are to be specified on the command line.
    Command syntax can be obtained by typing : 'nem_exe'
    with no argument. Longer help is obtained by typing :
    'nem_exe -h'.

    The program takes as input  a set of objects described, 
    on one hand, by quantitative variables (the 'observed data'), 
    and on the other hand (optionnally), by their neighborhood 
    relationships in the geographic space. It outputs  
    a partition of the objects ( fuzzy or hard ) and class 
    parameters. The outputs are 'locally optimal' according to 
    a 'penalized likelihood' criterion (see below). Only 
    two parameters are required: the number of classes, 
    and the 'beta' coefficient that measures the desired 
    spatial smoothing. The processus is otherwise entirely 
    unsupervised. Other optional parameters may be specified 
    to fit better to the problem at hand.

    The 'penalized likelihood' criterion is designed in order to
    produce a partition and class parameters that fit well to 
    the data (likelihood term), and a partition that is spatially
    homogeneous (penalization term). It is possible to give
    more or less importance to the spatial information by
    increasing or decreasing the 'beta' coefficient.

    The algorithm finds a 'local' maximum of the criterion
    through an iterative processus, akin to gradient ascent : 
    starting from an initial arbitrary solution, it computes 
    successive solutions that gradually improve the criterion ; 
    thence, the final partition depends on the initial partition.
    This program can either read a given initial partition,
    compute it by thresholding one variable's histogram,
    or generate random initial partitions.
\*/

/*\
    NEM_EXE.C

Vers-mod  Date         Who  Description

1.03-a    03-NOV-1996  MD   Add prototype of ReadLabelFile() in main()
1.03-b    03-NOV-1996  MD   Save Pk, Vk, and Ck in don.mf
1.03-c    22-AUG-1997  MD   Unknown labels : memberships set to 0
1.03-d    22-AUG-1997  MD   Save criteria to don.mf file
1.03-e    30-SEP-1997  MD   Save criterion L to don.mf file
1.03-f    30-SEP-1997  MD   Save algorithm parameters into don.mf file
1.04-a    04-OCT-1997  MD   Any name for initial partition file
1.04-b    09-OCT-1997  MD   If "-o -" classification saved to standard output
1.04-c    09-OCT-1997  MD   If "-s f -" initial classification from stdin
1.04-d    10-OCT-1997  MD   Random seed set by option
1.04-e    10-OCT-1997  MD   Visit order ("-O random")
1.04-f    13-OCT-1997  MD   Read reference classification
1.04-g    05-NOV-1997  MD   Call to srand48() replaced by srandom()
1.04-h    02-DEC-1997  MD   Add "GEM" to AlgoDesC
1.05-a    12-JAN-1998  MD   Count missing data
1.05-b    17-JAN-1998  MD   Program exit values changed
1.05-c    25-JAN-1998  MD   main() now mainfunc()
1.05-d    26-JAN-1998  MD   DataP->LabelV init to NULL by default
1.05-e    26-JAN-1998  MD   Free allocated data after processing
1.05-f    26-JAN-1998  MD   GenAlloc/GenFree instead of malloc/calloc/free
1.05-g    30-JAN-1998  MD   mainfunc() prototyped in mainfunc.h
1.05-h    30-JAN-1998  MD   Fixed bug : str file was opened and not closed
1.05-i    05-FEB-1998  MD   my_strupr returns void
1.05-j    05-FEB-1998  MD   save M (Markov pseudo-likelihood) into don.mf
1.06-a    17-JUN-1998  MD   use new StatModel structure (old NoisModel)
1.06-b    23-JUN-1998  MD   LEN_LINE -> nem_typ.h
1.06-c    20-SEP-1998  MD   change format of ref / label file (!!!)
1.06-d    20-SEP-1998  MD   struct ErrInfo added to CriterT
1.06-e    21-SEP-1998  MD   TieRule copied in Errinfo
1.07-a    26-FEB-1999  MD   Add Bernoulli family
1.07-b    26-FEB-1999  MD   Add "\n" at end of final classification file
\*/

#include "mainfunc.h"   /* Prototype of exported mainfunc() */

#include "nem_typ.h"    /* DataT, ... */
#include "nem_ver.h"    /* NemVersionStrC */
#include "nem_arg.h"    /* NemArg, ... */
#include "nem_alg.h"    /* ClassifyByNem, ... */
#include "nem_rnd.h"    /* RandomInteger, ... */
#include "lib_io.h"     /* ReadOpeningComments, ... */
#include "genmemo.h"    /* GenAlloc, ... */

#include <stdio.h>      /* printf, ... */
#ifdef __TURBOC__
#include <alloc.h>      /* coreleft, ... */ 
#endif
#include <string.h>     /* strncpy, ... */
#include <math.h>       /* sqrt, ... */



static const char *TypeDesC[ TYPE_NB ] = { "Spatial", "Image", "NoSpatial" } ;
static const char *AlgoDesC[ ALGO_NB ] = { "NEM" , "NCEM (C-step)" , 
					   "GEM (Monte-Carlo at E-step)" } ;
static const char *BetaDesVC[ BETA_NB ] = { "fixed", 
					    "pseudo-likelihood gradient",
					    "heuristic Hathaway crit",
					    "heuristic mixture likelihood"} ;
static const char *FamilyDesVC[ FAMILY_NB ] = { "Normal", "Laplace", 
						"Bernoulli" } ;
static const char *DisperDesVC[ DISPER_NB ] = { "S__", "SK_", "S_D", "S_KD" } ;
static const char *ProporDesVC[ PROPOR_NB ] = { "P_", "Pk" } ;
static const char *InitDesC[ INIT_NB ] = { "Sort_variables" ,
					   "Random_centers" , 
					   "Mixture_initial" ,
					   "Mixture_fixed" ,
					   "Full_labels" , 
					   "Partial_labels" } ;    /*V1.03-f*/

static const char *OrderDesC[ ORDER_NB ] = { "Direct_order" , "Random_order" };




/* ==================== LOCAL FUNCTION PROTOTYPING =================== */


/* Called by ClassifyByNem */

    static int GetInputPara 
        ( 
          int           Argc,
          const char    *Argv[],
          DataT         *DataP,         /* O and allocated */
          NemParaT      *NemParaP,      /* O */
          SpatialT      *SpatialP,      /* O and allocated */
          StatModelT    *StatModelP,    /* O and allocated */
	  ErrinfoT      *ErrinfoP,      /* O and allocated */
	  ErrcurT       *ErrcurP,       /* O and allocated */
          float         **ClassifMP     /* O and allocated */
        ) ;

    static int SaveResults
        ( 
	  const int          Argc,
	  const char*        Argv[],
          const int          Npt,                   /* I */
          const int          Nd,                    /* I */
          const float*       ClassifM,              /* I */
          const SpatialT*    SpatialP,              /* I */
          const NemParaT*    NemParaP,              /* I */
          const StatModelT*  StatModelP,            /* I */
	  const CriterT*     CriterP                /* I */ /*V1.03-d*/
        ) ;

static void FreeAllocatedData
(
 DataT*       DataP,       /* O and deallocated */
 SpatialT*    SpatialP,    /* O and deallocated */
 ModelParaT*  ModelParaP,  /* O and deallocated */  /*V1.06-a*/
 CriterT*     CriterP,     /* O and deallocated */  /*V1.06-d*/
 float*       ClassifM     /* O and deallocated */
) ;




/* Called by GetInputPara */


    static int SetVisitOrder  /*V1.04-e*/
        ( 
	     int         Npt,          /* I */
	     OrderET     VisitOrder,   /* I */
	     int**       SiteVisitVP   /* O and allocated (Npt) */
	) ;
                                     
    static int  ReadMatrixFile                     
         (                             
             const char  *FileName,     /* I */      /*V1.04-a*/
             int         Nl,            /* I */
             int         Nc,            /* I */
             float       **MatP         /* O and allocated (size : Nl * Nc) */
         ) ;                             
                                     
    static int  ReadLabelFile  /*V1.03-a*/
         ( 
             const char  *LabelName,    /* I */
             int         Npt,           /* I */
	     int         *KfileP,       /* O : file # classes */ /*V1.06-c*/
             int         **LabelVP,     /* O and allocated (Npt) */
             float       **ClassifMP    /* O and allocated (Npt*Nk) */
         );

    static int  MakeErrinfo
         ( 
             const char* RefName,       /* I : filename of reference class */
             int         N,             /* I : number of objects */
             int         Kc,            /* I : user number of classes */
	     TieET       TieRule,       /* I : specified MAP tie rule */
	     ErrinfoT*   ErrinfoP,      /* O and allocated */
	     ErrcurT*    ErrcurP        /* O and allocated */
         ) ;

    static int  ReadNeiFile                     
         (                              
             const char  *BaseName,     /* I */
             int         NbPts,         /* I */
             NeighET     NeighSpec,     /* I */
             char*       NeiDescS,      /* O [LEN_LINE+1] */
             SpatialT    *SpatialP      /* O and allocated */
         ) ;



/* Called by ReadNeiFile */

    static int  ReadPtsNeighs
                (
                    FILE        *Fnei,          /* I/O */
                    int         NbPts,          /* I */
                    int         *MaxNeiP,       /* O */
                    NeighDataT  *NeighDataP     /* O and allocated */
                ) ;

    static int  ReadImageNeigh
                (
                    FILE        *Fnei,          /* I/O */
                    NeighDataT  *NeighDataP     /* O and allocated */
                ) ;

    static int  SetImageNeigh
                (   
                    NeighET     NeighSpec,          /* I */
                    char*       NeiDescS,           /* O [LEN_LINE+1] */
                    NeighDataT* NeighDataP          /* O and allocated */
                ) ;


/* Called by MakeErrinfo */

static int factorial(int n);

static int compute_permutations    /* ret 0 if OK, -1 if memory error */
(
 const int Startval,         /* I : start value of integer suite */
 const int K,                /* I : size of integer suite to permute > 0 */
 int*      perms_Kfact_K_p[] /* O : matrix to store permuted values */
) ;


/* Called by compute_permutations */

static int rec_permutations        /* ret 0 if OK, -1 if memory error */
(
 const int array_A[],       /* I : remaining array to permute */
 const int A,               /* I : length of the remaining array : 0..K */
 const int K,               /* I : length of original array */
 int       offset,          /* I : first line of storage to use */
 int       perms_Kfact_K[]  /* O : matrix to store permuted values, use
			       lines :   offset -> offset + A! - 1
			       columns : K - A  -> K - 1 */
) ;


/* Called by */




/* ==================== GLOBAL FUNCTION DEFINITION =================== */


/* ------------------------------------------------------------------- */
int mainfunc( int argc, const char * * argv )
/*\
    NEM main program.
\*/
/* ------------------------------------------------------------------- */
{
    StatusET                err ;
    static  DataT           Data = {0} ;
    static  NemParaT        NemPara = {0} ;
    static  SpatialT        Spatial = {{{0}}} ;
    static  StatModelT      StatModel = {{0}} ;
    static  float           *ClassifM = 0;
    CriterT                 Criteria = {0} ; /*V1.03-d*/


    /* main program algorithm :
       - read all necessary data into memory
       - compute partition by NEM method
       - save results to file
    */

    fprintf( stderr , " * * * NEM (spatial data clustering) v%s * * *\n" ,
             NemVersionStrC ) ;

#ifdef __TURBOC__
    fprintf( stderr, "\n Initial free memory : %lu bytes\n\n", 
             (unsigned long) coreleft() );
#endif

    /* RandomSeedByTime( ) ;  V1.04-d*/

    if ( ( err = GetInputPara ( argc, argv, &Data, &NemPara, &Spatial, 
                                &StatModel, &Criteria.Errinfo, 
				&Criteria.Errcur, &ClassifM ) )
	 == STS_OK )
    {
#ifdef __TURBOC__
      srand( (unsigned) NemPara.Seed ) ;
#else
      srandom( NemPara.Seed ) ;  /*V1.04-g*/
#endif

        if ( ( err = ClassifyByNem( &NemPara, &Spatial, &Data, 
                                    &StatModel, ClassifM, 
				    &Criteria ) ) == STS_OK )
        {
            fprintf( stderr, "Saving results ...\n" ) ;
            SaveResults( argc, argv, Data.NbPts , Data.NbVars , ClassifM, 
                         &Spatial, &NemPara, &StatModel, &Criteria ) ;
        }

	FreeAllocatedData( &Data, &Spatial, &StatModel.Para, 
			   &Criteria, ClassifM ) ;
    }

    switch( err )
    {
        case STS_OK :
	  if ( strcmp( NemPara.OutBaseName , "-" ) != 0 )  /*V1.04-b*/
	    {
	      fprintf( stderr, "NEM completed, classification in %s\n",
		       NemPara.OutName ) ;
	      fprintf( stderr, " criteria and parameters in %s%s\n",
		       NemPara.OutBaseName, EXT_MFNAME ) ;

	      if ( NemPara.DoLog )
		{
		  fprintf( stderr, "Log of detailed running in %s\n",
			   NemPara.LogName );
		}
	    }
	  else
	    {
	      fprintf( stderr, "NEM completed, classification to screen\n") ;
	    }
	  return EXIT_OK ;

        case STS_W_EMPTYCLASS :
             fprintf( stderr, "*** NEM warning status : empty class\n" );
             return EXIT_W_RESULT ;
             
        case STS_I_DONE :
             return EXIT_OK ;

        case STS_E_ARG :
             PrintUsage( argv[ 0 ] ) ;
             return EXIT_E_ARGS ;

        case STS_E_MEMORY :
             fprintf( stderr, "*** NEM error status : not enough memory\n" );
#ifdef __TURBOC__
             fprintf( stderr, "\n Memory left : %lu bytes\n", 
                      (unsigned long) coreleft() );
#endif
             return EXIT_E_MEMORY ;

        case STS_E_FILEIN :
             fprintf( stderr, "*** NEM error status : can't open a file\n" );
             return EXIT_E_FILE ;

        case STS_E_FILE :
             fprintf( stderr, "*** NEM error status : wrong file format\n" );
             return EXIT_E_FILE ;

        case STS_E_FUNCARG :
             fprintf( stderr, "*** NEM internal error : bad arguments\n" );
             return EXIT_E_BUG ;

        default :
             fprintf( stderr, "*** NEM unknown error status (%d)\n", err );
             return EXIT_E_BUG ;
    } 
} /* end of mainfunc() */



/* ==================== LOCAL FUNCTION DEFINITION =================== */


/* ------------------------------------------------------------------- */
static int GetInputPara 
    (
          int           Argc,
          const char    *Argv[],
          DataT         *DataP,         /* O and allocated */
          NemParaT      *NemParaP,      /* O */
          SpatialT      *SpatialP,      /* O and allocated */
          StatModelT    *StatModelP,    /* O and allocated */
	  ErrinfoT      *ErrinfoP,      /* O and allocated */
	  ErrcurT       *ErrcurP,       /* O and allocated */
          float         **ClassifMP     /* O and allocated */
    ) 
/* ------------------------------------------------------------------- */
{
    StatusET    err ;
    const char* func = "GetInputPara" ;
    char        fname[ LEN_FILENAME + 1 ] ;
    char        namedat[ LEN_FILENAME + 1 ] ;   /*V1.04-a*/
    char        datadescS[ LEN_LINE + 1 ] ;
    char        neidescS[ LEN_LINE + 1 ] ;
    int         klabelfile ;


    /* Start of GetInputPara() */

    if ( ( err = NemArgs( Argc, Argv, fname, StatModelP, NemParaP,
			  datadescS, DataP, &SpatialP->NeighData.Image,
			  &SpatialP->Type ) ) 
        != STS_OK )
      {
        return err ;
      }

    /* Read points */
    fprintf( stderr, "Reading points ...\n" ) ;
    strncpy( namedat , fname , LEN_FILENAME ) ;  /*V1.04-a*/
    strncat( namedat , ".dat", LEN_FILENAME ) ;


    if ( ( err = ReadMatrixFile( namedat,        /*V1.04-a*/
                                 DataP->NbPts, 
                                 DataP->NbVars, 
                                 &DataP->PointsM ) ) != STS_OK )
       return err ;

    /* Count missing data */ /*V1.05-a*/
    {
      int i, j ;

      DataP->NbMiss = 0 ;
      for ( i = 0 ; i < DataP->NbPts ; i ++ )
	{
	  for ( j = 0 ; j < DataP->NbVars ; j ++ )
	    {
	      if ( isnan( DataP->PointsM[ i * DataP->NbVars + j ] ) )
		DataP->NbMiss ++ ;
	    }
	}
    }

    /* Allocate and set sites visit order */
    if ( ( err = SetVisitOrder( DataP->NbPts,   /*V1.04-e*/
				NemParaP->VisitOrder,
				& DataP->SiteVisitV ) ) != STS_OK )
      return err ;


    /* Allocate and eventually read initial classification matrix */
    DataP->LabelV = NULL ; /*V1.05-d*/
    switch( NemParaP->InitMode )
    {
    case INIT_FILE:
        fprintf( stderr, "Reading initial partition ...\n" ) ;
        if ( ( err = ReadMatrixFile( NemParaP->StartName,     /*V1.04-a*/
				     DataP->NbPts,
                                     StatModelP->Spec.K, 
                                     ClassifMP ) ) != STS_OK )
            return err ;
        break ;

    case INIT_SORT:    /* allocate partition (will be initialized later) */
    case INIT_RANDOM:  /* allocate partition (will be initialized later) */
    case INIT_MIXINI:
    case INIT_MIXFIX:

        /* Allocate classification matrix */
        if ( ( (*ClassifMP) = 
	       GenAlloc( DataP->NbPts * StatModelP->Spec.K, sizeof( float ),
			 0, func, "(*ClassifMP)" ) ) == NULL )
	  return STS_E_MEMORY ;
        break ;

    case INIT_LABEL:
      fprintf( stderr, "Reading known labels file ...\n" ) ;
      if ( ( err = ReadLabelFile( NemParaP->LabelName, DataP->NbPts, 
				  & klabelfile,
                                  & DataP->LabelV,
                                  ClassifMP ) ) != STS_OK )
        return err ;

      if ( klabelfile != StatModelP->Spec.K ) {
	fprintf( stderr, 
		 "Error : label file %d classes, command line %d classes\n",
		 klabelfile, StatModelP->Spec.K ) ;
	return STS_E_FILE ;
      }
      break;

    default: /* error */
        fprintf( stderr, "Unknown initialization mode (%d)\n", 
                 NemParaP->InitMode );
        return STS_E_FUNCARG ;
    }

    /* Eventually read reference class file */  /*V1.04-f*/
    if ( ( err = MakeErrinfo( NemParaP->RefName, DataP->NbPts, 
			      StatModelP->Spec.K, NemParaP->TieRule, 
			      ErrinfoP, ErrcurP ) ) != STS_OK )
      return err ;

    /* Read neighborhood file */
    if ( SpatialP->Type != TYPE_NONSPATIAL )
    {
        fprintf( stderr, "Reading neighborhood information ...\n" ) ;
        if ( ( err = ReadNeiFile( fname, 
                                  DataP->NbPts, 
                                  NemParaP->NeighSpec ,
                                  neidescS,
                                  SpatialP ) ) != STS_OK )
                return err ;
    }
    else
    {
        StatModelP->Para.Beta = 0.0 ; /*V1.06-a*/
        SpatialP->MaxNeighs   = 0 ;
    }

    fprintf( stderr, "\nData : " ) ;
    if ( strcmp( datadescS, "" ) ) fprintf( stderr, "%s\n", datadescS ) ;
    else                          fprintf( stderr, "\n" ) ;

    fprintf( stderr, "  file names =  %10s   |   nb points   = %10d\n", 
             fname, DataP->NbPts ) ;
    fprintf( stderr, "  type       =  %10s   |   dim         = %10d\n",
             TypeDesC[ SpatialP->Type ], DataP->NbVars ) ;
    if ( SpatialP->Type == TYPE_IMAGE )
    {
    fprintf( stderr, "  image size =  (%4d,%4d)\n", 
             SpatialP->NeighData.Image.Nl, SpatialP->NeighData.Image.Nc ) ;
    }
    if ( DataP->NbMiss > 0 ) /*V1.05-a*/
    {
    fprintf( stderr, "  %d missing values / %d\n",
	     DataP->NbMiss, DataP->NbPts * DataP->NbVars ) ;
    }

    if ( SpatialP->Type != TYPE_NONSPATIAL )
    {
    fprintf( stderr, "Neighborhood system :\n  max neighb =  %10d\n",
             SpatialP->MaxNeighs ) ;
    fprintf( stderr, "%s\n", neidescS ) ;
    }

    fprintf( stderr, "\n" ) ;
    fprintf( stderr, "NEM parameters :\n" ) ;
    fprintf( stderr, "Type of algorithm : '%s'\n" , AlgoDesC[ NemParaP->Algo ] ) ;
    fprintf( stderr, 
	     "  beta       =  %10.2f   |   nk                    = %3d\n", 
             StatModelP->Para.Beta, StatModelP->Spec.K) ;
    fprintf( stderr, 
	     "                %10s   |   model                 = %s, %s %s\n",
             " ", 
	     FamilyDesVC[ StatModelP->Spec.ClassFamily ],
	     ProporDesVC[ StatModelP->Spec.ClassPropor ],
	     DisperDesVC[ StatModelP->Spec.ClassDisper ]) ;
    fprintf( stderr, "\n" ) ;


    return STS_OK ;

} /* end of GetInputPara() */


/* ------------------------------------------------------------------- */
static int SetVisitOrder   /*V1.04-e*/
        ( 
	     int         Npt,          /* I */
	     OrderET     VisitOrder,   /* I */
	     int**       SiteVisitVP   /* O and allocated (Npt) */
	) 
/* ------------------------------------------------------------------- */
{
  int ivis ;


  if (( *SiteVisitVP = GenAlloc( Npt, sizeof( int ), 0, 
				 "SetVisitOrder", "SiteVisitVP" ) ) == NULL )
      return STS_E_MEMORY ;

  for ( ivis = 0 ; ivis < Npt ; ivis ++ )
    (*SiteVisitVP)[ ivis ] = ivis ;

  if ( VisitOrder == ORDER_RANDOM )
    {
      RandomPermutationAlgo( (*SiteVisitVP) , Npt ) ;
    }

  return STS_OK ;
}


/* ------------------------------------------------------------------- */
static int ReadMatrixFile
         (
             const char  *FileName,    /* I */      /*V1.04-a*/
             int         Nl,           /* I */
             int         Nc,           /* I */
             float       **MatP        /* O and allocated (size : Nl * Nc) */
         ) 
/* ------------------------------------------------------------------- */
{
    float       *Mat ;
    FILE        *fp ;
    int         il ;
    int         ic ;
    float       x ;


    if ( ( *MatP = GenAlloc( Nl *  Nc, sizeof( float ),
			     0, "ReadMatrixFile", FileName ) ) == NULL )
      return STS_E_MEMORY ;
    Mat = *MatP ;

    if ( strcmp( FileName , "-" ) != 0 )   /*V1.04-c*/
      {
	if ( ( fp = fopen( FileName, "r" ) ) == NULL )
	  {
	    fprintf( stderr, "File %s does not exist\n", FileName ) ;
	    GenFree( Mat ) ;
	    return STS_E_FILEIN ;
	  }
      }
    else
      fp = stdin ;

    for ( il = 0, ic = 0 ; ( il < Nl ) && (! feof( fp )) ; il ++ )
    {
        /* printf( "\n pt %d : ", il ) ; */
        for ( ic = 0 ; ( ic < Nc ) && (! feof( fp )) ; ic ++ )
        {
            if ( fscanf( fp, "%f", &x) == 1 )
            {
                Mat[ ( il * Nc ) + ic ] = x ;
                /* printf( " %g", x ) ; */
            }
            else    ic -- ;
        }
    }

    if ( strcmp( FileName , "-" ) != 0 )   /*V1.04-c*/
      fclose( fp ) ;

    if ( ( il < Nl ) || ( ic < Nc ) )
    {
        if ( ic == 0 )
        {
            ic = Nc ;
            il -- ;
        }
        fprintf( stderr, "%s : short file (%d/%d lines and %d/%d columns)\n", 
                 FileName, il, Nl, ic, Nc ) ;
        GenFree( Mat ) ;
        return STS_E_FILE ;
    }
    else
        return STS_OK ;

}  /* end of ReadMatrixFile() */




/* ------------------------------------------------------------------- */
static int  ReadLabelFile
         ( 
             const char  *LabelName,    /* I */
             int         Npt,           /* I */
	     int         *KfileP,       /* O : file # classes */ /*V1.06-c*/
             int         **LabelVP,     /* O and allocated (Npt) */
             float       **ClassifMP    /* O and allocated (Npt*Nk) */
         )
/* ------------------------------------------------------------------- */
{
    FILE  *fp;
    int   pt ;

    if ( ( fp = fopen( LabelName, "r" ) ) == NULL )
    {
        fprintf( stderr, "File %s does not exist\n", LabelName ) ;
        return STS_E_FILEIN ;
    }

    /* Read number of classes */ /*V1.06-c*/
    fscanf( fp , "%d" , KfileP ) ;

    /* Allocate classification matrix and vector */
    if ( ( *ClassifMP = 
	   GenAlloc( Npt *  (*KfileP), sizeof( float ),
		     0, "ReadLabelFile", LabelName ) ) == NULL )
        return STS_E_MEMORY ;

    if ( ( *LabelVP = 
	   GenAlloc( Npt, sizeof( int ),
		     0, "ReadLabelFile", LabelName ) ) == NULL )
        return STS_E_MEMORY ;

    /* Read labels and store in matrix and vector */
    for ( pt = 0 ; ( pt < Npt ) && (! feof( fp )) ; pt ++ ) {
        /* Read label indication */
        fscanf( fp , "%d" , & (*LabelVP)[ pt ] ) ;

        /* If label is known and valid */
        if ( ( 0 < (*LabelVP)[ pt ] ) && ( (*LabelVP)[ pt ] <= (*KfileP) ) ) 
	  LabelToClassVector( (*KfileP), (*LabelVP)[ pt ] - 1, 
			      & (*ClassifMP)[ pt * (*KfileP) ] ) ;
	/* Else (label unknown or invalid) */
        else {
            /*V1.03-c*/
	    (*LabelVP)[ pt ] = 0 ; /* enforce label to 0 if it is invalid */
	    LabelToClassVector( (*KfileP), (*LabelVP)[ pt ] - 1, 
				& (*ClassifMP)[ pt * (*KfileP) ] ) ;
	}
      }

    fclose( fp ) ;


    if ( pt < Npt )
    {
        fprintf( stderr, "%s : short file (%d/%d labels)\n", 
                 LabelName, pt , Npt ) ;
        GenFree( (*ClassifMP) ) ; (*ClassifMP) = NULL ;
        GenFree( (*LabelVP) ) ;   (*LabelVP)   = NULL ;

        return STS_E_FILE ;
    }
    else
        return STS_OK ;

}  /* end of ReadLabelFile() */



/* ------------------------------------------------------------------- */
static int  MakeErrinfo
         ( 
             const char* RefName,       /* I : filename of reference class */
             int         N,             /* I : number of objects */
             int         Kc,            /* I : user number of classes */
	     TieET       TieRule,       /* I : specified MAP tie rule */
	     ErrinfoT*   ErrinfoP,      /* O and allocated */
	     ErrcurT*    ErrcurP        /* O and allocated */
         )
/* ------------------------------------------------------------------- */
{
  StatusET err ;
  int      *tmpV ;
  int      ipt ;

  if ( strcmp( RefName, "" ) != 0 ) {
    ErrinfoP->Kc = Kc ;
    if ( ( err = ReadLabelFile( RefName, N, 
				& ErrinfoP->Kr,
				& tmpV,
				& ErrinfoP->Refclas_N_Kr ) ) != STS_OK )
      return err ;

    /* Check all reference labels ok */
    for ( ipt = 0, err = STS_OK ; 
	  ( ipt < N ) && ( err == STS_OK ) ; ipt ++ ) {
      if ( ( tmpV[ ipt ] <= 0 ) || ( tmpV[ ipt ] > ErrinfoP->Kr ) ) {
	fprintf( stderr, 
		 "Reference class for point %d not in 1..%d \n", 
		 ipt + 1, ErrinfoP->Kr ) ;
	err = STS_E_FILE ;
      }
    }
    GenFree( tmpV ) ; tmpV = NULL ;
    if ( err != STS_OK )
      return err ;

    /* Compute greatest number of classes and permutations of classes */
    ErrinfoP->Km = ( ErrinfoP->Kc > ErrinfoP->Kr ) ? 
      ErrinfoP->Kc : ErrinfoP->Kr ;

    ErrinfoP->Kmfac = factorial( ErrinfoP->Km ) ;

    ErrinfoP->TieRule = TieRule ;

    compute_permutations( 0, ErrinfoP->Km, & ErrinfoP->Perm_Kmfac_Km ) ;

    /* Allocate and initialize later computed stuff */
    if ( ( ErrcurP->Agree_Km_Km = 
	   GenAlloc( ErrinfoP->Km * ErrinfoP->Km, sizeof( float ),
		     0, "MakeErrinfo", "Agree_Km_Km" ) ) == NULL )
        return STS_E_MEMORY ;

    if ( ( ErrcurP->Loclas_N_Kc = 
	   GenAlloc( N * ErrinfoP->Kc, sizeof( float ),
		     0, "MakeErrinfo", "Loclas_N_Kc" ) ) == NULL )
        return STS_E_MEMORY ;

    ErrcurP->Ibestpermut = -1 ;
    ErrcurP->Errorrate   = -2.0 ;

    return STS_OK ;
  }
  else {
    ErrinfoP->Kr = 0 ;
    ErrinfoP->Refclas_N_Kr = NULL ;
    ErrcurP->Errorrate     = -1.0 ;
    return STS_OK ;
  }

}  /* end of MakeErrinfo() */





/* ------------------------------------------------------------------- */
static int factorial(int n)
/* ------------------------------------------------------------------- */
{
  int result = 1;

  for (; n>0; n--)
    result *= n;

  return result;
}


/* ------------------------------------------------------------------- */
static int compute_permutations    /* ret 0 if OK, -1 if memory error */
(
 const int Startval,         /* I : start value of integer suite */
 const int K,                /* I : size of integer suite to permute > 0 */
 int*      perms_Kfact_K_p[] /* O : matrix to store permuted values */
)
/* ------------------------------------------------------------------- */
{
  int   ik ;      /* current index and value of integer suite : 0..K-1 */
  int*  array_K ; /* integer suite to permute : [ 0 ... (K-1) ] */
  int   Kfact ;   /* K! */
  int   err ;     /* 0 if OK, -1 if memory error */

  if ( K <= 0 )
    return 1 ;

  Kfact = factorial( K ) ;

  if( ( (*perms_Kfact_K_p) = malloc( Kfact * K * sizeof( int ) ) ) == NULL )
    return -1 ;

  if( ( array_K = malloc( K * sizeof( int ) ) ) == NULL ) {
    free( (*perms_Kfact_K_p) ) ;
    (*perms_Kfact_K_p) = NULL ;
    return -1 ;
  }

  for ( ik = 0 ; ik < K ; ik ++ )
    array_K[ ik ] = Startval + ik ;

  err = rec_permutations( array_K, K, K, 0, (*perms_Kfact_K_p) ) ;

  free( array_K ) ;

  return err ;
}


/* ------------------------------------------------------------------- */
static int rec_permutations        /* ret 0 if OK, -1 if memory error */
(
 const int array_A[],       /* I : remaining array to permute */
 const int A,               /* I : length of the remaining array : 0..K */
 const int K,               /* I : length of original array */
 int       offset,          /* I : first line of storage to use */
 int       perms_Kfact_K[]  /* O : matrix to store permuted values, use
			       lines :   offset -> offset + A! - 1
			       columns : K - A  -> K - 1 */
)
/* ------------------------------------------------------------------- */
{
  int    err ;         /* 0 if currently OK */
  int    ia ;          /* array element currently removed : 0..A-1 */
  int    ja ;          /* array element currently copied : 0..A-1 (neq ia) */
  int    am1fact ;     /* (A-1)! */
  int    iam1fact ;    /* current perms line (without offset) : 0..am1fact-1 */
  int*   redarr_am1 ;  /* array with one element removed (A-1) */

  am1fact = factorial( A - 1 ) ;

  /* Check against out of bounds offset */
  if ( ( offset < 0 ) || ( factorial( K ) < ( offset + A * am1fact ) ) )
    return 1 ;
  
  if ( ( redarr_am1 = malloc( ( A - 1 ) * sizeof( int ) ) ) == NULL ) 
    return -1;

  /* For each element of given array */
  for ( ia = 0, err = 0 ; ( ia < A ) && ( err == 0 ) ; ia ++ ) {

    /* Copy (A-1)! times this element into the column ( K - A ) of perms, 
       starting from line ( offset + ia * (A-1)! ) (skip previous ia's) */
    for ( iam1fact = 0 ; iam1fact < am1fact ; iam1fact ++ ) 
      perms_Kfact_K[ ( offset + ia * am1fact + iam1fact ) * K + ( K - A ) ] =
	array_A[ ia ] ;

    /* Copy array into reduced array without this element */
    for ( ja = 0 ; ja < ia ; ja ++ )
      redarr_am1[ ja ] = array_A[ ja ] ;
    for ( ja = ia + 1 ; ja < A ; ja ++ )
      redarr_am1[ ja - 1 ] = array_A[ ja ] ;

    /* Recursive call with reduced array */
    err = rec_permutations( redarr_am1, A-1, K, offset + ia * am1fact,
			    perms_Kfact_K ) ;
  }

  free( redarr_am1 ) ;
  return err ;
}




/* ------------------------------------------------------------------- */
static int  ReadNeiFile
         ( 
             const char  *BaseName,          /* I */
             int         NbPts,              /* I */
             NeighET     NeighSpec,          /* I */
	     char*       NeiDescS,           /* O [LEN_LINE+1] */
             SpatialT    *SpatialP           /* I/O and allocated */
         ) 
/* ------------------------------------------------------------------- */
{
    char        infname[ LEN_FILENAME + 1 ] ;
    FILE        *fnei ;
    StatusET    err = STS_OK ;

    /* Read neighbourhood data depending on spatial type */
    if ( SpatialP->Type == TYPE_SPATIAL )
    {
        strncpy( infname, BaseName, LEN_FILENAME ) ;
        strncat( infname, ".nei" , LEN_FILENAME ) ;

        if ( ReadOpeningComments( infname , "#" , LEN_LINE , 
                                  & fnei , NeiDescS ) == -1 )
        {
            fprintf( stderr, "File %s does not exist\n", infname ) ;
            return STS_E_FILEIN ;
        }

        err = ReadPtsNeighs( fnei, NbPts, 
                             &SpatialP->MaxNeighs, 
                             &SpatialP->NeighData ) ;
        fclose( fnei ) ;
    }
    else    /* Type is image */
    {
        if ( NeighSpec == NEIGH_FILE )
        {
            strncpy( infname, BaseName, LEN_FILENAME ) ;
            strncat( infname, ".nei" , LEN_FILENAME ) ;
            if ( ReadOpeningComments( infname , "#" , LEN_LINE , 
		                      & fnei , NeiDescS ) == -1 )
            {
                fprintf( stderr, "File %s does not exist\n", infname ) ;
                return STS_E_FILEIN ;
            }

            err = ReadImageNeigh( fnei, &SpatialP->NeighData ) ;
            fclose( fnei ) ;
        }
        else
        {
            err = SetImageNeigh( NeighSpec, NeiDescS, &SpatialP->NeighData ) ;
        }

        SpatialP->MaxNeighs = SpatialP->NeighData.Image.NbNeigh ;
    }


    return err ;

}  /* end of ReadNeiFile() */


/* ------------------------------------------------------------------- */
static int      ReadPtsNeighs
                (
                    FILE        *Fnei,          /* I/O */
                    int         NbPts,          /* I */
                    int         *MaxNeiP,       /* O */
                    NeighDataT  *NeighDataP     /* O and allocated */
                ) 
/* ------------------------------------------------------------------- */
{
    const char* func = "ReadPtsNeighs" ;
    int       weighted ;
    int       l ;
    int       ipt ;
    StatusET  err = STS_OK ;
    PtNeighsT *ptsneighsV ;
    int       nmax ;


    /* Indicator of weights presence */
    fscanf( Fnei, "%d", &weighted ) ;

    /* Allocate structure of all points' neighbours */
    if ( ( ptsneighsV = GenAlloc( NbPts, sizeof( PtNeighsT ), 
				  0, func, "ptsneighsV" ) )
           == NULL )
    {
        fprintf( stderr, "Cannot allocate list of neighbours of %d points\n" ,
                 NbPts ) ;
        err = STS_E_MEMORY ;
        return err ;
    }
    NeighDataP->PtsNeighsV = ptsneighsV ;

    /* For each point, default nb of neighbours is zero */
    for ( ipt = 0 ; ipt < NbPts ; ipt ++ )
    {
        ptsneighsV[ ipt ].NbNeigh = 0 ;
    }

    /* Read neighbours of each point */
    nmax = 0 ;
    for( l = 0 ; ( err == STS_OK ) && ( ! feof( Fnei ) ) ; l ++ )
    {
        int     ipt ;

        if ( fscanf( Fnei, "%d", &ipt ) == 1 )
        {
            int nbv ;   /* Declared nb of neighbours */
            int nv ;    /* Effective nb of neighbours */

            if ( fscanf( Fnei, "%d", &nbv ) == 1 )
            {
                NeighT   *neighsV ;
                int      iv ;

                if ( ( neighsV = GenAlloc( nbv, sizeof( NeighT ),
					   0, func, "neighsV" ) ) == NULL )
                {
                    fprintf( stderr, "Can't allocate %d neighb. for pt %d\n" ,
                             nbv, ipt ) ;
                    err = STS_E_MEMORY ;
                    return  err ;
                }

                ptsneighsV[ ipt - 1 ].NeighsV = neighsV ;

                for ( iv = 0, nv = 0 ; 
                      ( iv < nbv ) && ( ! feof( Fnei ) ) ; 
                      iv ++ )
                {
                    int     iptv ;

                    if ( fscanf( Fnei, "%d", &iptv ) == 1 )
                    {
                        if ( ( 1 <= iptv ) && ( iptv <= NbPts ) )
                        {
                            neighsV[ nv ].Index = iptv - 1 ;
                            nv ++ ;
                        }
                    }
                    else
                    {
                        fprintf( stderr, 
                                 "Error in neighb. file l.%d : neighbor %d\n",
                                 l, iv ) ;
                        err = STS_E_FILE ;
                    }
                } /* end  for ( iv ... ) */

                if ( weighted )
                {
                    for ( iv = 0, nv = 0 ; 
                         ( iv < nbv ) && ( ! feof( Fnei ) ) ; 
                         iv ++ )
                    {
                        float   weight ;

                        if ( fscanf( Fnei, "%g", &weight) == 1 )
                        {
                            if ( weight != 0.0 )
                            {
                                neighsV[ nv ].Weight = weight ;
                                nv ++ ;
                            }
                        }
                        else
                        {
                          fprintf( stderr, 
                                  "Error in neighb. file l.%d : weight %d\n",
                                  l, iv ) ;
                            err = STS_E_FILE ;
                        }
                    } /* end  for ( iv ... ) */
                } /* end if weighted */
                else
                {
                    for ( iv = 0 ; iv < nv ; iv ++ )
                        neighsV[ iv ].Weight = 1.0 ;
                }

                ptsneighsV[ ipt - 1 ].NbNeigh = nv ;

                if ( nv > nmax )     nmax = nv ;

            } /* end if ( fscanf( ... , &nbv ) == 1 ) */
            else 
                 {}
        } /* end if ( fscanf( ..., &ipt ) == 1 ) */
        else
            {}
    } /* end for( l = 0 ; ( ! feof( Fnei ) ) ; l ++ ) */

    *MaxNeiP = nmax ;  

    return err ;

}   /* end of ReadPtsNeighs() */


/* ------------------------------------------------------------------- */
static int      ReadImageNeigh
                (
                    FILE        *Fnei,          /* I/O */
                    NeighDataT  *NeighDataP     /* O and allocated */
                ) 
/* ------------------------------------------------------------------- */
{
    StatusET  err = STS_OK ;
    int       dlmin, dlmax, dcmin, dcmax ;
    int       dl, dc ;
    int       tnb ;     /* window size */
    int       nv ;      /* nb of non-zero weights */
    INeighT   *tneighV;  /* to be allocated */

    fscanf( Fnei, "%d %d %d %d", &dlmin, &dlmax, &dcmin, &dcmax ) ;
    tnb = ( dlmax - dlmin + 1 ) * ( dcmax - dcmin + 1 ) ;

    if ( ( tneighV = GenAlloc( tnb, sizeof( INeighT ),
			       0, "ReadImageNeigh", "tneighV" ) ) == NULL )
        return STS_E_MEMORY ;

    NeighDataP->Image.NeighsV = tneighV ;

    nv = 0 ;
    for ( dl = dlmin ; dl <= dlmax ; dl ++ )
    {
        for ( dc = dcmin ; dc <= dcmax ; dc ++ )
        {
            float   weight ;

            if ( fscanf( Fnei, "%g", &weight ) == 1 )
            {
                if ( weight != 0.0 )
                {
                    tneighV[ nv ].Dl = dl ;
                    tneighV[ nv ].Dc = dc ;
                    tneighV[ nv ].Weight = weight ;
                    nv ++ ;
                }
                else /* else this neighbour has no weight : skip it */
                { 
                }
            }
            else
            {
                fprintf( stderr, "Neighbors file error (dl = %d, dc = %d)\n",
                         dl, dc ) ;
                err = STS_E_FILE ;
            }
        } /* for dc */
    } /* for dl */

    NeighDataP->Image.NbNeigh = nv ;

    return err ;

}   /* end of ReadImageNeighs() */


/* ------------------------------------------------------------------- */
static int  SetImageNeigh
            ( 
                NeighET     NeighSpec,          /* I */
                char*       NeiDescS,           /* O [LEN_LINE+1] */
                NeighDataT* NeighDataP          /* O and allocated */
            ) 
/* ------------------------------------------------------------------- */
{
    INeighT   *neighV;  /* to be allocated */

    switch( NeighSpec )
    {
    case NEIGH_FOUR :
        if ( ( neighV = GenAlloc( 4, sizeof( INeighT ),
				  0, "SetImageNeigh", "neighV" ) ) == NULL )
        {
            fprintf( stderr, "Could not allocate %d image neighbours\n", 
                     4 ) ;
            return STS_E_MEMORY ;
        }
        NeighDataP->Image.NeighsV = neighV ;
        NeighDataP->Image.NbNeigh = 4 ;

        neighV[ 0 ].Dl = -1 ;
        neighV[ 0 ].Dc = 0 ;
        neighV[ 0 ].Weight = 1.0 ;

        neighV[ 1 ].Dl = 0 ;
        neighV[ 1 ].Dc = -1 ;
        neighV[ 1 ].Weight = 1.0 ;

        neighV[ 2 ].Dl = 0 ;
        neighV[ 2 ].Dc = +1 ;
        neighV[ 2 ].Weight = 1.0 ;

        neighV[ 3 ].Dl = +1 ;
        neighV[ 3 ].Dc = 0 ;
        neighV[ 3 ].Weight = 1.0 ;

        strncpy( NeiDescS , 
                 "  Default 1st-order neighbors (horizontal and vertical)\n",
                 LEN_LINE ) ;
        break ;

    default :
        fprintf( stderr, "Unknown neighborhood type (%d)\n", NeighSpec ) ;
        return STS_E_FUNCARG ;
    }

    return STS_OK ;
}   /* end of SetNeigh() */


/* ------------------------------------------------------------------- */
static int SaveResults
        ( 
	  const int          Argc,
	  const char*        Argv[],
          const int          Npt,                   /* I */
          const int          Nd,                    /* I */
          const float*       ClassifM,              /* I */
          const SpatialT*    SpatialP,              /* I */
          const NemParaT*    NemParaP,              /* I */
          const StatModelT*  ModelP,                /* I */
	  const CriterT*     CriterP                /* I */ /*V1.03-d*/
        ) 
/* ------------------------------------------------------------------- */
{
    FILE*       fout ;
    FILE*       fmf ;
    char        mfname[ LEN_FILENAME ] ;
    int         iarg ;          /* arg counter   : 0..Argc-1 */
    int         ipt ;           /* point counter : 0..Npt */
    int         k ;             /* class counter : 0..nk-1 */
    int         d ;             /* dimension counter : 0..Nd-1 */
    int         nk = ModelP->Spec.K ;
    StatusET    err = STS_OK ;


    /* SaveResults algorithm :
       - save classification matrix, in hard or fuzzy partition format
       - save means of each class
    */

    if ( strcmp( NemParaP->OutBaseName, "-" ) != 0 )  /*V1.04-b*/
      {
	if ( ( fout = fopen( NemParaP->OutName, "w" ) ) == NULL )
	  {
	    fprintf( stderr, "Could not open file '%s' in write mode\n", 
		     NemParaP->OutName ) ;
	    return STS_E_FILEOUT ;
	  }
      }
    else
      fout = stdout ;

    if ( NemParaP->Format == FORMAT_HARD )
      {
        int*   kmaxesV;  /* classes having same maximum probabilites */

        if ( ( kmaxesV = GenAlloc( nk, sizeof( int ),
				   0, "SaveResults", "kmaxesV" ) ) == NULL ) 
            return STS_E_MEMORY ;

        /* For each record */
        for ( ipt = 0 ; ( ipt < Npt ) && ( err == STS_OK ) ; ipt ++ )
          {
            int   kmax ;

            /* Compute its MAP class */
            kmax = ComputeMAP( ClassifM, ipt, nk, NemParaP->TieRule, 
			       kmaxesV ) ;

            /* Save the record's class (add 1 to have k in 1..Nk) */
            if ( fprintf( fout, "%d ", kmax + 1 ) == EOF )
              {
                fprintf( stderr, "Cannot write Ci i = %d\n",
                            ipt + 1 ) ;
                err = STS_E_FILEOUT ;
              }

            if ( ( SpatialP->Type == TYPE_IMAGE ) &&
                 ( ((ipt + 1) % SpatialP->NeighData.Image.Nc) == 0 ) )
              {
                fprintf( fout, "\n" ) ;
              }
          } /* end for ( ipt ... ) */

	fprintf( fout, "\n" ) ; /*V1.07-b*/

        GenFree( kmaxesV ) ;

      } /* end if Format == HARD */
    else
      {
        for ( ipt = 0 ; ( ipt < Npt ) && ( err == STS_OK ) ; ipt ++ )
          {
            for ( k = 0 ; ( k < nk ) && ( err == STS_OK ) ; k ++ )
              {
                if ( fprintf( fout, " %5.3f ", 
                              ClassifM[ (ipt * nk) + k ] ) == EOF )
                  {
                    fprintf( stderr, "Cannot write Uik i = %d, k =%d\n",
                            ipt + 1, k + 1 ) ;
                    
                    err = STS_E_FILEOUT ;
                  }
              } /* end for ( k = 0 ... ) */
            fprintf( fout, "\n" ) ;
          } /* end for ( ipt ... ) */
      } /* end else Format != HARD */

    if ( fout != stdout )
      fclose( fout ) ;

    if ( strcmp( NemParaP->OutBaseName , "-" ) != 0 )  /*V1.04-b*/
      {
	strncpy( mfname, NemParaP->OutBaseName, LEN_FILENAME ) ;
	strncat( mfname, EXT_MFNAME, LEN_FILENAME ) ;
	if ( ( fmf = fopen( mfname, "w" ) ) == NULL )
	  {
	    fprintf( stderr, "Could not open file '%s' in write mode\n", 
		     mfname ) ;
	    return STS_E_FILEOUT ;
	  }
      }
    else
      fmf = stderr ;

    /*V1.03-f*/
    fprintf( fmf, "Command line : \n\n  " ) ;
    for ( iarg = 0 ; iarg < Argc ; iarg ++ )
      fprintf( fmf, "%s " , Argv[ iarg ] ) ;
    fprintf( fmf, "\n\n" ) ;

    /*V1.03-d*/ /*V1.03-e*/ /*V1.05-j*/
    fprintf( fmf, 
	     "Criteria U=NEM, D=Hathaway, L=mixture, M=markov ps-like, error\n\n" );
    fprintf( fmf, "  %g    %g    %g    %g   %g\n\n", 
	     CriterP->U, CriterP->D, CriterP->L, CriterP->M, 
	     CriterP->Errcur.Errorrate ) ; 

    fprintf( fmf, "Beta (%s)\n", BetaDesVC[ ModelP->Spec.BetaModel ] );
    fprintf( fmf, "  %6.4f\n", ModelP->Para.Beta ) ;

    /* GG Cumulate U criteria in file cumul-u.txt */
    
    FILE *CumFile;
    if ( (CumFile = fopen(  "cumul.u.txt", "a" )) == NULL ) {
         fprintf( stderr, "Cannot write in file %s \n", "cumul-u.txt" ) ;
    } else {
      fprintf( CumFile, " %g \n",  CriterP->U);
      
      close( CumFile );
    }

    /*V1.06-a*/
    switch( ModelP->Spec.ClassFamily )
      {
      case FAMILY_NORMAL:
	fprintf( fmf, "Mu (%d), Pk, and sigma (%d) of the %d classes\n\n", 
		 Nd, Nd, nk ) ;
	break;

      case FAMILY_LAPLACE:
	fprintf( fmf, "Mu (%d), Pk, and lambda (%d) of the %d classes\n\n", 
		 Nd, Nd, nk ) ;
	break;

      case FAMILY_BERNOULLI: /*V1.07-a*/
	fprintf( fmf, "Mu (%d), Pk, and disp (%d) of the %d classes\n\n", 
		 Nd, Nd, nk ) ;
	break;

      default:
	fprintf( fmf, "Mu (%d), Pk, and disp (%d) of the %d classes\n\n", 
		 Nd, Nd, nk ) ;
      }

    for ( k = 0 ; k < nk ; k ++ )
      {
        /* Save mean, proportion, volume, and normalized variance matrix */
        /* V1.03-b*/    /*V1.06-a*/
        for ( d = 0 ; d < Nd ; d ++ )
        {
            fprintf( fmf, " %10.3g ", ModelP->Para.Center_KD[ (k*Nd) + d ] ) ;
        }
        fprintf( fmf, "  %5.3g  ", ModelP->Para.Prop_K[ k ] ) ;
        for ( d = 0 ; d < Nd ; d ++ )
        {
	  switch( ModelP->Spec.ClassFamily )
	    {
	    case FAMILY_NORMAL:
	      fprintf( fmf, " %10g ", 
		       sqrt( ModelP->Para.Disp_KD[ (k*Nd) + d ] ) ) ;
	      break;

	    case FAMILY_LAPLACE:
	      fprintf( fmf, " %10g ", 
		       ModelP->Para.Disp_KD[ (k*Nd) + d ] ) ;
	      break;

	    case FAMILY_BERNOULLI: /*V1.07-a*/
	      fprintf( fmf, " %10g ", 
		       ModelP->Para.Disp_KD[ (k*Nd) + d ] ) ;
	      break;

	    default:
	      fprintf( fmf, " %10g ", 
		       ModelP->Para.Disp_KD[ (k*Nd) + d ] ) ;
	    }
        }
        fprintf( fmf, "\n" ) ;
      }

    if ( fmf != stderr )
      fclose( fmf ) ;

    
    return err ;

}  /* end of SaveResults() */



/* ------------------------------------------------------------------- */
static void FreeAllocatedData
(
 DataT*       DataP,       /* O and deallocated */
 SpatialT*    SpatialP,    /* O and deallocated */
 ModelParaT*  ModelParaP,  /* O and deallocated */  /*V1.06-a*/
 CriterT*     CriterP,     /* O and deallocated */  /*V1.06-d*/
 float*       ClassifM     /* O and deallocated */
)
/* ------------------------------------------------------------------- */
{
  int ipt ;  /* counter 0..Npt-1 to free each point neighbors */


  /* Free components of DataP */
  GenFree( DataP->PointsM ) ;     DataP->PointsM    = NULL ;
  GenFree( DataP->LabelV ) ;      DataP->LabelV     = NULL ;
  GenFree( DataP->SiteVisitV ) ;  DataP->SiteVisitV = NULL ;
  GenFree( DataP->SortPos_ND ) ;  DataP->SortPos_ND = NULL ;

  /* Free components of SpatialP */
  switch( SpatialP->Type )
    {
    case TYPE_SPATIAL: 
      /* deallocate each point's neighbors */
      for ( ipt = 0; ipt < DataP->NbPts ; ipt ++ )
	GenFree( SpatialP->NeighData.PtsNeighsV[ ipt ].NeighsV ) ;

      /* deallocate array of array of neighbors */
      GenFree( SpatialP->NeighData.PtsNeighsV ) ;

      break ;

    case TYPE_IMAGE:
      /* deallocate neighborhood window */
      GenFree( SpatialP->NeighData.Image.NeighsV ) ;
      break ;

    default: /* non spatial : no allocated spatial info */
      ;
    }

  /* Free components of ModelParaP */
  GenFree( ModelParaP->Center_KD ) ;
  GenFree( ModelParaP->Disp_KD ) ;
  GenFree( ModelParaP->Prop_K ) ;

  GenFree( ModelParaP->NbObs_K ) ;
  GenFree( ModelParaP->NbObs_KD ) ;
  GenFree( ModelParaP->Iner_KD ) ;


  /* Free components of CriterP */
  GenFree( CriterP->Errinfo.Refclas_N_Kr  ) ;
  GenFree( CriterP->Errinfo.Perm_Kmfac_Km ) ;
  GenFree( CriterP->Errcur.Agree_Km_Km    ) ;
  GenFree( CriterP->Errcur.Loclas_N_Kc    ) ;

  CriterP->Errinfo.Refclas_N_Kr  = NULL ;
  CriterP->Errinfo.Perm_Kmfac_Km = NULL ;
  CriterP->Errcur.Agree_Km_Km    = NULL ;
  CriterP->Errcur.Loclas_N_Kc    = NULL ;

  /* Free classification matrix */
  GenFree( ClassifM ) ;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
