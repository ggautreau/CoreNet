./._nem_alg.h                                                                                       000755  000765  000765  00000000312 11541743323 012657  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_alg.h                                                                                           000755  000765  000024  00000002151 11541743323 012600  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         #include "nem_typ.h"    /* DataT, ... */

/*
  V1.06-a   17-JUN-1998  NoiseModel -> StatModel
  V1.06-b   30-JUN-1998  ClassifyByNem const NemParaP, no const DataP (sort)
  V1.06-c   20-SEP-1998  Add LabelToClassVector
  V1.06-d   21-SEP-1998  Add TieRule arg to ComputeMAP
*/

    int ClassifyByNem
        ( 
          const NemParaT      *NemParaP,        /* I */
          const SpatialT      *SpatialP,        /* I */
          DataT               *DataP,           /* I/O */
          StatModelT          *StatModelP,      /* I/O */
          float               *ClassifM,        /* I/O */
	  CriterT             *CriterP          /* O */  /*V1.03-f*/
        ) ;

    int ComputeMAP
        ( 
          const float* ClassifM,   /* I */
          int          Ipt,        /* I */
          int          Nk,         /* I */
	  TieET        TieRule,    /* I */ /*V1.06-d*/
          int*         kmaxesV     /* T [Nk] */
        ) ;


void LabelToClassVector
 ( 
  const int Nk,    /* I: number of classes */
  const int Label, /* I: class of interest 0..Nk-1 */
  float* Cout_K    /* O: classification vector [Nk] */
 ) ;
                                                                                                                                                                                                                                                                                                                                                                                                                       ./._nem_arg.c                                                                                       000755  000765  000765  00000000312 11541743324 012661  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_arg.c                                                                                           000755  000765  000024  00000125357 11541743324 012620  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\

    NEM_ARG.C

    Programme NEM (Neighborhood EM) : Traitement des arguments

    Van Mo DANG       Janvier 96

    Steps to add an option :
    1 - check unused option characters by running 'nem_exe' without args
    2 - add option fetching in NemArgs()
    3 - add field to NemParaT in "nem_typ.h"  : OptMeaning
    4 - add eventual enum in "nem_typ.h"      : OptET
    5 - add eventual string array             : OptStrVC
    6 - add option default value              : DEFAULT_OPT
    7 - add comment to this header            : 1.06-k blabla
    8 - add short help in PrintUsage()        : -o [%s ] blabla {%s|%s}
    9 - add long help/vers in "nem_user.txt"

Vers-mod  Date         Who Description

1.03-a    22-AUG-1997  MD  '-l' to output logfile (default now is NO log)
1.04-a    04-OCT-1997  MD  '-s f myfile.zz' instead of just '-s f'
1.04-b    05-OCT-1997  MD  '-B heu <bstep>' implemented
1.04-c    09-OCT-1997  MD  if "-o -" log name is input base name
1.04-d    10-OCT-1997  MD  '-c cvthres' implemented
1.04-e    10-OCT-1997  MD  '-S seed' implemented
1.04-f    10-OCT-1997  MD  '-O order' implemented
1.04-g    13-OCT-1997  MD  '-D' implemented
1.04-h    13-OCT-1997  MD  '-R' implemented (reference class)
1.04-i    13-OCT-1997  MD  '-C crit'
1.04-j    04-NOV-1997  MD  '-I eiter'
1.04-k    11-JAN-1997  MD  now '-h general/filein/fileout/examples/vers'
1.05-a    12-JAN-1998  MD  '-M miss' implemented
1.05-b    05-FEB-1998  MD  add 'M' in '-C crit'
1.06-a    17-JUN-1998  MD  ReadStrFile and alloc parameters <- from nem_exe.c
1.06-b    23-JUN-1998  MD  use new StatModel structure (old NoiseModel)
1.06-b    01-AUG-1998  MD  model switch now '-m norm/lapl p* s**'
1.06-c    03-AUG-1998  MD  add '-s mi' and '-s mf'
1.06-d    10-SEP-1998  MD  add '-U para|seq'
1.06-e    10-SEP-1998  MD  add '-t random|first'
1.06-f    15-SEP-1998  MD  add '-G nbit conv step'
1.06-g    21-SEP-1998  MD  change '-c 0.01' to '-c {none|[clas]|crit} 0.01'
1.06-h    05-OCT-1998  MD  default nb of random inits = D * 5 
1.06-i    05-OCT-1998  MD  default crit = CRIT_M
1.07-a    08-APR-1999  MD  Add Bernoulli family
\*/

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
#define DEFAULT_CRIT         CRIT_M             /*V1.04-h*/
#define DEFAULT_CVTHRES      0.01               /*V1.04-d*/
#define DEFAULT_CVTEST       CVTEST_CLAS        /*V1.06-g*/
#define DEFAULT_FAMILY       FAMILY_NORMAL      /*V1.06-b*/
#define DEFAULT_DISPER       DISPER___          /*V1.06-b*/
#define DEFAULT_PROPOR       PROPOR__           /*V1.06-b*/
#define DEFAULT_FORMAT       FORMAT_HARD
#define DEFAULT_INIT         INIT_SORT
#define DEFAULT_MISSING      MISSING_REPLACE    /*V1.05-a*/
#define DEFAULT_SORTEDVAR    0
#define DEFAULT_NEIGHSPEC    NEIGH_FOUR
#define DEFAULT_NBITERS      100
#define DEFAULT_NBEITERS     1
#define DEFAULT_NBRANDINITS  5                  /*V1.06-h*/
#define DEFAULT_ORDER        ORDER_DIRECT       /*V1.04-f*/
#define DEFAULT_UPDATE       UPDATE_SEQ         /*V1.06-d*/
#define DEFAULT_TIE          TIE_RANDOM         /*V1.06-e*/


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


static const char *TypeStrC[ TYPE_NB ] = { "S", "I", "N" } ; /*V1.06-a*/

static const char *AlgoStrVC[ ALGO_NB ] = { "nem" , "ncem" , "gem" } ;
static const char *BtaStrVC[ BETA_NB ] = { "fix" , "psgrad", 
					   "heu_d" , "heu_l" } ;
       const char *CritStrVC[ CRIT_NB ] = { "U" , "M", "D" , "L" } ;
static const char *CvTestStrVC[ CVTEST_NB ] = { "none", "clas" , "crit" } ;
static const char *FormatStrVC[ FORMAT_NB ] = { "hard", "fuzzy" } ;
static const char *HelpStrVC[ HELP_NB ] = {
  "general",
  "options",
  "examples",
  "filein",
  "fileout",
  "versions"
} ; /* V1.04-k*/

static const char *InitStrVC[ INIT_NB ] = { "s", "r", "mi", "mf", "f", "l" } ;
/*V1.05-a*/
static const char *MissStrVC[ MISSING_NB ] = { "replace" , "ignore" } ;
/*V1.06-b*/
static const char *FamilyStrVC[ FAMILY_NB ] = { "norm", "lapl", "bern" } ;
static const char *DisperStrVC[ DISPER_NB ] = { "s__", "sk_", "s_d", "skd" } ;
static const char *ProporStrVC[ PROPOR_NB ] = { "p_", "pk" } ;
static const char *NeighStrVC[ NEIGH_NB ] = { "4", "f" } ;
static const char *OrderStrVC[ ORDER_NB ] = { "direct", "random" } ;/*V1.04-f*/
static const char *UpdateStrVC[ UPDATE_NB ] = { "seq", "para" } ;/*V1.06-d*/
static const char *TieStrVC[ TIE_NB ] = { "random", "first" } ;/*V1.06-e*/


  static int  ReadStrFile  /*V1.06-a*/
         (
             const char*  BaseName,     /* I */
             char*        Comment,      /* O */
             DataT*       DataP,        /* O */
             ImageNeighT* ImageP,       /* O */
             TypeET*      TypeP         /* O */
         );                             

  static char GetOptionSwitch( const char *Arg ) ;
  static int  GetEnum( const char* S, const char* SV[], int SizeV );/*V1.04-b*/
  static StatusET PrintHelp( const char* HelpType ) ;
  static void my_strupr( char *s ) ;




/* ==================== LOCAL FUNCTION PROTOTYPING =================== */

/* Called by NemArgs */

StatusET GetMixPara        /* STS_OK or STS_E_ARG */
 (
  const FamilyET  Family,     /* I : distribution family */
  const int       K,          /* I : number of classes */
  const int       D,          /* I : number of variables */
  const char*     Opts_Q[],   /* I : array of options */
  const int       Q,          /* I : number of options */
  int*            IoptP,      /* I/O : index of current option */
  ModelParaT*     ParaP       /* O : read parameters */
 ) ;



/* ==================== GLOBAL FUNCTION DEFINITION =================== */


/* ------------------------------------------------------------------- */
int NemArgs
        (
          int           Argc,
          const char    *Argv[],
          char          *Fname,         /* O */
          StatModelT    *StatModelP,    /* O */
          NemParaT      *NemParaP,      /* O */
	  char*         DatadescS,      /* O */
	  DataT*        DataP,          /* O */
	  ImageNeighT*  ImageP,         /* O */
	  TypeET*       SpatialTypeP    /* O */	  
        ) 
/* ------------------------------------------------------------------- */
{
  StatusET     err = STS_OK ;
  const char*  func = "NemArgs" ;
  int          nk ;
  int          nbopt ;
  int          iopt ;
  const char** opts;        

  switch ( Argc ) /*V1.04-k*/
  {
  case 1:  /* No args : print usage */
    PrintUsage( Argv[ 0 ] ) ;
    return STS_I_DONE ;

  case 2:  /* 1 arg : only possibility is -v to ask version information */
    if ( GetOptionSwitch( Argv[ 1 ] ) == 'v' )
      {
	PrintVersions( stdout ) ;
	return STS_I_DONE ;
      }
    else
      return STS_E_ARG ;

  case 3:  /* 2 args : -h help or file nk */
    switch( GetOptionSwitch( Argv[ 1 ] ) )
      {
      case '?':
      case 'h':
	return PrintHelp( Argv[2] ) ;

      default: /* continue processing args */
             break;
      }

  default:  /* 3 or more args : continue processing args */
    break ;
  }

  strncpy( Fname, Argv[1], LEN_FILENAME ) ;
  sscanf( Argv[2], "%d", &nk ) ;
  StatModelP->Spec.K = nk ;
  if ( nk <= 0 )
  {
      fprintf( stderr, "Nb of classes must be > 0 (here %s)\n",
               Argv[2] ) ;
      return STS_E_ARG ;
  }


  /* !!! Read structure file */ /*V1.06-a*/
  if ( ( err = ReadStrFile( Fname, 
			    DatadescS,
			    DataP, 
			    ImageP,
			    SpatialTypeP ) ) != STS_OK )
    return err ;

  /* !!! Allocate model parameters */ /*V1.06-a*/
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
  StatModelP->Spec.ClassFamily = DEFAULT_FAMILY ;
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
  NemParaP->NbRandomInits = DEFAULT_NBRANDINITS*DataP->NbVars ;  /*V1.06-h*/
  NemParaP->Seed          = time( NULL ) ;          /*V1.04-e*/
  NemParaP->Format        = DEFAULT_FORMAT ;
  NemParaP->InitMode      = DEFAULT_INIT ;
  NemParaP->SortedVar     = DEFAULT_SORTEDVAR ;
  NemParaP->NeighSpec     = DEFAULT_NEIGHSPEC ;
  NemParaP->VisitOrder    = DEFAULT_ORDER ;         /*V1.04-f*/
  NemParaP->SiteUpdate    = DEFAULT_UPDATE ;        /*V1.06-d*/
  NemParaP->TieRule       = DEFAULT_TIE ;           /*V1.06-e*/
  NemParaP->Debug         = FALSE ;                 /*V1.04-g*/
  strncpy( NemParaP->OutBaseName, Fname, LEN_FILENAME ) ;
  strncpy( NemParaP->NeighName, Fname, LEN_FILENAME ) ;
  strncat( NemParaP->NeighName, ".nei", LEN_FILENAME ) ;
  strncpy( NemParaP->RefName, "", LEN_FILENAME ) ;
  
  /* Treat options */
  nbopt = Argc - 3 ;
  opts  = & ( Argv[ 3 ] ) ;

  for ( iopt = 0, err = STS_OK ; 
        ( iopt < nbopt ) && ( err == STS_OK ) ;
        iopt ++ )
    {
      char copt = GetOptionSwitch( opts[ iopt ] ) ;

      if ( copt != '\0' ) 
        switch( copt )
        {
        case 'a' :    /* next arg is type of algorithm */
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->Algo = GetEnum( opts[ iopt ], AlgoStrVC, ALGO_NB ) ;
              if ( NemParaP->Algo == -1 )
                {
                  fprintf( stderr, " Unknown type of algorithm %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'b' :    /* next arg is beta */
          iopt ++ ;
          if ( iopt < nbopt )
            {
              StatModelP->Para.Beta = atof( opts[ iopt ] ) ;
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'c' :    /* next arg : conv test */  /*V1.04-d*/ /*V1.06-g*/
	  iopt ++ ;
          if ( iopt < nbopt ) /* conv test given */ {
	    NemParaP->CvTest=GetEnum( opts[ iopt ], CvTestStrVC, CVTEST_NB );
	    if ( NemParaP->CvTest == -1 ) {
	      fprintf( stderr, " Unknown convergence test %s\n", opts[iopt] ) ;
	      err = STS_E_ARG ;
	    }
	    else if ( NemParaP->CvTest != CVTEST_NONE ) /* get threshold */ {
	      iopt ++ ;
	      if ( iopt < nbopt ) /* threshold given */ {
		NemParaP->CvThres = atof( opts[ iopt ] ) ;
		if ( NemParaP->CvThres <= 0 ) {
		  fprintf( stderr, " Conv threshold must be > 0 (here %s)\n",
			   opts[ iopt ] ) ;
		  err = STS_E_ARG ;
		} /* else threshold > 0 : OK */
	      } 
	      else /* threshold not given */ {
		fprintf( stderr, " Expected threshold for conv test %s\n", 
			 opts[ iopt - 1 ] ) ;
		err = STS_E_ARG ;
	      }
            } /* threshold checked */
	  } /* conv test given */
          else {
	    fprintf( stderr, " Expected value after switch %s\n", 
		     opts[ iopt - 1 ] ) ;
	    err = STS_E_ARG ;
	  }
          break ;


        case 'f' :    /* next arg is partition input/output format */
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->Format=GetEnum( opts[ iopt ], FormatStrVC, FORMAT_NB );
              if ( NemParaP->Format == -1 )
                {
                  fprintf( stderr, " Unknown format %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'i' :    /* next arg is number of iterations */
          iopt ++ ;
          if ( iopt < nbopt )
            {
              NemParaP->NbIters = atoi( opts[ iopt ] ) ;
              if ( NemParaP->NbIters < 0 )
                {
                  fprintf( stderr, "Nb iterations must be >= 0 (here %s)\n",
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'l' :    /* next arg is 'y' if log requested */  /*V1.03-a*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              if ( opts[ iopt ][ 0 ] == 'y' )
      	  NemParaP->DoLog = TRUE ;
      	else
      	  NemParaP->DoLog = FALSE ;
            }
          else
            {
              fprintf( stderr, " Expected y or n after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'm' :    /* next 3 args is stat model */ /*V1.06-b*/
          if ( iopt + 3 < nbopt )
            {
	      iopt ++ ;
	      StatModelP->Spec.ClassFamily = 
		GetEnum( opts[ iopt ], FamilyStrVC, FAMILY_NB );
              if ( StatModelP->Spec.ClassFamily == -1 )
                {
                  fprintf( stderr, " Unknown family %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
	      iopt ++ ;
	      StatModelP->Spec.ClassPropor = 
		GetEnum( opts[ iopt ], ProporStrVC, PROPOR_NB );
              if ( StatModelP->Spec.ClassPropor == -1 )
                {
                  fprintf( stderr, " Unknown proportion %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
	      iopt ++ ;
	      StatModelP->Spec.ClassDisper = 
		GetEnum( opts[ iopt ], DisperStrVC, DISPER_NB );
              if ( StatModelP->Spec.ClassDisper == -1 )
                {
                  fprintf( stderr, " Unknown dispersion %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected 3 values after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'n' :    /* next arg is neighborhood specification */
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->NeighSpec=GetEnum( opts[ iopt ], NeighStrVC, NEIGH_NB);
              if ( NemParaP->NeighSpec == -1 )
                {
                  fprintf( stderr, " Unknown neighborhood specification %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'o' :    /* next arg is output file name */
          iopt ++ ;
          if ( iopt < nbopt )
            {
              strncpy( NemParaP->OutBaseName, opts[ iopt ], LEN_FILENAME ) ;
            }
          else
            {
              fprintf( stderr, " Expected file name after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 's' :    /* next arg is initialization mode */ 
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->InitMode = GetEnum( opts[ iopt ], InitStrVC, INIT_NB );

              switch( NemParaP->InitMode )  {

                case INIT_SORT:  /* next arg is sorted var */
                  iopt ++ ;
                  if ( iopt < nbopt )
                    {
                      NemParaP->SortedVar = atoi( opts[ iopt ] ) - 1 ;
                    } 
                  else
                    {
                      fprintf( stderr, 
                               " Expected variable index after switch %s\n", 
                               opts[ iopt - 1 ] ) ;
                      err = STS_E_ARG ;
                    }
                  break ;

                case INIT_RANDOM:  /* next arg is number of initializations */
                  iopt ++ ;
                  if ( iopt < nbopt )
                    {
                      NemParaP->NbRandomInits = atoi( opts[ iopt ] ) ;
                    } 
                  else
                    {
                      fprintf( stderr, 
                               " Expected nb of trials after switch %s\n", 
                               opts[ iopt - 1 ] ) ;
                      err = STS_E_ARG ;
                    }
                  break ;

                case INIT_FILE:  /* next arg is name of initial partition */
		  iopt ++ ;   /*V1.04-a*/
		  if ( iopt < nbopt )
		    {
		      strncpy( NemParaP->StartName, opts[ iopt ], 
			       LEN_FILENAME ) ;
		    }
		  else
		    {
		      fprintf( stderr, 
			       " Expected file name after switch %s\n", 
			       opts[ iopt - 1 ] ) ;
		      err = STS_E_ARG ;
		    }
                  break ;

		case INIT_LABEL:/* next arg is name of partial labels */
		  iopt ++ ;
		  if ( iopt < nbopt )
		    {
		      strncpy( NemParaP->LabelName, opts[ iopt ], 
			       LEN_FILENAME ) ;
		    }
		  else
		    {
		      fprintf( stderr, 
			       " Expected file name after switch %s\n", 
			       opts[ iopt - 1 ] ) ;
		      err = STS_E_ARG ;
		    }
		  break ;

		case INIT_MIXINI: /* next args are parameter values */
		case INIT_MIXFIX: /* next args are parameter values */
		  /*V1.06-c*/
		  err = GetMixPara( StatModelP->Spec.ClassFamily, nk, 
				    DataP->NbVars, opts, nbopt, 
				    & iopt, & StatModelP->Para ) ;
		  break ;

                default:
                  fprintf( stderr, " Unknown init mode %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected init mode after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 't' :    /* next arg is MAP tie rule */  /*V1.06-e*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->TieRule=GetEnum(opts[ iopt ],TieStrVC,TIE_NB);
              if ( NemParaP->TieRule == -1 )
                {
                  fprintf( stderr, " Unknown tie rule %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



        case 'B' :    /* next arg is beta estimation method */ /*V1.04-b*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              StatModelP->Spec.BetaModel = 
		GetEnum( opts[ iopt ], BtaStrVC, BETA_NB );
	      if ( StatModelP->Spec.BetaModel == -1 )
		{
                  fprintf( stderr, " Unknown type of beta estimation %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
		}
	      else {
		if ( ( StatModelP->Spec.BetaModel == BETA_HEUL ) |
		     ( StatModelP->Spec.BetaModel == BETA_HEUD ) )
		  NemParaP->Crit = CRIT_U ;  /* force to choose U criterion */
	      }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



        case 'C' :    /* next arg is crit to choose local max */ /*V1.04-i*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              NemParaP->Crit = GetEnum( opts[ iopt ], CritStrVC, CRIT_NB );
	      if ( NemParaP->Crit == -1 )
		{
                  fprintf( stderr, " Unknown criterion name %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
		}
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



	case 'G' :   /* next 4 args are parameters of beta gradient */
	  if ( iopt + 4 < nbopt )
	    {
	      NemParaP->BtaPsGrad.NbIter    = atof( opts[ iopt + 1 ] ) ;
	      NemParaP->BtaPsGrad.ConvThres = atof( opts[ iopt + 2 ] ) ;
	      NemParaP->BtaPsGrad.Step      = atof( opts[ iopt + 3 ] ) ;
	      NemParaP->BtaPsGrad.RandInit  = atoi( opts[ iopt + 4 ] ) ;
	      iopt += 4 ;
	    }
	  else
	    {
	      fprintf( stderr, " Expected 4 values after switch %s\n", 
		       opts[ iopt - 1 ] ) ;
	      err = STS_E_ARG ;
	    }
	  break ;




	case 'H' :   /* next 5 args are parameters of beta heuristic */
	  if ( iopt + 5 < nbopt )
	    {
	      NemParaP->BtaHeuStep = atof( opts[ iopt + 1 ] ) ;
	      NemParaP->BtaHeuMax  = atof( opts[ iopt + 2 ] ) ;
	      NemParaP->BtaHeuDDrop = atof( opts[ iopt + 3 ] ) ;
	      NemParaP->BtaHeuDLoss = atof( opts[ iopt + 4 ] ) ;
	      NemParaP->BtaHeuLLoss = atof( opts[ iopt + 5 ] ) ;
	      iopt += 5 ;
	    }
	  else
	    {
	      fprintf( stderr, " Expected 4 values after switch %s\n", 
		       opts[ iopt - 1 ] ) ;
	      err = STS_E_ARG ;
	    }
	  break ;


        case 'I' :    /* next arg is number of E-step iterations */ /*V1.04-j*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              NemParaP->NbEIters = atoi( opts[ iopt ] ) ;
              if ( NemParaP->NbEIters < 0 )
                {
                  fprintf( stderr, "Nb iterations must be >= 0 (here %s)\n",
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'M' :    /* next arg is missing data mode */     /*V1.05-a*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->MissMode=GetEnum(opts[ iopt ], MissStrVC, MISSING_NB);
              if ( NemParaP->MissMode == -1 )
                {
                  fprintf( stderr, " Unknown missing mode %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



        case 'O' :    /* next arg is visit order type */     /*V1.04-f*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->VisitOrder=GetEnum(opts[ iopt ], OrderStrVC, ORDER_NB);
              if ( NemParaP->VisitOrder == -1 )
                {
                  fprintf( stderr, " Unknown visit order type %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



        case 'R' :    /* next arg is reference class file */  /*V1.04-h*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              strncpy( NemParaP->RefName, opts[ iopt ], LEN_FILENAME ) ;
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;



        case 'S' :    /* next arg is random generator seed */  /*V1.04-e*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
              NemParaP->Seed = atol( opts[ iopt ] ) ;
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'T' :    /* set debug mode to be true */  /*V1.04-g*/
	  NemParaP->Debug = TRUE ;
          break ;


        case 'U' :    /* next arg is site update mode */  /*V1.06-d*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      NemParaP->SiteUpdate=GetEnum(opts[ iopt ],UpdateStrVC,UPDATE_NB);
              if ( NemParaP->SiteUpdate == -1 )
                {
                  fprintf( stderr, " Unknown site update mode %s\n", 
                           opts[ iopt ] ) ;
                  err = STS_E_ARG ;
                }
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;


        case 'v' :    /* print out version information */
          PrintVersions( stdout ) ;
          err = STS_I_DONE ;
          break ;


        case '?' :    /* long help requested */
        case 'h' :    /*V1.04-k*/
          iopt ++ ;
          if ( iopt < nbopt )
            {
	      err = PrintHelp( opts[ iopt ] ) ;
            }
          else
            {
              fprintf( stderr, " Expected value after switch %s\n", 
                       opts[ iopt - 1 ] ) ;
              err = STS_E_ARG ;
            }
          break ;

        default :
          fprintf( stderr, "Unknown option switch '%s'\n", opts[ iopt ] ) ;
          err = STS_E_ARG ;
        } /* end if ( copt != '\0' ) switch ( copt ) */
      else
        {
          fprintf( stderr, " Expected an option switch at arg %d (%s)\n",
                   iopt + 3 , opts[ iopt ] ) ;
          err = STS_E_ARG ;
        }
    } /* end for ( iopt ... ) */

  strncpy( NemParaP->OutName, NemParaP->OutBaseName, LEN_FILENAME ) ;
  strncat( NemParaP->OutName, 
           NemParaP->Format == FORMAT_HARD ? EXT_OUTNAMEHARD : EXT_OUTNAMEFUZZY,
           LEN_FILENAME ) ;

  if ( NemParaP->DoLog )
  {
    if ( strcmp( NemParaP->OutBaseName , "-" ) != 0 ) /*V1.04-c*/
      {
	strncpy( NemParaP->LogName, NemParaP->OutBaseName, LEN_FILENAME ) ;
	strncat( NemParaP->LogName, EXT_LOGNAME, LEN_FILENAME ) ;
      }
    else
      {
	strncpy( NemParaP->LogName, Fname , LEN_FILENAME ) ;
	strncat( NemParaP->LogName, EXT_LOGNAME, LEN_FILENAME ) ;	
      }
  }
  else
  {
      strcpy( NemParaP->LogName, "" ) ;
  }

  return err ;

} /* end of NemArgs() */



/* ==================== LOCAL FUNCTION DEFINITION =================== */



/* ------------------------------------------------------------------- */
static int  ReadStrFile
            (
                const char*  BaseName,          /* I */
                char*        CommentS,          /* O */
                DataT*       DataP,             /* O */
                ImageNeighT* ImageP,            /* O */
                TypeET*      TypeP              /* O */
            )
/*\
  Read structure file
       expected format : 
         type size nbvars
           type : N | I | S
           size : (if type == I) lx ly (if type == N or S) nbpts
           nbvars : integer
\*/
/* ------------------------------------------------------------------- */
{
    char    infname[ LEN_FILENAME + 1 ] ;
    FILE    *fstr ;
    int     nbelts ;
    char    type[100] ;  /* N | I | S */
    int     n1, n2, n3 ;


    strncpy( infname, BaseName, LEN_FILENAME ) ;
    strncat( infname, ".str" , LEN_FILENAME ) ;

    if ( ReadOpeningComments( infname , "#" , LEN_LINE , 
                              & fstr , CommentS ) == -1 )
    {
        fprintf( stderr, "File %s does not exist\n", infname ) ;
        return STS_E_FILEIN ;
    }


    /* Now we expect a line having structure : "I nl nc d" or "S n d" or
       "N n d"
    */
    nbelts = fscanf( fstr, "%s %d %d %d", type, &n1, &n2, &n3 ) ;
    fclose( fstr ) ;  /*V1.05-h*/

    if ( nbelts < 3 )
    {
        fprintf( stderr, "Structure file (%s) not enough fields\n",
                 infname ) ;
        return STS_E_FILE ;        
    }

    my_strupr( type ) ;

    if ( ! strcmp( type, TypeStrC[ TYPE_NONSPATIAL ] ) )
    {
        /* No spatial information */
        DataP->NbPts  = n1 ;
        DataP->NbVars = n2 ;
        *TypeP        = TYPE_NONSPATIAL ;
    }
    else if ( ! strcmp( type, TypeStrC[ TYPE_SPATIAL ] ) )
    {
        /* Spatial irregular */
        DataP->NbPts  = n1 ;
        DataP->NbVars = n2 ;
        *TypeP        = TYPE_SPATIAL ;
    }
    else if ( ! strcmp( type, TypeStrC[ TYPE_IMAGE ] ) )
    {
        /* Image */
        if ( nbelts < 4 )
        {
            fprintf( stderr, "Structure file (%s) not enough fields\n",
                     infname ) ;
            return STS_E_FILE ;        
        }             
        ImageP->Nl    = n1 ;
        ImageP->Nc    = n2 ;
        DataP->NbPts  = ImageP->Nl * ImageP->Nc ;
        DataP->NbVars = n3 ;
        *TypeP = TYPE_IMAGE ;
    }
    else
    {
        /* default : Error in structure file */
        fprintf( stderr, "Data type %s unknown in file %s\n", 
                 type, infname ) ;
        return STS_E_FILE ;
    }

    return STS_OK ;

}  /* end of ReadStrFile() */




/* ------------------------------------------------------------------- */
static char GetOptionSwitch( const char *Arg ) 
/* ------------------------------------------------------------------- */
{
    if ( ( Arg[ 0 ] == '-' ) || ( Arg[ 0 ] == '/' ) )
      {
        return Arg[ 1 ] ;
      }
    else
      {
        return '\0' ;
      }
} /* end of GetOptionSwitch() */


/* ------------------------------------------------------------------- */
static int GetEnum( const char* S, const char* SV[], int SizeV ) /*V1.04-b*/
/*\

    Returns the position of string S in the array SV consisting of
    SizeV strings (indiced from 0 to SizeV - 1).  Returns -1 if S is
    not found in SV.

\*/
/* ------------------------------------------------------------------- */
{
  int pos ;
  int found ;

  for ( pos = 0 , found = FALSE ;
	( pos < SizeV ) && ( ! found ) ;
	pos ++ )
    {
      if ( ! strcmp( S , SV[ pos ] ) )
	found = TRUE ;
    }

  if ( found )
    return ( pos - 1 ) ;  /* '- 1' because for() adds 1 after finding */
  else
    return -1 ;
}





/* ------------------------------------------------------------------- */
StatusET GetMixPara        /* STS_OK or STS_E_ARG */
 (
  const FamilyET  Family,     /* I : distribution family */
  const int       K,          /* I : number of classes */
  const int       D,          /* I : number of variables */
  const char*     Opts_Q[],   /* I : array of options */
  const int       Q,          /* I : number of options */
  int*            IoptP,      /* I/O : index of current option */
  ModelParaT*     ParaP       /* O : read parameters */
 )
/* ------------------------------------------------------------------- */
{
  StatusET  sts = STS_OK ;
  int       k ;   /* class counter : 0..K-1 */
  int       d ;   /* variable counter : 0..D-1 */
  float     pK ;  /* remaining proportion for class K */
	    
  int       nbpara = K - 1 + K * D + K * D ;

  if ( (*IoptP) + nbpara >= Q ) {
    fprintf( stderr, "Need at least %d parameters (%d provided)\n",
	     nbpara, Q - (*IoptP) ) ;

    return STS_E_ARG ;
  }

  /* Read proportions */
  for ( k = 0, pK = 1 ; 
	k < K - 1 ; 
	k ++ ) {
    (*IoptP) ++ ;
    ParaP->Prop_K[ k ] = atof( Opts_Q[ (*IoptP) ] ) ;
    pK = pK - ParaP->Prop_K[ k ] ;
  }
  ParaP->Prop_K[ K - 1 ] = pK ;
  if ( pK <= 0.0 ) {
    fprintf( stderr, "Last class has pK = %5.2f <= 0\n", pK ) ;
    return STS_E_ARG ;
  }


  /* Read centers */
  for ( k = 0 ; k < K ; k ++ ) {

    for ( d = 0 ; d < D ; d ++ ) {

      (*IoptP) ++ ;
      ParaP->Center_KD[ k * D + d ] = atof( Opts_Q[ (*IoptP) ] ) ;
    }
  }

  /* Read dispersions */
  for ( k = 0 ; k < K ; k ++ ) {

    for ( d = 0 ; d < D ; d ++ ) {

      (*IoptP) ++ ;
      if ( Family == FAMILY_NORMAL ) {
	float x = atof( Opts_Q[ (*IoptP) ] ) ;
	ParaP->Disp_KD[ k * D + d ] = x * x ;
      }
      else
	ParaP->Disp_KD[ k * D + d ] = atof( Opts_Q[ (*IoptP) ] ) ;

      if ( ParaP->Disp_KD[ k * D + d ] <= 0 ) {

	fprintf( stderr, "Dispersion(k=%d, d=%d) = %5.3f <= 0\n", 
		 k+1, d+1, ParaP->Disp_KD[ k * D + d ] ) ;
	sts = STS_E_ARG ;
      }
    }
  }

  return sts ;

}   /* end of GetMixPara() */



/* ------------------------------------------------------------------- */
void PrintUsage( const char* CmdName )
/* ------------------------------------------------------------------- */
{
    fprintf( stderr, "Usage :   " ) ;
    fprintf( stderr, "  %s   file  nk  [ option1 option2 ... ]\n\n",
            CmdName ) ;
    fprintf( stderr, "  file       base name of input and output files\n" ) ;
    fprintf( stderr, "  nk         number of classes\n\n" ) ;

    fprintf( stderr, "Options :   [ default  ]\n\n" ) ;

    fprintf( stderr, "  -a algo   [ %-8s ]   type of algorithm { ",
            AlgoStrVC[ DEFAULT_ALGO ] ) ;
    {
      AlgoET ialgo ;
      
      for ( ialgo = 0 ; ialgo < ALGO_NB ; ialgo ++ )
        fprintf( stderr, "%s ", AlgoStrVC[ ialgo ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -b beta   [ %-8g ]   spatial coefficient, >= 0.0\n",
            DEFAULT_BETA ) ;

    fprintf( stderr, "  -c wh thr [%4s %5.3g]   convergence test { ",
            CvTestStrVC[ DEFAULT_CVTEST ], DEFAULT_CVTHRES ) ;
    {
      CvemET icvtest ;
      
      for ( icvtest = 0 ; icvtest < CVTEST_NB ; icvtest ++ )
        fprintf( stderr, "%s ", CvTestStrVC[ icvtest ] ) ;
    }
    fprintf( stderr, "} + threshold (>0)\n" ) ; /*V1.04-d*//*V1.06-g*/

    fprintf( stderr, "  -f format [ %-8s ]   partition format, { ",
            FormatStrVC[ DEFAULT_FORMAT ] ) ;
    {
      FormET ipart ;
      
      for ( ipart = 0 ; ipart < FORMAT_NB ; ipart ++ )
        fprintf( stderr, "%s ", FormatStrVC[ ipart ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -i itmax  [ %-8d ]   number of NEM iterations (>= 0)\n",
             DEFAULT_NBITERS ) ;

    fprintf( stderr, "  -l dolog  [ n        ]   log file or not { y n }\n" ) ;  /*V1.03-a*/

    fprintf( stderr, "  -m f p d  [%-4s %-2s %-3s]  model  {",
	     FamilyStrVC[ DEFAULT_FAMILY ],
	     ProporStrVC[ DEFAULT_PROPOR ],
	     DisperStrVC[ DEFAULT_DISPER ] ) ;
    {
      int imodel ;
      
      for ( imodel = 0 ; imodel < FAMILY_NB - 1 ; imodel ++ )
        fprintf( stderr, "%s ", FamilyStrVC[ imodel ] ) ;
      fprintf( stderr, "%s", FamilyStrVC[ FAMILY_NB - 1 ] ) ;
      fprintf( stderr, "} {" ) ;
      for ( imodel = 0 ; imodel < PROPOR_NB - 1 ; imodel ++ )
        fprintf( stderr, "%s ", ProporStrVC[ imodel ] ) ;
      fprintf( stderr, "%s", ProporStrVC[ PROPOR_NB - 1 ] ) ;
      fprintf( stderr, "} {" ) ;
      for ( imodel = 0 ; imodel < DISPER_NB - 1 ; imodel ++ )
        fprintf( stderr, "%s ", DisperStrVC[ imodel ] ) ;
      fprintf( stderr, "%s", DisperStrVC[ DISPER_NB - 1 ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -n neigh  [ %-8s ]   neighbour specification (image), { ",
            NeighStrVC[ DEFAULT_NEIGHSPEC ] ) ;
    {
      NeighET inei ;
      
      for ( inei = 0 ; inei < NEIGH_NB ; inei ++ )
        fprintf( stderr, "%s ", NeighStrVC[ inei ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -o fout   [ file     ]   output files basename\n" ) ;

    fprintf( stderr, "  -s init   [ s 1      ]   init : " ) ;
    fprintf( stderr, "s <v> = sort var <v>\n" ) ;
    fprintf( stderr, "%33s f <ini.uf> = from <ini.uf>\n", " " ) ; /*V1.04-a*/
    fprintf( stderr, "%33s r <n> = <n> random inits\n", " " ) ;
    fprintf( stderr, "%33s l <file> = use known labels from <file>\n", " " ) ;
    fprintf( stderr, "%33s mi/mf <para> = initial/fixed parameters\n", " " ) ;

    fprintf( stderr, "--- Press return for more options ---" ) ;
    getchar( ) ;

    fprintf( stderr, "  -t tie    [ %-8s ]   MAP tie rule  { ",
            TieStrVC[ DEFAULT_TIE ] ) ;
    {
      TieET itie ;
      
      for ( itie = 0 ; itie < TIE_NB ; itie ++ )
        fprintf( stderr, "%s ", TieStrVC[ itie ] ) ;
    }
    fprintf( stderr, "}\n" ) ;


    fprintf( stderr, "  -B bmod   [ %-8s ]   b estimation mode { ",
	     BtaStrVC[ DEFAULT_BTAMODE ] ) ;
    {
      BetaET bmod ;
      for (bmod = 0; bmod < BETA_NB; bmod ++)
	fprintf( stderr, "%s ", BtaStrVC[ bmod ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -C crit   [ %-8s ]   local maximum criterion { ",
	     CritStrVC[ DEFAULT_CRIT ] ) ;
    {
      CritET  crit ;
      for (crit = 0; crit < CRIT_NB; crit ++)
	fprintf( stderr, "%s ", CritStrVC[ crit ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -G nit conv step rand [ %2d %5.3f %3.1f %1d ]  %s\n",
	     DEFAULT_BTAGRADNIT, DEFAULT_BTAGRADCVTH, DEFAULT_BTAGRADSTEP,
	     DEFAULT_BTAGRADRAND,
	     "parameters of beta gradient estimation") ;

    fprintf( stderr, "  -H bstep bmax ddrop dloss lloss [ %3.2f %3.1f %3.1f %3.1f %3.2f ]\n %33s %s\n",
	     DEFAULT_BTAHEUSTEP, DEFAULT_BTAHEUMAX, DEFAULT_BTAHEUDDROP,
	     DEFAULT_BTAHEUDLOSS, DEFAULT_BTAHEULLOSS, " ",
	     "parameters of beta heuristic") ;

    fprintf( stderr, "  -I eiter  [ %-8d ]   number of E-step iterations (>= 1)\n",
             DEFAULT_NBEITERS ) ;/*V1.04-j*/


    fprintf( stderr, "  -M miss   [ %-8s ]   missing data processing { ",
            MissStrVC[ DEFAULT_MISSING ] ) ;
    {
      MissET imissing ;
      
      for ( imissing = 0 ; imissing < MISSING_NB ; imissing ++ )
        fprintf( stderr, "%s ", MissStrVC[ imissing ] ) ;
    }
    fprintf( stderr, "}\n" ) ;


    fprintf( stderr, "  -O order  [ %-8s ]   order of site visit { ",
            OrderStrVC[ DEFAULT_ORDER ] ) ;
    {
      OrderET iorder ;
      
      for ( iorder = 0 ; iorder < ORDER_NB ; iorder ++ )
        fprintf( stderr, "%s ", OrderStrVC[ iorder ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "  -S seed   [ <time>   ]   random generator seed \n" ) ;     /*V1.04-e*/

    fprintf( stderr, "  -T                       print debugging information\n" ) ;     /*V1.04-f*/

    fprintf( stderr, "  -U update [ %-8s ]   site update scheme  { ",
            UpdateStrVC[ DEFAULT_UPDATE ] ) ;
    {
      UpdET iupdate ;
      
      for ( iupdate = 0 ; iupdate < UPDATE_NB ; iupdate ++ )
        fprintf( stderr, "%s ", UpdateStrVC[ iupdate ] ) ;
    }
    fprintf( stderr, "}\n" ) ;

    fprintf( stderr, "\n\nYou may also just type arguments : \n" ) ;
    fprintf( stderr, "  -v                      versions information\n" ) ;
    fprintf( stderr, "  -h help_topic           longer help - help topics are\n {" ) ;
    {
      HelpET ihelp ;

      for ( ihelp = 0 ; ihelp < HELP_NB ; ihelp ++ )
	fprintf( stderr, " %s", HelpStrVC[ ihelp ] ) ;
    }
    fprintf( stderr, " } \n" ) ;

} /* end of PrintUsage() */


/* ------------------------------------------------------------------- */

static StatusET PrintHelp( const char* HelpType ) 

/* ------------------------------------------------------------------- */
{
  StatusET err = STS_I_DONE ;

  switch( GetEnum( HelpType, HelpStrVC, HELP_NB) )
    {
    case HELP_GENERAL:
      PrintHelpGeneral( stdout ) ;
      break ;

    case HELP_OPTIONS:
      PrintHelpOptions( stdout ) ;
      break ;

    case HELP_EXAMPLES:
      PrintHelpExamples( stdout ) ;
      break ;

    case HELP_FILEIN:
      PrintHelpFileIn( stdout ) ;
      break ;

    case HELP_FILEOUT:
      PrintHelpFileOut( stdout ) ;
      break ;

    case HELP_VERSIONS:
      PrintHelpVersions( stdout ) ;
      break ;

    default:
      fprintf( stderr, " Unknown help type %s\n", 
	       HelpType ) ;
      err = STS_E_ARG ;
    }  

  return err ;
}


/*V1.06-a*/
/* ------------------------------------------------------------------- */
static void my_strupr( char *s )
/* ------------------------------------------------------------------- */
{
    if ( s == NULL ) return ;
    for ( ; (*s) != '\0' ; s ++ )
    {
        if ( ( (*s) > 'a' ) && ( (*s) < 'z' ) )
        {
            (*s) = (*s) + ( 'A' - 'a' ) ;
        }
    }
}


/* ======================================================================= */
/*V1.06-c*/
/* start test of NemArgs 
int main( int argc, const char* argv[] )
{
  char        fname[ LEN_FILENAME + 1 ] ;
  char        datadescS[ LEN_LINE + 1 ] ;
  NemParaT    NemPara ;
  SpatialT    Spatial ;
  StatModelT  StatModel ;
  DataT       Data ;

  StatusET  err = 
    NemArgs( argc, argv, fname, & StatModel, & NemPara,
	     datadescS, & Data, & Spatial.NeighData.Image,
	     & Spatial.Type ) ;

  return err ;
}
 end test of NemArgs */
                                                                                                                                                                                                                                                                                 ./._nem_arg.h                                                                                       000755  000765  000765  00000000312 11541743325 012667  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_arg.h                                                                                           000755  000765  000024  00000001113 11541743325 012605  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         #include "nem_typ.h"    /* DataT, ... */

/*
    V1.06-a   17-JUN-1998  NoiseModel -> StatModel
 */

int NemArgs
    (
          int           Argc,
          const char    *Argv[],
          char          *Fname,         /* O */
          StatModelT    *StatModelP,    /* O */
          NemParaT      *NemParaP,      /* O */
	  char*         DatadescS,      /* O */
	  DataT*        DataP,          /* O */
	  ImageNeighT*  ImageP,         /* O */
	  TypeET*       SpatialTypeP    /* O */	  
    ) ;


void PrintUsage( const char* CmdName ) ;

extern const char *CritStrVC[ CRIT_NB ] ;
                                                                                                                                                                                                                                                                                                                                                                                                                                                     ./._nem_exe.c                                                                                       000755  000765  000765  00000000312 11541743327 012674  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_exe.c                                                                                           000755  000765  000024  00000140125 11541743327 012621  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\

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
                                                                                                                                                                                                                                                                                                                                                                                                                                           ./._nem_hlp.c                                                                                       000755  000765  000765  00000000312 11541743331 012671  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_hlp.c                                                                                           000755  000765  000024  00000111162 11541743331 012615  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\

    nem_hlp.c

    Project NEM : display of long help

    Thu Apr  8 19:24:47 1999

\*/

#include "nem_hlp.h"  /* Exported prototypes */

void PrintHelpGeneral( FILE* F )
{
    fprintf( F , "\n" ) ;
    fprintf( F , "Goal\n" ) ;
    fprintf( F , "====\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    This program computes a partition of a given set of objects \n" ) ;
    fprintf( F , "    described by one or several numeric variables and by their \n" ) ;
    fprintf( F , "    spatial relationships, using the 'Neighborhood EM' algorithm \n" ) ;
    fprintf( F , "    (NEM). This algorithm is derived from the EM algorithm applied\n" ) ;
    fprintf( F , "    to a hidden Markov random field model. Its new feature consists in\n" ) ;
    fprintf( F , "    taking into account some spatial interdependance between the \n" ) ;
    fprintf( F , "    objects.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    It may be used for:\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - unsupervised segmentation of color or gray-level images \n" ) ;
    fprintf( F , "      (points = pixel values, geographic position = pixel coordinates) ;\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - clustering of spatial data like socio-economical activities of\n" ) ;
    fprintf( F , "      neighbouring counties, etc.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    The algorithm takes as input an objects-variables table, \n" ) ;
    fprintf( F , "    and a specification of the neighborhood relationship between \n" ) ;
    fprintf( F , "    the objects. It produces as output a fuzzy or a hard partition \n" ) ;
    fprintf( F , "    of the objects.  The main algorithm is described in the following\n" ) ;
    fprintf( F , "    paper:\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Ambroise, C., Dang, V.M. and Govaert, G. (1997). Clustering of spatial\n" ) ;
    fprintf( F , "  data by the EM algorithm, in A.~Soares, J.~G\'omez-Hernandez and\n" ) ;
    fprintf( F , "  R.~Froidevaux, eds, `geoENV I - Geostatistics for Environmental\n" ) ;
    fprintf( F , "  Applications', Vol. 9 of `Quantitative Geology and Geostatistics', \n" ) ;
    fprintf( F , "  Kluwer Academic Publisher, pp.~493--504.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Changing default parameters\n" ) ;
    fprintf( F , "===========================\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    The behaviour of this clustering algorithm can be adjusted in many\n" ) ;
    fprintf( F , "    ways to fit a particular problem. The main possibilities\n" ) ;
    fprintf( F , "    are described below.  \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - The assumed degree of spatial interdependance is controlled by\n" ) ;
    fprintf( F , "    the value of the 'beta' coefficient (option '-b beta_value'). The higher\n" ) ;
    fprintf( F , "    it is, the smoother the partition will look in the geographic\n" ) ;
    fprintf( F , "    space, but the less it will fit to the data. The default value (1.0)\n" ) ;
    fprintf( F , "    seems to work well in most image segmentation problems where the\n" ) ;
    fprintf( F , "    patches are supposed to be spatially smooth. For beta's lowest\n" ) ;
    fprintf( F , "    value (0.0), the algorithm is the same as Dempster et al's EM\n" ) ;
    fprintf( F , "    algorithm (1974), and does a 'spatially blind' segmentation.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - The algorithm 'Neighborhood Classification EM' (NCEM) can be used\n" ) ;
    fprintf( F , "    instead of NEM (option '-a ncem'). The principle of NCEM consists in\n" ) ;
    fprintf( F , "    'hardening' the classification matrix at each iteration (C-step\n" ) ;
    fprintf( F , "    after the E-step). Practically, NCEM converges faster than NEM,\n" ) ;
    fprintf( F , "    but gives a poorer segmentation on data containing a high level of\n" ) ;
    fprintf( F , "    noise.  \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - In some applications, the class may be already known for a part of\n" ) ;
    fprintf( F , "    the sample. Such a knowledge can be taken into account by the\n" ) ;
    fprintf( F , "    Neighborhood EM algorithm, and may improve considerably the resulting\n" ) ;
    fprintf( F , "    classification (option '-s l file.ck').  \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Incomplete observations are taken into account in the\n" ) ;
    fprintf( F , "    probabilistic model and the program.  It is simply assumed that\n" ) ;
    fprintf( F , "    the missingness occurs at random, i.e. it does not depend on the\n" ) ;
    fprintf( F , "    missing value itself nor on the unobserved class.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Notice\n" ) ;
    fprintf( F , "======\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    The program is only provided to make it easier to test the\n" ) ;
    fprintf( F , "    behavior of the Neighborhood EM clustering algorithm and compare\n" ) ;
    fprintf( F , "    it with other algorithms.    Although I have tried to write and \n" ) ;
    fprintf( F , "    test the program as carefully as possible, it is not guaranteed \n" ) ;
    fprintf( F , "    to be error-free.  Please contact me if you have\n" ) ;
    fprintf( F , "    any questions or problems in using it.  Please also mention its\n" ) ;
    fprintf( F , "    origin if you use it for a published work. Finally I would be\n" ) ;
    fprintf( F , "    interested to know for what kind of problem you have found\n" ) ;
    fprintf( F , "    this program to be of use.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Van M� Dang\n" ) ;
    fprintf( F , "Van.Mo.Dang@utc.fr\n" ) ;
    fprintf( F , "http://www.hds.utc.fr/~mdang\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
} /* end of PrintHelpGeneral() */


void PrintHelpOptions( FILE* F )
{
    fprintf( F , "\n" ) ;
    fprintf( F , "Command Syntax\n" ) ;
    fprintf( F , "==============\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " Usage :    nem_exe   file  K  [ option1 option2 ... ]\n" ) ;
    fprintf( F , " ------- \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " Arguments :\n" ) ;
    fprintf( F , " -----------\n" ) ;
    fprintf( F , "   file       base name of input files ___.str and ___.dat\n" ) ;
    fprintf( F , "   K          number of classes\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " Options :   [ default  ] { possible values }\n" ) ;
    fprintf( F , " --------- \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " \n" ) ;
    fprintf( F , "  -a algo   [ nem      ]   { nem ncem gem }\n" ) ;
    fprintf( F , "     Algorithm to compute classification at E-step of each iteration : \n" ) ;
    fprintf( F , "      ncem = crisp classification by ICM procedure\n" ) ;
    fprintf( F , "      nem  = fuzzy classification by mean field approximation\n" ) ;
    fprintf( F , "      gem  = fuzzy classification by Gibbs sampling (Monte-Carlo simulations)\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -b beta   [ 1        ]   (0.0 <= beta <= 4)\n" ) ;
    fprintf( F , "     Coefficient of spatial smoothing to apply. This matches the \n" ) ;
    fprintf( F , "     Potts random field strength of interaction with 4-neighbor contexts.\n" ) ;
    fprintf( F , "     Notice that b = 0 is equivalent to EM for a mixture model.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -c wh thr [clas  0.04]   { none clas crit } and (>0)\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "       none     = no convergence test, i.e. do all specified iterations.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "       clas thr = stop the iterations when the largest difference\n" ) ;
    fprintf( F , "                  between previous and current classification matrix\n" ) ;
    fprintf( F , "                  is <= threshold. A threshold 0.04 is usually optimal.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "       crit thr = stop the iterations when \n" ) ;
    fprintf( F , "                  | (current_crit - last_crit)/current_crit | < threshold.\n" ) ;
    fprintf( F , "                  A threshold of 0.001 is usually best.  The test uses \n" ) ;
    fprintf( F , "                  the criterion selected by the option -C.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -f format [ hard     ]   { hard fuzzy }\n" ) ;
    fprintf( F , "     Format of output partition. Hard = N integers having values from\n" ) ;
    fprintf( F , "     1 to nk : for each observation, give the number of the class where\n" ) ;
    fprintf( F , "     it has highest grade of membership. Fuzzy = N x K reals between\n" ) ;
    fprintf( F , "     0 and 1 : for each observation, give its grade of membership in\n" ) ;
    fprintf( F , "     each class.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -i itmax  [ 100      ]   (>= 0)\n" ) ;
    fprintf( F , "     Maximum number of NEM iterations.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -l dolog  [ n        ]   { y n }\n" ) ;
    fprintf( F , "     Produce a log file or not to see the results of each iteration.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -m f p d  [norm p_ s__]  { norm lapl } { p_ pk } { s__ sk_ s_d skd }\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "     Mixture model assumption to use\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "      norm/lapl/bern : normal, Laplace or Bernoulli distributions\n" ) ;
    fprintf( F , "                       Bernoulli distributions are for binary data\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "      p_/pk     : clusters have equal / varying proportions\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "      s__/...   : variance model\n" ) ;
    fprintf( F , "                    s__ : same variance in all clusters and variables\n" ) ;
    fprintf( F , "                    sk_ : one variance per cluster, same in all variables\n" ) ;
    fprintf( F , "                    s_d : one variance per variable, same in all clusters\n" ) ;
    fprintf( F , "                    skd : one variance per cluster and variable\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "                  the variables are assumed independent within a cluster\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -n neigh  [ 4        ]   { 4 f }\n" ) ;
    fprintf( F , "     Neighborhood specification to use in the case of an image.\n" ) ;
    fprintf( F , "     4 = default 4-nearest neighbor system. f = specify neighborhood\n" ) ;
    fprintf( F , "     window in file.nei.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -o fout   [ file     ]\n" ) ;
    fprintf( F , "     Output files basename ___.cf, ___.mf and ___.log. Default is \n" ) ;
    fprintf( F , "     to use input file basename. Specify '-' to output the\n" ) ;
    fprintf( F , "     classification to standard output; useful to pipe the result to an \n" ) ;
    fprintf( F , "     'nem_exe -s f -' session.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -s init   [ s 1      ]   { s <v> | f <ini.uf> | r <n> | l <file> } \n" ) ;
    fprintf( F , "     Initialization mode.\n" ) ;
    fprintf( F , "     -s s <v>   \n" ) ;
    fprintf( F , "        Sort observations by variable <v>, then divide them\n" ) ;
    fprintf( F , "        in K quantiles of equal size to get initial partition.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "     -s f <ini.uf>\n" ) ;
    fprintf( F , "        Read initial fuzzy classification from file <ini.uf>. \n" ) ;
    fprintf( F , "        '-s f -' reads from standard input.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "     -s r <n>\n" ) ;
    fprintf( F , "        Start <n> times from random parameters (means chosen at\n" ) ;
    fprintf( F , "        random among the observations), then keep result with highest\n" ) ;
    fprintf( F , "        criterion.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "     -s l <file>\n" ) ;
    fprintf( F , "        Use partially known labels given in <file> to compute initial \n" ) ;
    fprintf( F , "        parameters. Those labels remain fixed throughout the\n" ) ;
    fprintf( F , "        clustering process.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "     -s mi <para> \n" ) ;
    fprintf( F , "     -s mf <para>\n" ) ;
    fprintf( F , "        Use specified parameters at beginning (mi) or throughout the\n" ) ;
    fprintf( F , "        clustering process (mf). Parameters syntax :\n" ) ;
    fprintf( F , "        p_1 ... p_{K-1}   m_11 .. m_1D m21 ... m_KD   s_11 ... s_KD.\n" ) ;
    fprintf( F , "        The s_kd are the standard errors for normal distributions, or\n" ) ;
    fprintf( F , "        the scale parameters for the Laplace distributions.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -t tie    [ random   ]   { random first }\n" ) ;
    fprintf( F , "     How to choose the class with highest probability when several\n" ) ;
    fprintf( F , "     classes have same maximum probability in MAP\n" ) ;
    fprintf( F , "     classification. 'random' draws uniformly between ex-aequo\n" ) ;
    fprintf( F , "     classes, 'first' chooses class with lowest index.\n" ) ;
    fprintf( F , " \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -B bmod   [ fix      ]   { fix psgrad heu_d heu_l }\n" ) ;
    fprintf( F , "     Procedure to estimate beta automatically :\n" ) ;
    fprintf( F , "      fix   = no estimation of beta, use beta given by option '-b'\n" ) ;
    fprintf( F , "      psgrad = pseudo-likelihood gradient ascent\n" ) ;
    fprintf( F , "      heu_d = heuristic using drop of fuzzy within cluster inertia\n" ) ;
    fprintf( F , "      heu_l = heuristic using drop of mixture likelihood\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -C crit   [ U        ]   local maximum criterion { U M D L }\n" ) ;
    fprintf( F , "     Criterion used to select the best local solution from random starts :\n" ) ;
    fprintf( F , "      U = fuzzy spatial clustering criterion\n" ) ;
    fprintf( F , "      M = fuzzy pseudo-likelihood\n" ) ;
    fprintf( F , "      D = fuzzy within cluster inertia\n" ) ;
    fprintf( F , "      L = likelihood of mixture parameters\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -G nit conv step rand [  1 0.001 0.0 0 ]  \n" ) ;
    fprintf( F , "     Parameters of beta gradient estimation\n" ) ;
    fprintf( F , "      nit = number of gradient iterations\n" ) ;
    fprintf( F , "      conv = threshold to test convergence (|g'|<conv*N is tested)\n" ) ;
    fprintf( F , "      step = > 0 for fixed step, <= 0 for Newton step = 1/g''\n" ) ;
    fprintf( F , "      rand = in -s r  init mode, initial beta random (1) or fixed by -b (0)\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -H bstep bmax ddrop dloss lloss [ 0.10 2.0 0.8 0.5 0.02 ]  \n" ) ;
    fprintf( F , "     Parameters of beta D and L heuristics :\n" ) ;
    fprintf( F , "      bstep = step of beta increase\n" ) ;
    fprintf( F , "      bmax  = maximal value of beta to test\n" ) ;
    fprintf( F , "      ddrop = threshold of allowed D drop (higher = less detection)\n" ) ;
    fprintf( F , "      dloss = threshold of allowed D loss (higher = less detection)\n" ) ;
    fprintf( F , "      lloss = threshold of allowed L loss (higher = less detection)\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -I eiter  [ 1        ]   (>= 1)\n" ) ;
    fprintf( F , "     Number of internal E-step iterations (nem and ncem algorithms),\n" ) ;
    fprintf( F , "     i.e. number of sweeps through whole dataset to compute\n" ) ;
    fprintf( F , "     classification at each iteration.  For gem algorithm, indicates\n" ) ;
    fprintf( F , "     number of sweeps through the dataset to compute the average\n" ) ;
    fprintf( F , "     frequency of class occurrence -> a large value is recommended\n" ) ;
    fprintf( F , "     for the gem algorithm (typically 50).\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -M miss   [ replace  ]    { replace ignore }\n" ) ;
    fprintf( F , "     How to deal with missing data.  Replace by expected value (EM) or\n" ) ;
    fprintf( F , "     ignore when computing mean (maximize fuzzy clustering criterion).\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -O order  [ direct   ]   order of site visit { direct random }\n" ) ;
    fprintf( F , "     Order in which to classify the observations at E-step.\n" ) ;
    fprintf( F , "      direct = given order 1..N\n" ) ;
    fprintf( F , "      random = random permutation of 1..N\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -S seed   [ <time>  ]    (integer)\n" ) ;
    fprintf( F , "     Specify a seed for the random number generator. Default uses\n" ) ;
    fprintf( F , "     current system clock.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -T test   [ n       ]    { y n }\n" ) ;
    fprintf( F , "     Print some debugging information or not.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -U update [ seq      ]   { seq para }\n" ) ;
    fprintf( F , "     Update the class of the sites in a sequential or parallel manner.\n" ) ;
    fprintf( F , "     'seq' works best. 'para' is more grounded, because it is EM with\n" ) ;
    fprintf( F , "     mean field approximation, but it requires a low spatial smoothing.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "You may also just type arguments : \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  -v                      versions information\n" ) ;
    fprintf( F , "  -h help_topic           longer help - help topics are\n" ) ;
    fprintf( F , "     general\n" ) ;
    fprintf( F , "     options\n" ) ;
    fprintf( F , "     examples\n" ) ;
    fprintf( F , "     filein\n" ) ;
    fprintf( F , "     fileout\n" ) ;
    fprintf( F , "     versions\n" ) ;
    fprintf( F , " \n" ) ;
    fprintf( F , " \n" ) ;
} /* end of PrintHelpOptions() */


void PrintHelpExamples( FILE* F )
{
    fprintf( F , "\n" ) ;
    fprintf( F , "    Examples :\n" ) ;
    fprintf( F , "    ----------\n" ) ;
    fprintf( F , " nem_exe  myimg  3 -b 0.5 -s r 10 -o Res/myimg3r_05 >&! Res/myimg3r_05.out &\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    This unix command clusters data set myimg (files myimg.dat, myimg.str)\n" ) ;
    fprintf( F , "    into 3 classes ; spatial coefficient is 0.5 ; do 10 random starts ;\n" ) ;
    fprintf( F , "    save results in files Res/myimg3r_05.* (___.cf ___.log ___.mf) .\n" ) ;
    fprintf( F , "    Last part of the command is unix-specific : it saves screen output\n" ) ;
    fprintf( F , "    to file Res/myimg3r_05.out and executes the program in background.\n" ) ;
    fprintf( F , "\n" ) ;
} /* end of PrintHelpExamples() */


void PrintHelpFileIn( FILE* F )
{
    fprintf( F , "\n" ) ;
    fprintf( F , "Input Files\n" ) ;
    fprintf( F , "===========\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    Input files are in ASCII format.\n" ) ;
    fprintf( F , "    2 or 3 input files are required : file.str   file.dat   [ file.nei ]\n" ) ;
    fprintf( F , "    2 optional input files :          file.u0    file.ck\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " 1) file.str\n" ) ;
    fprintf( F , " -----------\n" ) ;
    fprintf( F , "    This gives the structure of the data : \n" ) ;
    fprintf( F , "    type of spatial repartition (image, spatial or non-spatial),\n" ) ;
    fprintf( F , "    number of objects and variables. This file may start with\n" ) ;
    fprintf( F , "    comment lines to describe the dataset (lines beginning with #).\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Ex 1 : Color image (each pixel is described by 3 variables) \n" ) ;
    fprintf( F , "             of 200 lines (height) and 300 columns (width)\n" ) ;
    fprintf( F , "    # RGB biomedical coloscopic image. Look for 3 or 4 classes.\n" ) ;
    fprintf( F , "    I 200 300 3\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Ex 2 : Spatial dataset of 3000 objects described by 4 variables\n" ) ;
    fprintf( F , "    # Economic data on counties (region of Centre). About 5 classes.\n" ) ;
    fprintf( F , "    S 3000 4\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Ex 3 : Non-spatial dataset of 3000 objects described by 4 variables\n" ) ;
    fprintf( F , "    # Companies characteristics for risk assessment. About 6 classes.\n" ) ;
    fprintf( F , "    N 3000 4\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " 2) file.dat\n" ) ;
    fprintf( F , " -----------\n" ) ;
    fprintf( F , "    Contains the objects-variables table (only the non-spatial\n" ) ;
    fprintf( F , "    variables).  Missing data are specified by NaN.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    If the dataset is an image, the pixels must be listed \n" ) ;
    fprintf( F , "    line-by-line first, i.e. : \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    x(1,1) x(1,2) ... x(1,nc)  \n" ) ;
    fprintf( F , "    x(2,1) x(2,2) ... x(2,nc)\n" ) ;
    fprintf( F , "    ...\n" ) ;
    fprintf( F , "    x(nl,1) x(nl,2) ... x(nl,nc)\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    where nl = number of lines, nc = number of columns, and \n" ) ;
    fprintf( F , "    x(i,j) = values of pixel at line i and column j (a set of 3\n" ) ;
    fprintf( F , "    numbers for a color image).\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Ex 1 : Color image\n" ) ;
    fprintf( F , "    50 100 120\n" ) ;
    fprintf( F , "    51  99 122\n" ) ;
    fprintf( F , "    ...\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Ex 2 : Spatial data (4 variables)\n" ) ;
    fprintf( F , "    0.31 200 41 1200 \n" ) ;
    fprintf( F , "    0.28 202 43 1180\n" ) ;
    fprintf( F , "    ...\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " 3) file.nei\n" ) ;
    fprintf( F , " -----------\n" ) ;
    fprintf( F , "    Specifies the spatial relationships between the objects. \n" ) ;
    fprintf( F , "    * For an image, this file is optional and allows to specify\n" ) ;
    fprintf( F , "    a particular neighborhood system ; if no file is specified,\n" ) ;
    fprintf( F , "    default neighborhood system taken is 4-nearest neighbours.\n" ) ;
    fprintf( F , "    * For other spatial data, this file is required.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    Format is different for an image and other spatial data\n" ) ;
    fprintf( F , "    (see examples below).\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Ex 1 : 4 nearest-neighbours in an image (default)\n" ) ;
    fprintf( F , "    -1 1           /* at most 1 pixel on the left and 1 on the right */\n" ) ;
    fprintf( F , "    -1 1           /* at most 1 pixel up and 1 down */\n" ) ;
    fprintf( F , "    0 1 0          \n" ) ;
    fprintf( F , "    1 0 1          /* 4 equally weighted neighbors : */\n" ) ;
    fprintf( F , "    0 1 0	   /* up,left,right,down */            \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Ex 2 : other spatial data, unweighted neighborhood graph\n" ) ;
    fprintf( F , "    0              /* 0 = no weight specified (use default weight = 1.0) */\n" ) ;
    fprintf( F , "    1   3   2 5 7  /* object 1 has 3 neighbors : objects 2, 5 and 7 */\n" ) ;
    fprintf( F , "    2   2   1 3\n" ) ;
    fprintf( F , "    ...\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Ex 2 : other spatial data, weighted neighborhood graph\n" ) ;
    fprintf( F , "    1                           /* 1 = weights are specified */\n" ) ;
    fprintf( F , "    1   3   2 5 7  0.5 0.6 0.8  /* object 1 has 3 neighbors : 2, 5 and 7 */\n" ) ;
    fprintf( F , "    2   2   1 3    0.3 0.5\n" ) ;
    fprintf( F , "    ...\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " 4.a) file.ck  (option -s l file.ck)\n" ) ;
    fprintf( F , " ------------\n" ) ;
    fprintf( F , "    Gives the class of the points for which the label is already\n" ) ;
    fprintf( F , "    known. Expected format is N + 1 integers (N being the total\n" ) ;
    fprintf( F , "    number of objects in the sample) :\n" ) ;
    fprintf( F , "    - the first number is the number of classes K ;\n" ) ;
    fprintf( F , "    - the N following integers, in the same order as the\n" ) ;
    fprintf( F , "    objects in file.dat, indicate the class to which the\n" ) ;
    fprintf( F , "    corresponding object belongs.  If the value is 0 or greater than the\n" ) ;
    fprintf( F , "    number of classes, then the object is considered to have no known\n" ) ;
    fprintf( F , "    label (its label will be computed by the algorithm). \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    In the current implementation, each class must contain at least\n" ) ;
    fprintf( F , "    one observation with known label. Those pre-labeled observations\n" ) ;
    fprintf( F , "    are used to initialize the centers of the clusters.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    5     /* number of classes */\n" ) ;
    fprintf( F , "    0     /* object 1 has no known label */\n" ) ;
    fprintf( F , "    5     /* object 2 belongs to class 5 */\n" ) ;
    fprintf( F , "    1     /* object 3 belongs to class 1 */\n" ) ;
    fprintf( F , "    ...\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " 4.b) file.u0  (option -s f)\n" ) ;
    fprintf( F , " ------------\n" ) ;
    fprintf( F , "    Gives an initial fuzzy classification to start the algorithm. \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    0.9 0.1  /* object 1 : initial membership = 0.9/0.1 in class 1/2 */\n" ) ;
    fprintf( F , "    0.3 0.7  /* object 2 : initial membership = 0.3/0.7 in class 1/2 */\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
} /* end of PrintHelpFileIn() */


void PrintHelpFileOut( FILE* F )
{
    fprintf( F , "\n" ) ;
    fprintf( F , "Output Files\n" ) ;
    fprintf( F , "============\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    Output files are in ASCII format.\n" ) ;
    fprintf( F , "    2 files are output :              file.cf    file.mf\n" ) ;
    fprintf( F , "    1 optional file is output :       file.log\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " 1) file.cf (or file.uf for fuzzy partition)\n" ) ;
    fprintf( F , " -------------------------------------------\n" ) ;
    fprintf( F , "    Gives the partition found by the algorithm.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Ex 1 : Image segmented in 2 'hard' classes (linewise order, as the \n" ) ;
    fprintf( F , "           input data file)\n" ) ;
    fprintf( F , "    1        /* Object 1 belongs to class 1 */\n" ) ;
    fprintf( F , "    1        /* Object 2 belongs to class 1 */\n" ) ;
    fprintf( F , "    2        /* Object 3 belongs to class 2 */\n" ) ;
    fprintf( F , "    1\n" ) ;
    fprintf( F , "    2\n" ) ;
    fprintf( F , "    2\n" ) ;
    fprintf( F , "    ...\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "    - Ex 2 : Data segmented in 3 fuzzy classes\n" ) ;
    fprintf( F , "    0.2  0.7  0.1    /* Object 1 : class 1 with probability 0.2 , etc. */\n" ) ;
    fprintf( F , "    0.8  0.05 0.015\n" ) ;
    fprintf( F , "    ...\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " 2) file.mf\n" ) ;
    fprintf( F , " ----------\n" ) ;
    fprintf( F , "    \n" ) ;
    fprintf( F , "    This file gives the values of criteria optimized by the algorithm,\n" ) ;
    fprintf( F , "    and the parameters of the mixture calculated by the\n" ) ;
    fprintf( F , "    algorithm (means, proportions, scale parameter). Example: \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ START EXAMPLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" ) ;
    fprintf( F , "Command line : \n" ) ;
    fprintf( F , " \n" ) ;
    fprintf( F , "  nem_exe intro 2 -b 1.00 -f fuzzy -R intro.cr -o intro_b10 -s r 10 -C M \n" ) ;
    fprintf( F , " \n" ) ;
    fprintf( F , "Criteria U=NEM, D=Hathaway, L=mixture, M=markov ps-like, error\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "  56.3807    -2671.73    -1788.76    -2793.36   0.335938\n" ) ;
    fprintf( F , " \n" ) ;
    fprintf( F , "Beta (fixed)\n" ) ;
    fprintf( F , "  1.0000\n" ) ;
    fprintf( F , "Mu (4), Pk, and sigma (4) of the 2 classes\n" ) ;
    fprintf( F , " \n" ) ;
    fprintf( F , "   0.947  0.00456  0.0199  0.998   0.5    1.15226  1.15226  1.15226  1.15226 \n" ) ;
    fprintf( F , "   0.0614 1.03     0.864  -0.163   0.5    1.15226  1.15226  1.15226  1.15226 \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END EXAMPLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , " 3) file.log\n" ) ;
    fprintf( F , " -----------\n" ) ;
    fprintf( F , "    Details each iteration results (optimized criterion and \n" ) ;
    fprintf( F , "    class parameters).\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
} /* end of PrintHelpFileOut() */


void PrintHelpVersions( FILE* F )
{
    fprintf( F , "\n" ) ;
    fprintf( F , "History of modifications\n" ) ;
    fprintf( F , "========================\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Version 0.00  (01.02.1996)\n" ) ;
    fprintf( F , "------------\n" ) ;
    fprintf( F , "First version released on WEB. Initialization by sorting variable.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Version 1.00  (31.05.1996)  \n" ) ;
    fprintf( F , "------------\n" ) ;
    fprintf( F , "Random initializations. Image default neighborhood system. Long help.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Version 1.01  (27.06.1996)  \n" ) ;
    fprintf( F , "------------\n" ) ;
    fprintf( F , "Added an -a ncem option to implement the crisp clustering version of NEM.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Each observation is now updated sequentially in turn. The previous \n" ) ;
    fprintf( F , "parallel updating would produce chessboard-like images for high betas.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Version 1.02  (17.10.1996)\n" ) ;
    fprintf( F , "------------\n" ) ;
    fprintf( F , "Added the possibility to take into account a partial knowledge\n" ) ;
    fprintf( F , "of the classification into the clustering procedure \n" ) ;
    fprintf( F , "(option '-s l <file.ck>').\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Version 1.03  (02.10.1997)  \n" ) ;
    fprintf( F , "------------\n" ) ;
    fprintf( F , "If a partial knowledge of the classification is available, the\n" ) ;
    fprintf( F , "intial cluster centers are computed from the observations with\n" ) ;
    fprintf( F , "known labels.  \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "The log file is made optional (option '-l y').  The ___.mf file now\n" ) ;
    fprintf( F , "also contains the estimated cluster proportions, volumes and\n" ) ;
    fprintf( F , "covariance matrices and the values of the final criteria.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Version 1.04  (11.01.1998)  \n" ) ;
    fprintf( F , "------------\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Two heuristics are implemented to estimate the spatial smoothness\n" ) ;
    fprintf( F , "coefficient beta.  The first heuristic is based on detecting a sharp\n" ) ;
    fprintf( F , "drop in the fuzzy within-cluster inertia D, or a sufficient decrease\n" ) ;
    fprintf( F , "from its maximum value, when beta is slowly increased.  The second\n" ) ;
    fprintf( F , "heuristic is based on detecting a sufficient decrease of the\n" ) ;
    fprintf( F , "log-likelihood L of the mixture parameters from its maximum value.\n" ) ;
    fprintf( F , "The heuristics may be invoked with '-B heu_d' or '-B heu_l'. Their\n" ) ;
    fprintf( F , "default parameters may be changed with '-H ...'. \n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "The final partition may now be printed to standard ouput instead of to a\n" ) ;
    fprintf( F , "file (option '-o -').  The result can thus be redirected as an \n" ) ;
    fprintf( F , "initial partition to another nem_exe session's input.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "The initial partition may now be read from any file with option '-s f\n" ) ;
    fprintf( F , "foo.bar' (previously inputbasename.u0 was used).  In particular, the\n" ) ;
    fprintf( F , "partition may be read from standard input (option '-s f -').  This allows\n" ) ;
    fprintf( F , "to read the initial partition through a pipe from the result of a previous \n" ) ;
    fprintf( F , "nem_exe session.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "At each iteration, the fuzzy classification at the E-step may now be\n" ) ;
    fprintf( F , "computed by either of two methods :\n" ) ;
    fprintf( F , "- Neighborhood EM's fixed point technique (default, '-a nem')\n" ) ;
    fprintf( F , "- Gibbsian EM's Gibbs sampler technique ('-a gem').\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "A longer explanation is given for the options and successive versions\n" ) ;
    fprintf( F , "(option '-h helptopic').\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Other options have been added, to test the effect of alternative\n" ) ;
    fprintf( F , "parameter values.  Those options usually do not change considerably\n" ) ;
    fprintf( F , "the default behaviour of the algorithm :\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "- At the E-step of each iteration, the classification may now be\n" ) ;
    fprintf( F , "  updated in a random order instead of 1..N (option '-O random').\n" ) ;
    fprintf( F , "- The seed of the random number generator may now be given (option '-S seed').\n" ) ;
    fprintf( F , "- The number of E-step internal iterations may be changed (option '-I eiter').\n" ) ;
    fprintf( F , "- The convergence threshold may also be changed (option '-c cvthres').\n" ) ;
    fprintf( F , "- Another criterion than U may be chosen to select best result ('-C crit').\n" ) ;
    fprintf( F , "- Compute the classification error in two-class case ('-R refclass').\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "A few internal changes were also made to make the program portable to\n" ) ;
    fprintf( F , "the djgpp gcc compiler for MS-DOS (srand48 replaced by srandom).\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "Version 1.05  (09-APR-1998)  \n" ) ;
    fprintf( F , "------------\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "This version mainly adds the capability to deal with missing data.\n" ) ;
    fprintf( F , "This means that some of the N observation vectors may be incompletely\n" ) ;
    fprintf( F , "observed.  In the 'name.dat' file, use NaN to indicate unobserved\n" ) ;
    fprintf( F , "components of an observation vector.  Two slightly different\n" ) ;
    fprintf( F , "techniques are provided to deal with missing data.  The first and\n" ) ;
    fprintf( F , "default behaviour (invoked using switch '-M replace') roughly consists\n" ) ;
    fprintf( F , "in replacing any missing component with its expected value.  This\n" ) ;
    fprintf( F , "technique implements the EM procedure and finds parameters maximizing\n" ) ;
    fprintf( F , "the likelihood.  The alternative behaviour (invoked using switch '-M\n" ) ;
    fprintf( F , "ignore') consists in ignoring missing components.  This means that the\n" ) ;
    fprintf( F , "means are computed using only observed components.  This alternative\n" ) ;
    fprintf( F , "technique finds a classification matrix and parameters which maximize\n" ) ;
    fprintf( F , "the fuzzy classifying log-likelihood.  It appears to converge a bit\n" ) ;
    fprintf( F , "faster than the 'replace' mode.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "The iteration count is displayed in a more economic way now (all the\n" ) ;
    fprintf( F , "iteration numbers were displayed separately).\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "A new criterion is computed, the fuzzy pseudo-likelihood, named M.\n" ) ;
    fprintf( F , "Using this criterion in order to choose the best result \n" ) ;
    fprintf( F , "may prove less sensitive to the value of beta than using U.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "In the random start strategy, two initialization tactics have been\n" ) ;
    fprintf( F , "made more sensible.  The initial volumes are computed as whole volume\n" ) ;
    fprintf( F , "/ number of classes (the whole volume was used previously).  The means\n" ) ;
    fprintf( F , "are redrawn until all drawn means are different --- this avoids the\n" ) ;
    fprintf( F , "problem of artificially merging together two classes.\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "A few internal changes were also made in order to allow direct call to\n" ) ;
    fprintf( F , "the program from as a Matlab function.  This allowed to detect and\n" ) ;
    fprintf( F , "remove a few potential bugs that had gone unnoticed (a file not\n" ) ;
    fprintf( F , "closed, use of memory just after freeing it).\n" ) ;
    fprintf( F , "\n" ) ;
    fprintf( F , "\n" ) ;
} /* end of PrintHelpVersions() */


                                                                                                                                                                                                                                                                                                                                                                                                              ./._nem_hlp.h                                                                                       000755  000765  000765  00000000312 11541743331 012676  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_hlp.h                                                                                           000755  000765  000024  00000000517 11541743331 012623  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         #ifndef nem_hlp_H
#define nem_hlp_H

#include <stdio.h>  /* FILE */

extern void PrintHelpGeneral( FILE* F ) ;

extern void PrintHelpOptions( FILE* F ) ;

extern void PrintHelpExamples( FILE* F ) ;

extern void PrintHelpFileIn( FILE* F ) ;

extern void PrintHelpFileOut( FILE* F ) ;

extern void PrintHelpVersions( FILE* F ) ;

#endif
                                                                                                                                                                                 ./._nem_mod.c                                                                                       000755  000765  000765  00000000312 11541743333 012667  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_mod.c                                                                                           000755  000765  000024  00000233015 11541743333 012615  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\

    NEM_MOD.C

    Programme NEM (Neighborhood EM) : routines calcul parametres classes

    Van Mo DANG       Janvier 96


    Vers-mod  Date         Description

    1.03-a    01-NOV-1996  Add function SetIdMatrix + calls by ParaPk.. 
                           to output variance matrix in file don.mf
    1.04-a    08-JAN-1998  Fixed error in DensPkVkI log(vk) needs * Nd
    1.05-a    12-JAN-1998  MissMode in ParaP*V*I
    1.05-b    12-JAN-1998  Process missing data in DensPkVkI
    1.05-c    12-JAN-1998  Add ParaPkV*I_Missing() (empty)
    1.05-d    15-JAN-1998  Add CommonGaussDiagMissing() <= ParaPkV*I_Missing()
    1.05-e    26-JAN-1998  Call to GenAlloc
    1.06-a    29-JUN-1998  Reflect change to nem_mod.h
    1.06-b    29-JUN-1998  Reflect change to nem_typ.h
    1.06-c    29-JUN-1998  Compute Para.NbObs_KD
    1.06-d    02-JUL-1998  Bug: CommonGaus/replace, sum(x-m)^2 neq sum(x^2)-m^2
    1.06-e    02-JUL-1998  Bug: in VkI_Missing/replace, div inertia by nk*D
    1.06-f    02-JUL-1998  Add EstimParaLaplace and called subfunctions
    1.06-g    01-DEC-1998  FkP double* instead of *float in dens...
    1.07-a    26-FEB-1999  Add FAMILY_BERNOULLI in GetDensityFunc and EstimPara
    1.07-b    26-FEB-1999  Add DensBernoulli
    1.07-c    03-MAR-1999  Fix bug DensBernoulli: disp==0 may give nonzero dens
\*/

#include "genmemo.h"    /* GenAlloc */
#include "nem_typ.h"    /* DataT, ... */
#include "nem_mod.h"    /* ParaP_V_I, ... */
#include <stdio.h>      /* printf, ... */
#include <stdlib.h>     /* malloc, ... */
#include <string.h>     /* memcpy, ... */
#include <math.h>       /* exp, ... */
#include <values.h>     /* MAXFLOAT */


#define TWO_PI    2 * 3.14159 

/*V1.05-d*/
#define _IJ       ( ( i * D ) + j )    /* access a (.,D) matrix by (i,j) */
#define _HJ       ( ( h * D ) + j )    /* access a (.,D) matrix by (h,j) */
#define _IH       ( ( i * K ) + h )    /* access a (.,K) matrix by (i,h) */
#define sqr(x)    ((x)*(x))            /* macro for x^2 */



/* ==================== LOCAL FUNCTION PROTOTYPING =================== */


/* Indirect call by GetDensityFunc() */

static int DensNormalDiag      /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : model parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        ) ;


static int DensLaplaceDiag      /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : model parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        ) ;


static int DensBernoulli        /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : model parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        ) ;


/* Called by EstimPara() */

static StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimParaNormal 
(
  const float   *C_NK,      /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           Nk,         /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* O : estimated parameters */
) ;




static StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimParaLaplace
(
  const float   *C_NK,      /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           Nk,         /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* O : estimated parameters */
) ;




/* Called by EstimParaNormal() */

/*V1.05-d*/
static StatusET          /* Return status : OK, EMPTY or MEMORY */
CommonGaussDiag
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  int           N,             /* I : nb of observation vectors */
  int           D,             /* I : nb of variables */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  int           K,             /* I : number of clusters */
  MissET        Miss,          /* I : how to treat missing data */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        OldNewCen_KD,  /* I/O : old then updated means (K,D) */

  int*          EmptyK_P,      /* O : index of empty class (0 or 1 ..K) */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD,          /* O : size of a class and variable (K,D) */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */
 ) ; 


static void InerToDisp
(
 DisperET      DispType, /* I : dispersion model */
 int           N,        /* I : nb of observation vectors */
 int           Nk,       /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD,  /* O : dispersion in each class/variable */
 StatusET*     StsP      /* [O] : STS_E_FUNCARG if unknown DispType */
) ;



/* Called by InerToDisp() */

static void InerToDisp__
(
 int           N,        /* I : nb of observation vectors */
 int           Nk,       /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
) ;

static void InerToDispK_
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
) ;


static void InerToDisp_D
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
) ;


static void InerToDispKD
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
) ;



/* Called by EstimParaLaplace() */


/*V1.06-f*/
static StatusET          /* Return status : OK, EMPTY or MEMORY */
CommonLaplaceDiag
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  int           N,             /* I : nb of observation vectors */
  int           D,             /* I : nb of variables */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  int           K,             /* I : number of clusters */
  MissET        Miss,          /* I : how to treat missing data */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        OldNewCen_KD,  /* I/O : old then updated means (K,D) */

  int*          EmptyK_P,      /* O : index of empty class (0 or 1 ..K) */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD,          /* O : size of a class and variable (K,D) */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */
 ) ; 


#if 0 /* already declared above */

static void InerToDisp
(
 DisperET      DispType, /* I : dispersion model */
 int           N,        /* I : nb of observation vectors */
 int           Nk,       /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
) ;

#endif


/* Called by CommonLaplaceDiag() */

static void EstimSizes
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD           /* O : size of a class and variable (K,D) */  
 ) ;


static StatusET          /* Return status : OK, EMPTY */
EstimLaplaceCenters
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const MissET  Miss,          /* I : how to treat missing data */
  const float*  N_K,           /* I : size of a class (K) */
  const float*  N_KD,          /* I : size of a class and variable (K,D) */  
  const float*  OldCen_KD,     /* I : old centers */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        OldNewCen_KD,  /* I/O : eventually updated centers (K,D) */
  int*          EmptyK_P       /* O : index of empty class (0 or 1 ..K) */
 ) ;


static void EstimLaplaceIner
 ( 
  const float*  X_ND,          /* I : data matrix (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const MissET  Miss,          /* I : how to treat missing data */
  const float*  N_K,           /* I : size of a class (K) */
  const float*  N_KD,          /* I : size of a class and variable (K,D) */  
  const float*  OldCen_KD,     /* I : old centers */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */
  const float*  NewCen_KD,     /* I : new centers */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */  
 ) ;



/* Called by EstimLaplaceCenters() */

static void ComputeMedian
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const int     H,             /* I : current cluster */
  const int     J,             /* I : current variable */
  const float   totwei,        /* I : total weight */
  int*          ImedP,         /* O : index of median of observed values */
  float*        CumweiP,       /* O : cumulated weights until median obs. */
  float*        MedvalP        /* O : median value (eventually midway) */
 ) ;



static float FindMinInerLaplaceEM  /* ret: minimizer of expected inertium */
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const int     H,             /* I : current cluster */
  const int     J,             /* I : current variable */
  const int     Imed,          /* I : median position of observed data */
  const float   Cumwei,        /* I : cumulated weight until Xmed */
  const float   Cen0,          /* I : previous center */
  const float   Disp0,         /* I : previous dispersion */
  const float   Nhj,           /* I : number of observations */
  const float   Nmis           /* I : number of missing data */
 ) ;



/* Called by FindMinInerLaplaceEM() */

static float DerivInerDir      /* ret: +/- R'(Y), R expected inertia */ 
 (
  const float   Y,             /* I : point at which to compute R' */
  const float   Weidif,        /* I : "after+med-bef" or "bef+med-after" */
  const float   Intwei,        /* I : cumulated intermediate weight */
  const float   Nmis,          /* I : number of missing data */
  const float   Cen0,          /* I : previous center */
  const float   Disp0          /* I : previous dispersion */
 ) ;




/* ==================== GLOBAL FUNCTION DEFINITION =================== */


/* ------------------------------------------------------------------- */
int GetDensityFunc  /* STS_OK or STS_E_FUNCARG */
        (
            const ModelSpecT  *SpecP,           /* I */
            CompuDensFT**     CompuDensFP       /* O */
        )
/* ------------------------------------------------------------------- */
{

    switch( SpecP->ClassFamily )
    {
        case FAMILY_NORMAL:
	  *CompuDensFP = DensNormalDiag ;
	  return STS_OK ;

        case FAMILY_LAPLACE:
	  *CompuDensFP = DensLaplaceDiag ;
	  return STS_OK ;

        case FAMILY_BERNOULLI:  /*V1.07-a*/
	  *CompuDensFP = DensBernoulli ;
	  return STS_OK ;

        default :
	  *CompuDensFP = NULL ;
	  fprintf( stderr, "GetDensityFunc bad arg : family = %d\n",
                   SpecP->ClassFamily ) ;
	  return STS_E_FUNCARG ;
    }

  /*???*/
}   /* end of GetDensityFunc() */




/* ------------------------------------------------------------------- */
StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimPara 
(
  const float   *C_NK,      /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           Nk,         /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* I/O : previous and new estimated parameters */
)
/* ------------------------------------------------------------------- */
{
  StatusET    sts ;  /* return status */

  int         k ;

  /* Family dependent estimation method 
   */
  switch( SpecP->ClassFamily ) {
  case FAMILY_NORMAL:
    sts = EstimParaNormal( C_NK, DataP, Nk, MissMode, SpecP, 
			   EmptyK_P, ParaP ) ;
    break ;

  case FAMILY_LAPLACE:
    sts = EstimParaLaplace( C_NK, DataP, Nk, MissMode, SpecP, 
			    EmptyK_P, ParaP ) ;      
    break;

  case FAMILY_BERNOULLI: /*+++DANGER: MissMode forced to MISSING_IGNORE+++*/
    sts = EstimParaLaplace( C_NK, DataP, Nk, MISSING_IGNORE, SpecP, 
			    EmptyK_P, ParaP ) ;  /*V1.07-a*/
    break;

  default:
    sts = STS_E_FUNCARG ;
  }
  
  /* Estimate proportions */
  if ( SpecP->ClassPropor == PROPOR_K ) {

    for ( k = 0; k < Nk ; k ++ )
      ParaP->Prop_K[ k ] = ParaP->NbObs_K[ k ] / DataP->NbPts ;
  }
  else {
    
    for ( k = 0; k < Nk ; k ++ )
      ParaP->Prop_K[ k ] = 1.0 / Nk ;
  }

  return sts ;
  /*???*/
}   /* end of EstimPara() */




/* ==================== LOCAL FUNCTION DEFINITION =================== */





/* ------------------------------------------------------------------- */
static int DensNormalDiag      /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : noise parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        )
/* ------------------------------------------------------------------- */
{
  int     d ;      /* 0..Nd-1 : current variable */
  float   dk ;     /* sum_d [ log(2pi skd2) + (xd-mkd)^2/skd2 ] */
  int     nbobs ;  /* 0..Nd-1 : nb of observed variables */ /*V1.05-b*/
  int     nuldisp; /* TRUE if one of the dispersions is zero */

  /* Gaussian mixture density, diagonal models :

     Complete data :
     fk(x) = [ prod_d (2*pi * vkd)^(-0.5) ] * exp( -0.5 * dist_x_mk )
     where  dist_x_mk = sum[d=1,Nd]( (x(d) - mk(d))^2 ) / vkd

     Incomplete data :
     fk(x) = [ prod_{d in oi} (2*pi * vkd)^(-0.5) ] * exp( -0.5 * dist_x_mk )
     where   dist_x_mk = sum[d in oi]  (x(d) - mk(d))^2 / vk
     =>
     log fk(x) = - 0.5 sum_{d in oi}[ log(2*pi * vkd) + (xd-mkd)^2/vkd ]

  */

  for ( d = 0, 
          dk = 0.0, nbobs = 0, nuldisp = 0 ; 
        d < Nd ; 
        d ++ )
    {
      if ( ! isnan( XV[ d ] ) )   /*V1.05-b*/
        {
          float dif = XV[ d ] - ParaP->Center_KD[ (Ik * Nd) + d ] ;

          if ( ParaP->Disp_KD[ (Ik * Nd) + d ] > EPSILON )
            dk = dk 
	      + log( TWO_PI * ParaP->Disp_KD[ (Ik * Nd) + d ] ) 
	      + (dif * dif) / ParaP->Disp_KD[ (Ik * Nd) + d ] ;
          else
            nuldisp = 1 ;

          nbobs ++ ;
        }
  }

  if ( ! nuldisp )
    {
      /*V1.04-a*//*V1.05-b*/
      *LogFkP = -0.5 * dk ;
      *FkP = exp( *LogFkP ) ;
      return 0 ;
    }
  else
    {
      *LogFkP = - MAXFLOAT ;
      *FkP = 0.0 ;
      return -1 ;
    }


}   /* end of DensNormalDiag() */



/* ------------------------------------------------------------------- */
static int DensLaplaceDiag      /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : noise parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        )
/* ------------------------------------------------------------------- */
{
  int     d ;      /* 0..Nd-1 : current variable */
  float   dk ;     /* sum_d [ log(2 lkd) + |xd-mkd| / lkd ] */
  int     nbobs ;  /* 0..Nd-1 : nb of observed variables */ /*V1.05-b*/
  int     nuldisp; /* TRUE if one of the dispersions is zero */

  /* Laplace mixture density, diagonal models :

     Complete data :
     fk(x) = [ prod_d (2*vkd)^(-1) ] * exp( - dist_x_mk )
     where  dist_x_mk = sum[d=1,Nd] | x(d) - mk(d) | / vkd

     Incomplete data :
     fk(x) = [ prod_{d in oi} (2*vkd)^(-1) ] * exp( - dist_x_mk )
     where   dist_x_mk = sum[d in oi]  | x(d) - mk(d) | / vkd
     =>
     log fk(x) = - sum_{d in oi}[ log(2*vkd) + |xd-mkd| / vkd ]

  */

  for ( d = 0, 
          dk = 0.0, nbobs = 0, nuldisp = 0 ; 
        d < Nd ; 
        d ++ )
    {
      if ( ! isnan( XV[ d ] ) )   /*V1.05-b*/
        {
          float dif = XV[ d ] - ParaP->Center_KD[ (Ik * Nd) + d ] ;

          if ( ParaP->Disp_KD[ (Ik * Nd) + d ] > EPSILON )
            dk = dk 
      	+ log( 2 * ParaP->Disp_KD[ (Ik * Nd) + d ] ) 
      	+ fabs( dif ) / ParaP->Disp_KD[ (Ik * Nd) + d ] ;
          else
            nuldisp = 1 ;

          nbobs ++ ;
        }
  }

  if ( ! nuldisp )
    {
      /*V1.04-a*//*V1.05-b*/
      *LogFkP = - dk ;
      *FkP = exp( *LogFkP ) ;
      return 0 ;
    }
  else
    {
      *LogFkP = - MAXFLOAT ;
      *FkP = 0.0 ;
      return -1 ;
    }

}   /* end of DensLaplaceDiag() */


/* ------------------------------------------------------------------- */
static int DensBernoulli      /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,         /* I : point dimension */
            int                Ik,         /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP ,     /* I : noise parameters */
            const float*       XV,         /* I : point (dim d) */
            double*            FkP,        /* O : density for class Ik */
            float*             LogFkP      /* O : log of density */
        )
/* ------------------------------------------------------------------- */
{
  int     d ;      /* 0..Nd-1 : current variable */
  float   dk ;     /* sum_d [ -log(1-vkd) + log{(1-vkd)/vkd} |xd-mkd| ] */
  int     nbobs ;  /* 0..Nd-1 : nb of observed variables */ /*V1.05-b*/
  int     nuldens; /* TRUE if the probability is zero */

  /* Bernoulli density, diagonal models :

     Complete data :
     fk(x) = prod_d v_kd^|x_d - m_kd| (1-v_kd)^( 1 - |x_d - m_kd| )

     Incomplete data :
     fk(x) = prod_{d in oi} v_kd^|x_d - m_kd| (1-v_kd)^( 1 - |x_d - m_kd| )
     =>
     log fk(x) = - sum_{d in oi} [ -log(1-vkd) + |xd-mkd| * log((1-vkd)/vkd) ]

     if there is a d such that vkd = 0 and |xd-mkd| != 0, prob=0.

  */

  for ( d = 0, 
          dk = 0.0, nbobs = 0, nuldens = 0 ; 
        d < Nd ; 
        d ++ ) {

      if ( ! isnan( XV[ d ] ) )   /*V1.05-b*/
        {
	  float disp = ParaP->Disp_KD[ (Ik * Nd) + d ] ;
          int   absdif = 
	    abs( (int) ( XV[ d ] - ParaP->Center_KD[ (Ik * Nd) + d ] ) ) ;

          if ( disp > EPSILON )
            dk = dk + absdif * log( ( 1 - disp ) / disp ) - log( 1 - disp ) ;
          else  /* null dispersion */ {

	    if ( absdif != 0 ) /* prob(xid != center) = 0 */ {
	      nuldens = 1 ;
	    }
	    else /* prob(xid = center) = 1 => no change to dk */ {
	      dk = dk + 0.0 - 0.0 ;
	    }
	  }

          nbobs ++ ;
        }
  }

  if ( ! nuldens )
    {
      /*V1.04-a*//*V1.05-b*/
      *LogFkP = - dk ;
      *FkP = exp( *LogFkP ) ;
      return 0 ;
    }
  else
    {
      *LogFkP = - MAXFLOAT ;
      *FkP = 0.0 ;
      return -1 ;
    }

}   /* end of DensBernoulli() */


/* ------------------------------------------------------------------- */
static StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimParaNormal 
(
  const float   *C_NK,      /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           Nk,         /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* O : estimated parameters */
)
/* ------------------------------------------------------------------- */
{
  int      N = DataP->NbPts ;
  int      D = DataP->NbVars ;

  StatusET sts ;


  /* Common computation for any diagonal gaussian model with missing data */
  sts = CommonGaussDiag
    ( DataP->PointsM, N, D, C_NK, Nk, MissMode, ParaP->Disp_KD, 
      ParaP->Center_KD, 
      EmptyK_P, ParaP->NbObs_K, ParaP->NbObs_KD, ParaP->Iner_KD ) ;

  /* Dispersion is computed depending on dispersion model */
  InerToDisp( SpecP->ClassDisper, N, Nk, D, ParaP->NbObs_K, ParaP->NbObs_KD, 
	      ParaP->Iner_KD, MissMode, 
	      ParaP->Disp_KD, & sts ) ;

  return sts ;

}   /* end of EstimParaNormal() */



/* ------------------------------------------------------------------- */
static StatusET          /* Return status : OK, EMPTY or MEMORY */
CommonGaussDiag
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  int           N,             /* I : nb of observation vectors */
  int           D,             /* I : nb of variables */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  int           K,             /* I : number of clusters */
  MissET        Miss,          /* I : how to treat missing data */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        OldNewCen_KD,  /* I/O : old then updated means (K,D) */

  int*          EmptyK_P,      /* O : index of empty class (0 or 1 ..K) */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD,          /* O : size of a class and variable (K,D) */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */
 ) 
/*\

   This function does the common computation that is needed at the
   M-step, in all the diagonal Gaussian models with missing data.

   The routine takes as input a data matrix Xij_nd (possibly
   containing NaN values), a classification matrix Cih_nk, and the
   kind of missing data processing Miss (REPLACE or IGNORE) ; and also
   the old volumes OldVh_k, covariance matrices OldChjj_kdd, and means
   OldNewMhj_kd of the classes, to be used in the REPLACE mode.

   It produces as output the reestimated mean vectors in OldNewMhj_kd;
   reestimated proportions in NewP_k; the last encountered empty class
   in EmptyK_P (or 0 if no empty class); the {sum[i] cih} in Nh_k, the
   {sum[i] cih rij} in Nhj_kd (number of present observations in class h
   and variable j); and the "missing data dispersion" of each
   class/variable in Disphj_kd (i.e. the missing data equivalent of
   {sum[i] Cih (xij - mhj)^2}).

\*/
/* ------------------------------------------------------------------- */
{
  const char* func = "CommonGaussDiag" ;

  float*   sumdata_KD ;    /* sum[i] cih rij xij   (K,D), local alloc */
  float*   sumsquare_KD ;  /* sum[i] cih rij xij^2 (K,D), local alloc */
  float*   oldmean_KD ;    /* saved previous means (K,D), local alloc */

  int      h ;             /* current class    0..K-1 */
  int      j ;             /* current variable 0..D-1 */
  int      i ;             /* current object   0..N-1 */
  StatusET sts ;           /* return status */

  float    obsinerhj ;


  sts = STS_OK ;

  /* Allocate local structures */
  sumdata_KD   = GenAlloc( K * D, sizeof(float), 1, func, "sumdata" ) ;
  sumsquare_KD = GenAlloc( K * D, sizeof(float), 1, func, "sumsquare" );
  oldmean_KD   = GenAlloc( K * D, sizeof(float), 1, func, "oldmean" ) ;

  /* Save old mean value */
  memcpy( oldmean_KD , OldNewCen_KD , K * D * sizeof( float ) ) ;

  /* Set no empty class by default */
  (*EmptyK_P) = 0 ;

  /* For each class h */
  for ( h = 0 ; h < K ; h ++ ) {

      /* For each variable j */
      for ( j = 0 ; j < D ; j ++ ) {

	  /* Initialize the sum[i] quantities to 0 */
	  N_K[ h ]             = 0.0 ;
	  N_KD[ _HJ ]          = 0.0 ;
	  sumdata_KD[ _HJ ]    = 0.0 ;
	  sumsquare_KD[ _HJ ]  = 0.0 ;

	  /* For each object i */
	  for ( i = 0 ; i < N ; i ++ ) {

	      float cih = C_NK[ _IH ] ;
	      float xij = X_ND[ _IJ ] ;

	      /* Increment class size */
	      N_K[ h ] += cih ;

	      /* If this object's j variable is observed */
	      if ( ! isnan( xij ) )
		{
		  /* Increment all sum[i] of observed quantities */
		  N_KD[ _HJ ]          += cih ;
		  sumdata_KD[ _HJ ]    += cih * xij ;
		  sumsquare_KD[ _HJ ]  += cih * xij * xij ;
		}
	      /* Else xij is missing, don't change sum[i] of observed data */
	  }
	  /* End  For each object i (all sums[i] are now done) */

	  /* If this class size > 0 */
	  if ( N_K[ h ] > 0 ) {
	      /* 
	       * Compute means Mhj and dispersions Dhj from computed sums,
	       * depending on REPLACE or IGNORE missing mode 
	       */
	      if ( Miss == MISSING_REPLACE ) {
		  /* 
		   * new_mhj = 
		   *  ( sum[i] cih rij xij + ( nh - nhj ) old_mhj ) / nh 
		   *
		   * disp_hj = 
		   *  "observed dispersion"_hj +
		   *  ( nh - nhj ) [ (new_mhj - old_mhj)^2 + old_vh old_chjj ]
		   */

		  OldNewCen_KD[ _HJ ] = 
		    ( sumdata_KD[ _HJ ] + 
		      ( N_K[ h ] - N_KD[ _HJ ] ) * oldmean_KD[ _HJ ] )
		    / N_K[ h ] ;

		  /* V1.06-d
		   * "observed dispersion"_hj = 
		   *      sum[i] cih rij (xij - new_mhj)^2 =
		   *      sum[i] cih rij xij^2 - 
		   *       new_mhj * ( 2 * sum[i] cih rij xij - nhj new_mhj )
		   */
		  obsinerhj = sumsquare_KD[ _HJ ] - 
		    OldNewCen_KD[ _HJ ] * 
		    (2*sumdata_KD[ _HJ ] - N_KD[ _HJ ] * OldNewCen_KD[ _HJ ]) ;

		  Iner_KD[ _HJ ] = 
		    obsinerhj +
		    ( N_K[ h ] - N_KD[ _HJ ] ) *
		    ( sqr( OldNewCen_KD[ _HJ ] - oldmean_KD[ _HJ ] ) +
		      OldDisp_KD[ _HJ ] ) ;
	      }
	      else  { /* assume Miss == MISSING_IGNORE */

		  /* 
		   * new_mhj = 
		   *  . if  nhj <> 0,  ( sum[i] cih rij xij ) / nhj  
		   *  . if  nhj == 0,  old_mhj
		   *
		   * disp_hj = 
		   *  "observed dispersion"_hj 
		   */

		  if ( N_KD[ _HJ ] > 0 )
		    OldNewCen_KD[ _HJ ] = sumdata_KD[ _HJ ] / N_KD[ _HJ ];
		  else
		    OldNewCen_KD[ _HJ ] = oldmean_KD[ _HJ ] ;

		  /*
		   * "observed dispersion"_hj = 
		   *      sum[i] cih rij (xij - new_mhj)^2 =
		   *      sum[i] cih rij xij^2 - nhj new_mhj^2
		   */
		  obsinerhj = sumsquare_KD[ _HJ ] - 
		    N_KD[ _HJ ] * sqr( OldNewCen_KD[ _HJ ] ) ;

		  Iner_KD[ _HJ ] = obsinerhj ;
	      }

	  }
	  else /* then this class size == 0 */ {
	      /* signal empty class and which one */
	      sts = STS_W_EMPTYCLASS ;
	      (*EmptyK_P) = h + 1 ;
	  }

      }
      /* end  For each variable j */

  }
  /* end  For each class h */


  /* Free local structures */
  GenFree( oldmean_KD ) ;
  GenFree( sumsquare_KD ) ;
  GenFree( sumdata_KD  ) ;

  return sts ;

}   /* end of CommonGaussDiag() */



/* ------------------------------------------------------------------- */
static void InerToDisp
(
 DisperET      DispType, /* I : dispersion model */
 int           N,        /* I : nb of observation vectors */
 int           Nk,       /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD,  /* O : dispersion in each class/variable */
 StatusET*     StsP      /* [O] : STS_E_FUNCARG if unknown DispType */
)
/* ------------------------------------------------------------------- */
{
  switch( DispType ) {
  case DISPER___:  /* Common dispersion for all classes and variables */
    InerToDisp__( N, Nk, D, NbObs_K, NbObs_KD, Iner_KD, MissMode, 
		  Disp_KD ) ;
    break ;

  case DISPER_K_:  /* In each class, common dispersion for all variables */
    InerToDispK_( N, Nk, D, NbObs_K, NbObs_KD, Iner_KD, MissMode, 
		  Disp_KD ) ;
    break ;

  case DISPER__D:  /* In each class, common dispersion for all variables */
    InerToDisp_D( N, Nk, D, NbObs_K, NbObs_KD, Iner_KD, MissMode, 
		  Disp_KD ) ;
    break ;

  case DISPER_KD:  /* In each class, common dispersion for all variables */
    InerToDispKD( N, Nk, D, NbObs_K, NbObs_KD, Iner_KD, MissMode, 
		  Disp_KD ) ;
    break ;

  default:
    (*StsP) = STS_E_FUNCARG ;
  }
}   /* end of InerToDisp() */


/* ------------------------------------------------------------------- */
static void InerToDisp__
(
 int           N,        /* I : nb of observation vectors */
 int           Nk,       /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
)
/* ------------------------------------------------------------------- */
{
  int      k ;         /* current class    0..K-1 */
  int      d ;         /* current variable 0..D-1 */

  float    vol ;       /* common volume */
  float    nobs ;      /* number of effective observations : ]0,N*D] */

  /* Compute common volume : 
     sum dispersions over the classes and variables 
     then divide by all/observed data size */
  
  vol = 0.0 ;
  nobs = 0.0 ;

  for ( k = 0; k < Nk ; k ++ )
    {
      /* If not empty class */
      if ( NbObs_K[ k ] > 0 )
	{
	  for ( d = 0 ; d < D ; d ++ )
	    {
	      vol  += Iner_KD[ k * D + d ] ;
	      nobs += NbObs_KD[ k * D + d ] ; 
	    }
	}
    }

  if ( MissMode == MISSING_REPLACE )
    vol /= ( N * D ) ;
  else
    vol /= nobs ;  /* nobs must be = N * D - DataP->NbMiss */

  /* For each class k and variable d */
  for ( k = 0 ; k < Nk ; k ++ )
    for ( d = 0 ; d < D ; d ++ )
    {
      /* Set its volume to the common volume */
      Disp_KD[ k * D + d ] = vol ;
    }

}   /* end of InerToDisp__() */


/* ------------------------------------------------------------------- */
static void InerToDispK_
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
)
/* ------------------------------------------------------------------- */
{
  int      k ;         /* current class    0..K-1 */
  int      d ;         /* current variable 0..D-1 */

  float    sumd_nkd ;    /* size of obs class */
  float    sumd_inerkd ; /* obs inertia in class k */
  float    dispk ;       /* dispersion in class k */


  /* For each class k */
  for ( k = 0; k < K ; k ++ ) {

      /* If not empty class */
      if ( NbObs_K[ k ] > 0 ) {

	  /* Reestimate its volume: sum inertia over the variables 
	     then divide by all/observed data size in this class */

	  sumd_nkd    = 0.0 ;
	  sumd_inerkd = 0.0 ;  /* use first component as sum first */

	  for ( d = 0 ; d < D ; d ++ )
	    {
	      sumd_nkd     += NbObs_KD[ k * D + d ] ;
	      sumd_inerkd  += Iner_KD[ k * D + d ] ;
	    }

	  if ( MissMode == MISSING_REPLACE )
	    dispk = sumd_inerkd /= ( D * NbObs_K[ k ] ) ;  /*V1.06-e*/
	  else
	    dispk = sumd_inerkd / sumd_nkd ;

	  /* Copy same dispersion into all variables of this class */
	  for ( d = 0 ; d < D ; d ++ )
	    {
	      Disp_KD[ k * D + d ] = dispk ;
	    }
      }
      /* (Else class k is empty : keep old dispersion */ 

    }

  /*??*/

}   /* end of InerToDispK_() */




/* ------------------------------------------------------------------- */
static void InerToDisp_D
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
)
/* ------------------------------------------------------------------- */
{
  int      k ;         /* current class    0..K-1 */
  int      d ;         /* current variable 0..D-1 */

  float    sumk_nkd ;    /* obs size of variable d */
  float    sumk_inerkd ; /* obs inertia in variable d */
  float    dispd ;       /* dispersion of variable d */

  /* For each variable d */
  for ( d = 0; d < D ; d ++ ) {

    /* Reestimate its volume: sum inertia over the classes 
       then divide by all/observed data size in this variable */

    sumk_nkd    = 0.0 ;
    sumk_inerkd = 0.0 ;

    for ( k = 0 ; k < K ; k ++ ) {
      sumk_nkd    += NbObs_KD[ k * D + d ] ;
      sumk_inerkd += Iner_KD[ k * D + d ] ;
    }

    if ( MissMode == MISSING_REPLACE )
      dispd = sumk_inerkd / N ;
    else
      dispd = sumk_inerkd / sumk_nkd ;

    /* Copy same dispersion of this variable into all classes */
    for ( k = 0 ; k < K ; k ++ ) {
      Disp_KD[ k * D + d ] = dispd ;
    }
  }

  /*???*/

}   /* end of InerToDisp_D() */



/* ------------------------------------------------------------------- */
static void InerToDispKD
(
 int           N,        /* I : nb of observation vectors */
 int           K,        /* I : nb of classes */
 int           D,        /* I : nb of variables */
 const float*  NbObs_K,  /* I : nb of observation vectors in each class */
 const float*  NbObs_KD, /* I : nb of observations in each class/variable*/
 const float*  Iner_KD,  /* I : inertia of each class/variable */
 MissET        MissMode, /* I : how to treat missing data */ /*V1.05-b*/
 float*        Disp_KD   /* O : dispersion in each class/variable */
)
/* ------------------------------------------------------------------- */
{
  int      k ;         /* current class    0..K-1 */
  int      d ;         /* current variable 0..D-1 */

  /* For each class k */
  for ( k = 0 ; k < K ; k ++ ) {

    /* For each variable d */
    for ( d = 0; d < D ; d ++ ) {

      /* Reestimate its volume: divide inertia by all/observed data 
	 size in this variable/ class if not empty */

      if ( MissMode == MISSING_REPLACE ) {
	if ( NbObs_K[ k ] > EPSILON )
	  Disp_KD[ k * D + d ] = Iner_KD[ k * D + d ] / NbObs_K[ k ] ;
	/* else empty class : do not update dispersion */
      }
      else /* IGNORE */ {
	if ( NbObs_KD[ k * D + d ] > EPSILON )
	  Disp_KD[ k * D + d ] = Iner_KD[ k * D + d ] / NbObs_KD[ k * D + d ] ;
      }
    }
  }

  /*???*/

}   /* end of InerToDispKD() */




/* ------------------------------------------------------------------- */
static StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY */ /*V1.06-b*/
EstimParaLaplace
(
  const float   *C_NK,      /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           Nk,         /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/
  const ModelSpecT* SpecP,  /* I : model specification */ /*V1.06-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 */ /*V1.05-b*/
  ModelParaT    *ParaP      /* O : estimated parameters */
)
/* ------------------------------------------------------------------- */
{
  int      N = DataP->NbPts ;
  int      D = DataP->NbVars ;

  StatusET sts ;


  /* Common computation for any diagonal gaussian model with missing data */
  sts = CommonLaplaceDiag
    ( DataP->PointsM, DataP->SortPos_ND, N, D, C_NK, Nk, MissMode, 
      ParaP->Disp_KD, 
      ParaP->Center_KD, 
      EmptyK_P, ParaP->NbObs_K, ParaP->NbObs_KD, ParaP->Iner_KD ) ;

  /* Dispersion is computed depending on dispersion model */
  InerToDisp( SpecP->ClassDisper, N, Nk, D, ParaP->NbObs_K, ParaP->NbObs_KD, 
	      ParaP->Iner_KD, MissMode, 
	      ParaP->Disp_KD, & sts ) ;

  return sts ;

  /*???*/
}   /* end of EstimParaLaplace() */


/*V1.06-f*/
/* ------------------------------------------------------------------- */
static StatusET          /* Return status : OK, EMPTY or MEMORY */
CommonLaplaceDiag
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  int           N,             /* I : nb of observation vectors */
  int           D,             /* I : nb of variables */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  int           K,             /* I : number of clusters */
  MissET        Miss,          /* I : how to treat missing data */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        OldNewCen_KD,  /* I/O : old then updated means (K,D) */

  int*          EmptyK_P,      /* O : index of empty class (0 or 1 ..K) */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD,          /* O : size of a class and variable (K,D) */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */
 )
/* ------------------------------------------------------------------- */
{
  const char* func = "CommonLaplaceDiag" ;

  float*   oldcent_KD ;    /* saved previous centers (K,D), local alloc */

  StatusET sts ;           /* return status */


  /* Allocate local structures */
  oldcent_KD   = GenAlloc( K * D, sizeof(float), 1, func, "oldcent_KD" ) ;

  /* Save old mean value */
  memcpy( oldcent_KD , OldNewCen_KD , K * D * sizeof( float ) ) ;

  EstimSizes( X_ND, C_NK, N, D, K, 
	      N_K, N_KD ) ;

  sts = EstimLaplaceCenters( X_ND, Sort_ND, C_NK, N, D, K, Miss, 
			     N_K, N_KD, oldcent_KD, OldDisp_KD, 
			     OldNewCen_KD, EmptyK_P ) ;

  EstimLaplaceIner( X_ND, C_NK, N, D, K, Miss, 
		    N_K, N_KD, oldcent_KD, OldDisp_KD, OldNewCen_KD, 
		    Iner_KD ) ;

  /* Free local structures */
  GenFree( oldcent_KD ) ;

  return sts ;

  /*???*/
}   /* end of CommonLaplaceDiag() */


/* ------------------------------------------------------------------- */
static void EstimSizes
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  float*        N_K,           /* O : size of a class (K) */
  float*        N_KD           /* O : size of a class and variable (K,D) */  
 )
/* ------------------------------------------------------------------- */
{
  int      h ;             /* current class    0..K-1 */
  int      j ;             /* current variable 0..D-1 */
  int      i ;             /* current object   0..N-1 */

  
  /* For each class h and variable j */
  for ( h = 0 ; h < K ; h ++ ) {

      for ( j = 0 ; j < D ; j ++ ) {

	  /* Initialize the sum[i] quantities to 0 */
	  N_K[ h ]             = 0.0 ;
	  N_KD[ _HJ ]          = 0.0 ;

	  /* For each object i, increment class size and
	     eventually observed data size if xij not nan */
	  for ( i = 0 ; i < N ; i ++ ) {

	      float cih = C_NK[ _IH ] ;
	      float xij = X_ND[ _IJ ] ;

	      N_K[ h ] += cih ;

	      if ( ! isnan( xij ) ) {
		  N_KD[ _HJ ]          += cih ;
	      }
	  } /* End  For each object i (all sums[i] are now done) */
      }
  }  /* End   for each class h and variable j */

}   /* end of EstimSizes() */



/* ------------------------------------------------------------------- */
static StatusET          /* Return status : OK, EMPTY */
EstimLaplaceCenters
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const MissET  Miss,          /* I : how to treat missing data */
  const float*  N_K,           /* I : size of a class (K) */
  const float*  N_KD,          /* I : size of a class and variable (K,D) */  
  const float*  OldCen_KD,     /* I : old centers */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */

  float*        NewCen_KD,     /* O : updated centers (K,D) */
  int*          EmptyK_P       /* O : index of empty class (0 or 1 ..K) */
 )
/* ------------------------------------------------------------------- */
{
  StatusET sts ;           /* return status */

  int      h ;             /* current class    0..K-1 */
  int      j ;             /* current variable 0..D-1 */

  int      imed ;          /* index of median of observed values */
  float    medval ;        /* median value (eventually midway) */
  float    cumwei ;        /* cumulated weights until median observation */


  /* Set no empty class by default */
  (*EmptyK_P) = 0 ;
  sts         = STS_OK ;


  /* For each class h and variable j */
  for ( h = 0 ; h < K ; h ++ ) {

    for ( j = 0 ; j < D ; j ++ ) {

      /* If this class size > 0 */
      if ( N_K[ h ] > EPSILON ) {
	
	/* First, find median position of observed values 
	 */
	ComputeMedian( X_ND, Sort_ND, C_NK, N, D, K, h, j, N_KD[ _HJ ],
		       & imed, & cumwei, & medval ) ;

	/* 
	 * Compute center from observed median, depending on presence
	 * of missing data, REPLACE or IGNORE missing mode 
	 */

	if ( N_KD[ h * D + j ] == N_K[ h ] ) /* no missing data */ {
	  NewCen_KD[ h * D + j ] = medval ;
	}
	else if ( N_KD[ h * D + j ] >= EPSILON ) /* some missing data */ {
	   if ( Miss == MISSING_REPLACE ) {

	     NewCen_KD[ h * D + j ] = 
	       FindMinInerLaplaceEM( X_ND, Sort_ND, C_NK, N, D, K, h, j, 
				     imed, cumwei, OldCen_KD[ h * D + j ], 
				     OldDisp_KD[ h * D + j ],
				     N_KD[ h * D + j ], 
				     N_K[ h ] - N_KD[ h * D + j ] ) ;

	     /* 
		NewCen_KD[ h * D + j ] =  
		( N_KD[ h * D + j ] * medval + 
		( N_K[ h ] - N_KD[ h * D + j ] ) * OldCen_KD[ h * D + j ] )
		/ N_K[ h ] ;
		*/
	   }
	   else  /* assume Miss == MISSING_IGNORE */ {
	     NewCen_KD[ h * D + j ] = medval ;
	   }
	}
	else /* N_KD[ h * D + j ] < EPSILON : all data missing */ {
	  NewCen_KD[ h * D + j ] = OldCen_KD[ h * D + j ] ;
	}

      }
      else /* then this class size == 0 */ {
	NewCen_KD[ h * D + j ] = OldCen_KD[ h * D + j ] ;
	/* signal empty class and which one */
	sts = STS_W_EMPTYCLASS ;
	(*EmptyK_P) = h + 1 ;
      }

    }
  }
  /* end  For each class h and variable j */

  return sts ;
  /*???*/
}   /* end of EstimLaplaceCenters() */



/* ------------------------------------------------------------------- */
static void ComputeMedian
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const int     H,             /* I : current cluster */
  const int     J,             /* I : current variable */
  const float   totwei,        /* I : total weight */
  int*          ImedP,         /* O : index of median of observed values */
  float*        CumweiP,       /* O : cumulated weights until median obs. */
  float*        MedvalP        /* O : median value (eventually midway) */
 )
/* ------------------------------------------------------------------- */
{
  float    halfwei = totwei / 2 ;
  int      i ;             /* current index in ascending order  0..N-1 */
  int      pos ;           /* actual position of current object 0..N-1 */
  int      posmed ;        /* actual position of median 0..N-1 */
  int      posnext ;       /* actual position of next to median 0..N-1 */

  
  /* Scan the data in ascending order (skip NaN values), and find position 
     where the cumulated weights >= half total weight */
  for ( i = 0, 
	  (*CumweiP) = 0.0 ;
	( i < N ) &&
	  (*CumweiP) < halfwei ;
	i ++ ) {

    pos = Sort_ND[ i * D + J ] ;
    if ( !isnan( X_ND[ pos * D + J ] ) ) 
      (*CumweiP) += C_NK[ pos * K + H ] ;

  }  /* At this point, (*CumweiP) >= halfwei */

  (*ImedP) = i - 1 ;
  posmed = Sort_ND[ (*ImedP) * D + J ] ;

  /* Median value depends if = or > */
  if ( (*CumweiP) > halfwei + EPSILON ) /* > : median value is here*/ {

    (*MedvalP) = X_ND[ posmed * D + J ] ;
  }
  else  /* = : median value is midway to next non nan observation */ {

    for ( i = (*ImedP) + 1 ;
	  isnan( X_ND[ Sort_ND[ i * D + J ] * D + J ] ) ||
	    ( C_NK[ Sort_ND[ i * D + J ] * K + H ] < EPSILON ) ;
	  i ++ ) {
    }
    posnext = Sort_ND[ i * D + J ] ;
    (*MedvalP) = 0.5 * ( X_ND[ posmed * D + J ] + X_ND[ posnext * D + J ] ) ;
  }

}   /* end of ComputeMedian() */




/* ------------------------------------------------------------------- */
static float FindMinInerLaplaceEM  /* ret: minimizer of expected inertium */
 (
  const float*  X_ND,          /* I : data matrix (N,D) */
  const int*    Sort_ND,       /* I : sorted index / variable (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const int     H,             /* I : current cluster */
  const int     J,             /* I : current variable */
  const int     Imed,          /* I : median position of observed data */
  const float   Cumwei,        /* I : cumulated weight until Xmed */
  const float   Cen0,          /* I : previous center */
  const float   Disp0,         /* I : previous dispersion */
  const float   Nhj,           /* I : number of observations */
  const float   Nmis           /* I : number of missing data */
 )
/* ------------------------------------------------------------------- */
{
  float   xmed ;       /* observation at median position */
  float   weidif ;     /* (before + med) - after or (med + after) - before */
  float   weibefore ;  /* cumulated weight before median position */
  int     direction ;  /* 1 if xmed < cen0, -1 otherwise */

  float   a ;       /* lower bound of current interval */
  float   b ;       /* upper bound of current interval */
  float   intwei ;  /* cumulated intermediate weight */
  float   d ;       /* current derivative R'(.) */
  float   res ;     /* value to return (minimizes expected inertia) */
  int     i ;       /* current index in ascending order  Ifirst..Ilast */

  int     sts ;     /* 0 : R' < 0 until now
		       1 : R'(a) >= 0
		       2 : R'(a) < 0 and R'(b) >= 0 */

  xmed = X_ND[ Sort_ND[ Imed * D + J ] * D + J ] ;

  if ( xmed < Cen0 - EPSILON ) {
    direction =  1 ;
    weidif    =  Cumwei - ( Nhj - Cumwei ) ;
  }
  else if ( xmed >  Cen0 + EPSILON ) {
    direction =  -1 ;
    weibefore =  Cumwei - C_NK[ Sort_ND[ Imed * K + H ] * K + H ] ;
    weidif    =  (Nhj - weibefore) - weibefore ;
  }
  else   /* xmed == Cen0, no search necessary */ {
    return Cen0 ;
  }


  /* Start at median position */
  a      = xmed ;
  intwei = 0.0 ;
  sts    = 0 ; 

  if ( (d = DerivInerDir( a, weidif, intwei, Nmis, Cen0, Disp0 )) >= 0 ) 
    /* R'(a) >= 0,  minimum at a */ {
    sts = 1 ;  
  }
  else  /* R'(a) < 0, look at b */ {

    for ( i = Imed + direction ;

	  ( sts == 0 ) && 
	    ( (X_ND[ Sort_ND[ i * D + J ] * D + J ] - Cen0) * direction < 0 ) 
	    && ( i >= 0 ) && ( i < N ) ;      /* "safety belt" */

	  i = i + direction ) 
      /* For each intermediate point */ {

      if ( ( ! isnan( X_ND[ Sort_ND[ i * D + J ] * D + J ] ) )
	   && ( C_NK[ Sort_ND[ i * D + J ] * K + H ] > EPSILON ) ) {

	b = X_ND[ Sort_ND[ i * D + J ] * D + J ] ;
	if ( (d = DerivInerDir( b, weidif, intwei, Nmis, Cen0, Disp0 )) >= 0 ) 
	  /* R'(b) >= 0, and R'(a) < 0 */ {

	  sts = 2 ;
	}
	else /* R'(b) < 0 */ {

	  /* Start checking new interval */
	  a = b ;
	  intwei += 2 * C_NK[ Sort_ND[ i * D + J ] * K + H ] ;
	  if ( (d = DerivInerDir( a, weidif, intwei, Nmis, Cen0, Disp0 )) >= 0)
	    /* R'(a) >= 0,  minimum at a */ {

	    sts = 1 ;
	  }
	}
      }
    }    /* end for each intermediate point */
  }

  if ( sts == 0 ) 
    /* R'(b) < 0, R'(a) < 0 => check for b = Cen0 */ {
    b = Cen0 ;
    if ( (d = DerivInerDir( b, weidif, intwei, Nmis, Cen0, Disp0 )) >= 0 ) 
      /* R'(b) >= 0, and R'(a) < 0 */ {

      sts = 2 ;
    }
  }

  switch( sts ) {
  case 0:
    /* R'(a) < 0, R'(b) < 0 => Cen0 is the minimizer */
    res = Cen0 ;
    break ;

  case 1:
    /* R'(a) >= 0 => a is the minimizer */
    res = a ;
    break ;

  default:
    /* R'(a) < 0, R'(b) >= 0 => 
       y in [a,b] minimizing R is solution of :
       | y - Cen0 | = - Disp0 log( 1 - ( weidif - intwei ) / Nmis )
       */
    res = Cen0 + direction * Disp0 * log( 1 - ( weidif + intwei ) / Nmis ) ;
  }

  return res ;
  /*???*/
}   /* end of FindMinInerLaplaceEM() */


/* ------------------------------------------------------------------- */
static float DerivInerDir      /* ret: +/- R'(Y), R expected inertia */ 
 (
  const float   Y,             /* I : point at which to compute R' */
  const float   Weidif,        /* I : "after+med-bef" or "bef+med-after" */
  const float   Intwei,        /* I : cumulated intermediate weight */
  const float   Nmis,          /* I : number of missing data */
  const float   Cen0,          /* I : previous center */
  const float   Disp0          /* I : previous dispersion */
 )
/*\
     This function returns 
       R'(Y)   if Xmed < Cen0
     - R'(Y)   if Xmed > Cen0
     Testing the positiveness of this function thus tests R'(Y) > 0
     in the first case and R'(Y) < 0 in the second case.
\*/
/* ------------------------------------------------------------------- */
{

  float d = Weidif + Intwei 
    - Nmis * ( 1 - exp( - fabs( Y - Cen0 ) / Disp0 ) ) ;

  return d ;

  /*???*/
}   /* end of DerivIner() */




/* ------------------------------------------------------------------- */
static void EstimLaplaceIner
 ( 
  const float*  X_ND,          /* I : data matrix (N,D) */
  const float*  C_NK,          /* I : classification matrix (N,K) */
  const int     N,             /* I : nb of observation vectors */
  const int     D,             /* I : nb of variables */
  const int     K,             /* I : number of clusters */
  const MissET  Miss,          /* I : how to treat missing data */
  const float*  N_K,           /* I : size of a class (K) */
  const float*  N_KD,          /* I : size of a class and variable (K,D) */  
  const float*  OldCen_KD,     /* I : old centers */
  const float*  OldDisp_KD,    /* I : old dispersions (K,D) */
  const float*  NewCen_KD,     /* I : new centers */
  float*        Iner_KD        /* O : inertia of a class/variable (K,D) */  
 )
/* ------------------------------------------------------------------- */
{
  int      h ;             /* current class    0..K-1 */
  int      j ;             /* current variable 0..D-1 */
  int      i ;             /* current object   0..N-1 */

  
  /* For each class h and variable j */
  for ( h = 0 ; h < K ; h ++ ) {

      for ( j = 0 ; j < D ; j ++ ) {

	  /* Initialize the inertium to 0 */
	  Iner_KD[ _HJ ] = 0.0 ;

	  /* For each non nan object i, increment class inertium */
	  for ( i = 0 ; i < N ; i ++ ) {

	      float cih = C_NK[ _IH ] ;
	      float xij = X_ND[ _IJ ] ;

	      if ( ! isnan( xij ) ) {
		Iner_KD[ _HJ ] += cih * fabs( xij - NewCen_KD[ _HJ ] ) ;
	      }

	  } /* End  For each non-nan object i (all sums[i] are now done) */

	  
	  if ( ( N_KD[ _HJ ] < N_K[ h ] ) &&  ( Miss == MISSING_REPLACE ) )
	       /* If some missing data and replace mode */ {

	    Iner_KD[ _HJ ] += 
	      ( N_K[ h ] - N_KD[ _HJ ] ) *
	      ( fabs( OldCen_KD[ _HJ ] - NewCen_KD[ _HJ ] ) +
		OldDisp_KD[ _HJ ] * 
		exp( - fabs( OldCen_KD[ _HJ ] - NewCen_KD[ _HJ ] ) / 
		     OldDisp_KD[ _HJ ] ) ) ; 
	  }

      }
  }  /* End   for each class h and variable j */
  
  /*???*/
}   /* end of EstimLaplaceIner() */




/*+++*/

#if 0


/* ------------------------------------------------------------------- */
StatusET                    /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaP_V_I
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
  int      k ;
  StatusET sts ;


  sts = ParaPkV_I( Cih_nk, DataP, K, MissMode, EmptyK_P, NoiseParaP ) ;

  for ( k = 0 ; k < K ; k ++ )
    {
        NoiseParaP->Pk[ k ] = 1.0 / K ;
    }
  return sts ;
}

/* ------------------------------------------------------------------- */
StatusET ParaPkV_I          /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{

    if ( DataP->NbMiss == 0 )
      return ParaPkV_I_Complete( Cih_nk, DataP, K, MissMode, 
				 EmptyK_P, NoiseParaP ) ;
    else
      return ParaPkV_I_Missing( Cih_nk, DataP, K, MissMode, 
				EmptyK_P, NoiseParaP ) ;

}   /* end of ParaPkV_I() */


/* ------------------------------------------------------------------- */
static StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaPkV_I_Complete
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
    int         h ;         /* class counter : 0 .. K - 1 */
    float       sumh_trTh ; /* sum over h of class h's dispersion matrix's trace */
    float       V ;         /* common volume = sum[h] trTh / (N * Nd) */
    int         nd  = DataP->NbVars ;   /* Nd : nb of variables */
    int         npt = DataP->NbPts ;    /* N : nb of points */

    /* For each class */
    for ( h = 0 , (*EmptyK_P) = 0, sumh_trTh = 0.0 ; 
          h < K ; 
          h ++ )
    {
        int     ipt ;
        int     d ;
        float   sumi_uih ;
        float   sumi_uih_yi2 ;

        /* 'fuzzy' cardinal : sizeh = sum[i]( uih ) */
        /* proportion :       ph    = 1 / K */
        /* mean :             mh(d) = sum[i]( uih * yi(d) ) / sizeh */
        /* trace of its dispersion matrix : 
           trTh = sum[i]( uih * sum[d] yi(d)� ) - sizeh * sum[d] mh(d)� */
        /* volume :           Vh = V = sum[h] trTh / (N * Nd) */

        /* Initialize means to zero */
        for ( d = 0 ; d < nd ; d ++ )
        {
            NoiseParaP->Mk[ ( h * nd ) + d ] = 0 ;
        }

        /* Do cumulated sums over all points */
        for ( ipt = 0, sumi_uih = 0.0, sumi_uih_yi2 = 0.0 ; 
              ipt < npt ; 
              ipt ++ )
        {
            float   sumd_uih_yid2 ;
            float   uih                             = Cih_nk[ ( ipt * K ) + h ] ;

            for ( d = 0, sumd_uih_yid2 = 0.0 ; 
                  d < nd ; 
                  d ++ )
            {
                float   yid   = DataP->PointsM[ (ipt * nd) + d ] ;

                NoiseParaP->Mk[ ( h * nd ) + d ]    += uih * yid ;
                sumd_uih_yid2                       += uih * yid * yid ;
            }

            sumi_uih        += uih ;
            sumi_uih_yi2    += sumd_uih_yid2 ;
        }

        /* Normalize the sums by the cardinal of the class */
        if ( sumi_uih != 0.0 )
        {
            float   sumd_mhd2 ;
            float   trTh ;
            float   invsizeh = 1 / sumi_uih ;

            for ( d = 0, sumd_mhd2 = 0.0 ; 
                  d < nd ; 
                  d ++ )
            {
                float   mhd         = NoiseParaP->Mk[ (h * nd) + d ] ;

                mhd                 *= invsizeh ;
                sumd_mhd2           += mhd * mhd ;

                NoiseParaP->Mk[ (h * nd) + d ]  = mhd ;
            }

            trTh                    = sumi_uih_yi2 - sumi_uih * sumd_mhd2 ;
            sumh_trTh               += trTh ;
        }
        else
        {
            (*EmptyK_P) = h + 1 ;
        }

        NoiseParaP->Pk[ h ] = sumi_uih / npt ;

        /* Set normalized variance matrix to identity */
        SetIdMatrix( nd , & NoiseParaP->Ck[ h*nd*nd ] ) ;  /*V1.03-a*/

    } /* end for ( h = 0 , (*EmptyK_P) = 0 ; h < K ; h ++ ) */

    V = sumh_trTh / ( npt * nd ) ;

    for ( h = 0 ; h < K ; h ++ )
    {
        NoiseParaP->Vk[ h ] = V ;
    }

    if ( (*EmptyK_P) == 0 )
      return STS_OK ;
    else
      return STS_W_EMPTYCLASS ;

}   /* end of ParaPkV_I_Complete() */



/* ------------------------------------------------------------------- */
static StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaPkV_I_Missing
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
  const char* func = "ParaPkV_I_Missing" ;

  int      N = DataP->NbPts ;
  int      D = DataP->NbVars ;

  int      h ;         /* current class    0..K-1 */
  int      j ;         /* current variable 0..D-1 */

  float*   nh_k ;      /* size / class (K) local alloc (unused) */
  float*   nhj_kd ;    /* size / obs class & var (K,D) local alloc (unused)*/
  float*   disphj_kd ; /* dispersion / class & var (K,D) local alloc */
  float    vol ;       /* common volume */

  StatusET sts ;


  /* Allocate local data */
  nh_k      = GenAlloc( K,     sizeof( float ), 1, func, "nh_k" ) ;
  nhj_kd    = GenAlloc( K * D, sizeof( float ), 1, func, "nhj_kd" ) ;
  disphj_kd = GenAlloc( K * D, sizeof( float ), 1, func, "disphj_kd" ) ;

  /* Common computation for any diagonal gaussian model with missing data */
  sts = CommonGaussDiagMissing
    ( DataP->PointsM, DataP->NbPts, DataP->NbVars, Cih_nk, K, MissMode,
      NoiseParaP->Vk, NoiseParaP->Ck, 
      NoiseParaP->Mk, NoiseParaP->Pk, EmptyK_P, nh_k, nhj_kd, disphj_kd ) ;

  /* Compute common volume : 
     sum dispersions over the classes and variables 
     then divide by all/observed data size */
  
  vol = 0.0 ;

  for ( h = 0; h < K ; h ++ )
    {
      /* If not empty class */
      if ( nh_k[ h ] > 0 )
	{
	  for ( j = 0 ; j < D ; j ++ )
	    {
	      vol += disphj_kd[ _HJ ] ;
	    }
	}
    }

  if ( MissMode == MISSING_REPLACE )
    vol /= ( N * D ) ;
  else
    vol /= ( N * D - DataP->NbMiss ) ;

  /* For each class h */
  for ( h = 0; h < K ; h ++ )
    {
      /* Set its volume to the common volume */
      NoiseParaP->Vk[ h ] = vol ;

      /* Set its normalized covariance to identity */
      SetIdMatrix( D, & NoiseParaP->Ck[ h * D * D ] ) ;
    }
      
  /* Free local data */
  GenFree( disphj_kd ) ;
  GenFree( nhj_kd ) ;
  GenFree( nh_k ) ;

  return sts ;
}   /* end of ParaPkV_I_Missing() */



/* ------------------------------------------------------------------- */
StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaP_VkI
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
  StatusET sts ;
  int k ;

  sts = ParaPkVkI( Cih_nk, DataP, K, MissMode, EmptyK_P, NoiseParaP ) ;

  for ( k = 0 ; k < K ; k ++ )
    {
        NoiseParaP->Pk[ k ] = 1.0 / K ;
    }
  return sts ;
}


/* ------------------------------------------------------------------- */
StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaPkVkI
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{

    if ( DataP->NbMiss == 0 )
      return ParaPkVkI_Complete( Cih_nk, DataP, K, MissMode, 
				 EmptyK_P, NoiseParaP ) ;
    else
      return ParaPkVkI_Missing( Cih_nk, DataP, K, MissMode, 
				EmptyK_P, NoiseParaP ) ;

}   /* end of ParaPkV_I() */



/* ------------------------------------------------------------------- */
static StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaPkVkI_Complete
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
    int         k ;

    /* For each class */
    for ( k = 0 , (*EmptyK_P) = 0 ; k < K ; k ++ )
    {
        int     ipt ;                  
        int     id ;
        int     nd  = DataP->NbVars ;
        int     npt = DataP->NbPts ;    /* = n nb points */
        float   sizek = 0.0 ;            /* = nk fuzzy cardinal [0,n] */
        float   trWk = 0.0 ;

        /* Compute its 'fuzzy' cardinal : nk = sum[i=1,n]( uik ) */
        /* Compute its proportion : pk = nk / n */
        /* Compute its mean (mk) : mk = sum[i=1,n]( uik * Y(i,:) ) / nk */
        /* Compute its intra-class inertia matrix : 
           Wk = sum[i=1,n]( uik * ( Y(i,:) - mk ) * ( Y(i,:) - mk )' ) / nk */
        /* Compute its volume : Vk = tr(Wk) / d */

        /* Initialize means to zero */
        for ( id = 0 ; id < nd ; id ++ )
        {
            NoiseParaP->Mk[ (k * nd) + id ] = 0 ;
        }

        /* Do cumulated sums over all points */
        for ( ipt = 0 ; ipt < npt ; ipt ++ )
        {
            float   uik = Cih_nk[ (ipt * K) + k ] ;

            sizek += uik ;
            for ( id = 0 ; id < nd ; id ++ )
            {
                float   yid = DataP->PointsM[ (ipt * nd) + id ] ;
                float   mkd = NoiseParaP->Mk[ (k * nd) + id ] ;

                mkd         = mkd + uik * yid ;
                trWk        = trWk + uik * yid * yid ;

                NoiseParaP->Mk[ (k * nd) + id ] = mkd ;
            }
        }

        /* Normalize the sums by the cardinal of the class */
        if ( sizek != 0.0 )
        {
            trWk /= sizek ;

            for ( id = 0 ; id < nd ; id ++ )
            {
                float   mkd = NoiseParaP->Mk[ (k * nd) + id ] ;

                mkd         = mkd / sizek ;
                trWk        = trWk - ( mkd * mkd ) ;

                NoiseParaP->Mk[ (k * nd) + id ] = mkd ;
            }
            NoiseParaP->Vk[ k ] = trWk / nd ;
        }
        else
        {
            (*EmptyK_P) = k + 1 ;
        }

        NoiseParaP->Pk[ k ] = sizek / npt ;

        /* Set normalized variance matrix to identity */
        SetIdMatrix( nd , & NoiseParaP->Ck[ k*nd*nd ] ) ;  /*V1.03-a*/

    } /* end for ( k = 0 , (*EmptyK_P) = 0 ; k < K ; k ++ ) */

    if ( (*EmptyK_P) == 0 )
      return STS_OK ;
    else
      return STS_W_EMPTYCLASS ;

}   /* end of ParaPkVkI_Complete() */


/* ------------------------------------------------------------------- */
static StatusET             /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
ParaPkVkI_Missing
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  NoiseParaT    *NoiseParaP /* O : estimated parameters */
 ) 
/* ------------------------------------------------------------------- */
{
  const char* func = "ParaPkVkI_Missing" ;

  int      N = DataP->NbPts ;
  int      D = DataP->NbVars ;

  int      h ;         /* current class    0..K-1 */
  int      j ;         /* current variable 0..D-1 */

  float*   nh_k ;      /* size / class (K) local alloc */
  float*   nhj_kd ;    /* size / obs class and variable (K,D) local alloc */
  float*   disphj_kd ; /* dispersion / class and variable (K,D) local alloc */

  float    sumj_nhj ;  /* size of obs class */

  StatusET sts ;


  /* Allocate local data */
  nh_k      = GenAlloc( K,     sizeof( float ), 1, func, "nh_k" ) ;
  nhj_kd    = GenAlloc( K * D, sizeof( float ), 1, func, "nhj_kd" ) ;
  disphj_kd = GenAlloc( K * D, sizeof( float ), 1, func, "disphj_kd" ) ;

  /* Common computation for any diagonal gaussian model with missing data */
  sts = CommonGaussDiagMissing
    ( DataP->PointsM, DataP->NbPts, DataP->NbVars, Cih_nk, K, MissMode,
      NoiseParaP->Vk, NoiseParaP->Ck, 
      NoiseParaP->Mk, NoiseParaP->Pk, EmptyK_P, nh_k, nhj_kd, disphj_kd ) ;

  /* For each class h */
  for ( h = 0; h < K ; h ++ )
    {
      /* If not empty class */
      if ( nh_k[ h ] > 0 )
	{
	  /* Reestimate its volume: sum dispersions over the variables 
	     then divide by all/observed data size */
	  sumj_nhj = 0 ;
	  NoiseParaP->Vk[ h ] = 0.0 ;  /* use as sum first */

	  for ( j = 0 ; j < D ; j ++ )
	    {
	      sumj_nhj            += nhj_kd[ _HJ ] ;
	      NoiseParaP->Vk[ h ] += disphj_kd[ _HJ ] ;
	    }

	  if ( MissMode == MISSING_REPLACE )
	    NoiseParaP->Vk[ h ] /= nh_k[ h ] ;
	  else
	    NoiseParaP->Vk[ h ] /= sumj_nhj ;
	}
      else /* class h is empty */
	{
	  NoiseParaP->Vk[ h ] = 0 ;
	}

      /* Set its normalized covariance to identity */
      SetIdMatrix( D, & NoiseParaP->Ck[ h * D * D ] ) ;
    }

  /* Free local data */
  GenFree( disphj_kd ) ;
  GenFree( nhj_kd ) ;
  GenFree( nh_k ) ;

  return sts ;
}



/*V1.05-d*/
/* ------------------------------------------------------------------- */
static StatusET                /* Returns : OK, W_EMPTYCLASS or E_MEMORY */
CommonGaussDiagMissing
 (
  const float*  Xij_nd,        /* I : data matrix (N,D) */
  int           N,             /* I : nb of observation vectors */
  int           D,             /* I : nb of variables */
  const float*  Cih_nk,        /* I : classification matrix (N,K) */
  int           K,             /* I : number of clusters */
  MissET        Miss,          /* I : how to treat missing data */
  const float*  OldVh_k,       /* I : old volumes (K) */
  const float*  OldChjj_kdd,   /* I : old normalized covariances (K,D,D) */

  float*        OldNewMhj_kd,  /* I/O : old then updated means (K,D) */
  float*        NewPh_k,       /* O : updated proportions (K) */
  int*          EmptyK_P,      /* O : index of empty class (0 or 1 ..K) */
  float*        Nh_k,          /* O : size of a class (K) */
  float*        Nhj_kd,        /* O : size of a class and variable (K,D) */
  float*        Disphj_kd      /* O : dispersion of a class/variable (K,D) */
 )
/*\

   This function does the common computation that is needed at the
   M-step, in all the diagonal Gaussian models with missing data.

   The routine takes as input a data matrix Xij_nd (possibly
   containing NaN values), a classification matrix Cih_nk, and the
   kind of missing data processing Miss (REPLACE or IGNORE) ; and also
   the old volumes OldVh_k, covariance matrices OldChjj_kdd, and means
   OldNewMhj_kd of the classes, to be used in the REPLACE mode.

   It produces as output the reestimated mean vectors in OldNewMhj_kd;
   reestimated proportions in NewP_k; the last encountered empty class
   in EmptyK_P (or 0 if no empty class); the {sum[i] cih} in Nh_k, the
   {sum[i] cih rij} in Nhj_kd (number of present observations in class h
   and variable j); and the "missing data dispersion" of each
   class/variable in Disphj_kd (i.e. the missing data equivalent of
   {sum[i] Cih (xij - mhj)^2}).

\*/
/* ------------------------------------------------------------------- */
{
  const char* func = "CommonGaussDiagMissing" ;

  float*   sumdata_hj_kd ;    /* sum[i] cih rij xij   (K,D), local alloc */
  float*   sumsquare_hj_kd ;  /* sum[i] cih rij xij^2 (K,D), local alloc */
  float*   oldmean_hj_kd ;    /* saved previous means (K,D), local alloc */
  int      h ;                /* current class    0..K-1 */
  int      j ;                /* current variable 0..D-1 */
  int      i ;                /* current object   0..N-1 */
  StatusET sts ;              /* return status */


  sts = STS_OK ;

  /* Allocate local structures */
  sumdata_hj_kd   = GenAlloc( K * D, sizeof(float), 1, func, "sumdata" ) ;
  sumsquare_hj_kd = GenAlloc( K * D, sizeof(float), 1, func, "sumsquare" );
  oldmean_hj_kd   = GenAlloc( K * D, sizeof(float), 1, func, "oldmean" ) ;

  /* Save old mean value */
  memcpy( oldmean_hj_kd , OldNewMhj_kd , K * D * sizeof( float ) ) ;

  /* Set no empty class by default */
  (*EmptyK_P) = 0 ;

  /* For each class h */
  for ( h = 0 ; h < K ; h ++ )
    {
      /* For each variable j */
      for ( j = 0 ; j < D ; j ++ )
	{
	  /* Initialize the sum[i] quantities to 0 */
	  Nh_k[ h ]               = 0.0 ;
	  Nhj_kd[ _HJ ]           = 0.0 ;
	  sumdata_hj_kd[ _HJ ]    = 0.0 ;
	  sumsquare_hj_kd[ _HJ ]  = 0.0 ;
	  OldNewMhj_kd[ _HJ ]     = 0.0 ; /* use Mhj as sum first */

	  /* For each object i */
	  for ( i = 0 ; i < N ; i ++ )
	    {
	      float cih = Cih_nk[ _IH ] ;
	      float xij = Xij_nd[ _IJ ] ;

	      /* Increment class size */
	      Nh_k[ h ] += cih ;

	      /* If this object's j variable is observed */
	      if ( ! isnan( xij ) )
		{
		  /* Increment all sum[i] of observed quantities */
		  Nhj_kd[ _HJ ]           += cih ;
		  sumdata_hj_kd[ _HJ ]    += cih * xij ;
		  sumsquare_hj_kd[ _HJ ]  += cih * xij * xij ;
		}
	      /* Else xij is missing, don't change sum[i] of observed data */
	    }
	  /* End  For each object i (all sums[i] are now done) */

	  /* If this class size > 0 */
	  if ( Nh_k[ h ] > 0 )
	    {
	      /* 
	       * Compute means Mhj and dispersions Dhj from computed sums,
	       * depending on REPLACE or IGNORE missing mode 
	       */
	      if ( Miss == MISSING_REPLACE )
		{
		  /* 
		   * new_mhj = 
		   *  ( sum[i] cih rij xij + ( nh - nhj ) old_mhj ) / nh 
		   *
		   * disp_hj = 
		   *  "observed dispersion"_hj +
		   *  ( nh - nhj ) [ (new_mhj - old_mhj)^2 + old_vh old_chjj ]
		   */
		  float obsdisphj ;

		  OldNewMhj_kd[ _HJ ] = 
		    ( sumdata_hj_kd[ _HJ ] + 
		      ( Nh_k[ h ] - Nhj_kd[ _HJ ] ) * oldmean_hj_kd[ _HJ ] )
		    / Nh_k[ h ] ;

		  /*
		   * "observed dispersion"_hj = 
		   *      sum[i] cih rij xij^2 - nhj new_mhj^2 
		   */
		  obsdisphj = sumsquare_hj_kd[ _HJ ] - 
		    Nhj_kd[ _HJ ] * sqr( OldNewMhj_kd[ _HJ ] ) ;

		  Disphj_kd[ _HJ ] = 
		    obsdisphj +
		    ( Nh_k[ h ] - Nhj_kd[ _HJ ] ) *
		    ( sqr( OldNewMhj_kd[ _HJ ] - oldmean_hj_kd[ _HJ ] ) +
		      OldVh_k[ h ] * OldChjj_kdd[ h * D * D + j * D + j ] ) ;
		}
	      else  /* assume Miss == MISSING_IGNORE */
		{
		  float obsdisphj ;

		  /* 
		   * new_mhj = 
		   *  . if  nhj <> 0,  ( sum[i] cih rij xij ) / nhj  
		   *  . if  nhj == 0,  old_mhj
		   *
		   * disp_hj = 
		   *  "observed dispersion"_hj 
		   */

		  if ( Nhj_kd[ _HJ ] > 0 )
		    OldNewMhj_kd[ _HJ ] = sumdata_hj_kd[ _HJ ] / Nhj_kd[ _HJ ];
		  else
		    OldNewMhj_kd[ _HJ ] = oldmean_hj_kd[ _HJ ] ;

		  /*
		   * "observed dispersion"_hj = 
		   *      sum[i] cih rij xij^2 - nhj new_mhj^2 
		   */
		  obsdisphj = sumsquare_hj_kd[ _HJ ] - 
		    Nhj_kd[ _HJ ] * sqr( OldNewMhj_kd[ _HJ ] ) ;

		  Disphj_kd[ _HJ ] = obsdisphj ;
		}

	    } 
	  else /* then this class size == 0 */
	    {
	      /* signal empty class and which one */
	      sts = STS_W_EMPTYCLASS ;
	      (*EmptyK_P) = h + 1 ;
	    }
	}
      /* end  For each variable j */

      /* Update proportions as  ph = nh / n */
      NewPh_k[ h ] = Nh_k[ h ] / N ;
    }
  /* end  For each class h */


  /* Free local structures */
  GenFree( oldmean_hj_kd ) ;
  GenFree( sumsquare_hj_kd ) ;
  GenFree( sumdata_hj_kd  ) ;

  return sts ;

}  /*  end of CommonGaussDiagMissing() */





/* ------------------------------------------------------------------- */
void SetIdMatrix ( int Nd , float* M )
/* Set matrix M to identity */ /*V1.03-a*/
/* ------------------------------------------------------------------- */
{
  int l , c ;  /* line and column counters : 0..Nd-1 */

  for ( l = 0 ; l < Nd ; l ++ )
    {
      for ( c = 0 ; c < Nd ; c ++ )
        {
          if ( c == l )
            M[ (l * Nd) + c ] = 1.0 ;
          else
            M[ (l * Nd) + c ] = 0.0 ;
        }
    }
}


#endif

/* ------------------------------------------------------------------- */


/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ./._nem_mod.h                                                                                       000755  000765  000765  00000000312 11541743333 012674  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_mod.h                                                                                           000755  000765  000024  00000002123 11541743333 012614  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\
    Prototypes of nem_noi.c exported functions

    1.05-a    12-JAN-1997  MissMode in ParaP*V*I
    1.05-b    17-JAN-1997  EmptyK_P and StatusET return in ParaP*V*I
    1.06-a    28-JUN-1998  GetDensityFunc <- nem_alg.c and del DensPkVkI
    1.06-b    28-JUN-1998  New EstimPara and del ParaP*V*I
\*/

#include "nem_typ.h"    /* NoiseParaT, ... */

/*V1.06-a*/
int GetDensityFunc  /* STS_OK or STS_E_FUNCARG */
        (
            const ModelSpecT  *SpecP,           /* I */
            CompuDensFT**     CompuDensFP       /* O */
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


                                                                                                                                                                                                                                                                                                                                                                                                                                             ./._nem_nei.c                                                                                       000755  000765  000765  00000000312 11541743334 012664  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_nei.c                                                                                           000755  000765  000024  00000011355 11541743334 012613  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\
    nem_nei.c

    Fonctions de voisinage.
\*/
#include    "nem_nei.h"     /* prototypes */
#include    <stdio.h>       /* fprintf */
#include    <string.h>      /* memcpy */


static int      GetNeighNone                 /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) ;


static int      GetNeighImage                /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) ;


static int      GetNeighIrreg                /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) ;



/* ------------------------------------------------------------------- */
int GetSpatialFunc  /* STS_OK or STS_E_FUNCARG */
    (
     TypeET          SpatialType,    /* I */
     GetNeighFT**    GetNeighFP      /* O */
    ) 
/* ------------------------------------------------------------------- */
{
    switch( SpatialType )
    {
        case TYPE_NONSPATIAL :
             *GetNeighFP = GetNeighNone ;
        return STS_OK ;

        case TYPE_SPATIAL :
             *GetNeighFP = GetNeighIrreg ;
        return STS_OK ;

        case TYPE_IMAGE :
             *GetNeighFP = GetNeighImage ;
        return STS_OK ;

        default :
             *GetNeighFP = NULL ;
                fprintf( stderr, "GetSpatialFuncs bad arg : Type = %d\n",
                         SpatialType ) ;
        return STS_E_FUNCARG ;
    }
}   /* end of GetSpatialFuncs() */




/* ------------------------------------------------------------------- */
int         GetNeighNone                 /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) 
/* ------------------------------------------------------------------- */
{
    return 0 ;
}


/* ------------------------------------------------------------------- */
int         GetNeighImage                /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) 
/* ------------------------------------------------------------------- */
{
    int     in ;
    int     rnn ;       /* real number of neighbours */
    int     nbn     = NeighDataP->Image.NbNeigh ;
    INeighT *neiV   = NeighDataP->Image.NeighsV ;
    int     nl      = NeighDataP->Image.Nl ;
    int     nc      = NeighDataP->Image.Nc ;
    int     line    = Ipt / nc ;
    int     col     = Ipt % nc ;

    if ( nbn > PtNeighsP->NbNeigh )
    {
        nbn = PtNeighsP->NbNeigh ;
    }

    for ( in = 0, rnn = 0 ; in < nbn ; in ++ )
    {
        int l = line + neiV[ in ].Dl ;
        int c = col + neiV[ in ].Dc ;

        if ( ( 0 <= l ) && ( l < nl ) && ( 0 <= c ) && ( c < nc ) )
        {
            PtNeighsP->NeighsV[ rnn ].Index  = l * nc + c ;
            PtNeighsP->NeighsV[ rnn ].Weight = neiV[ in ].Weight ;
            rnn ++ ;
        }
    }

    return rnn ;
}

/* ------------------------------------------------------------------- */
int         GetNeighIrreg                /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) 
/* ------------------------------------------------------------------- */
{
    PtNeighsT*  ptnP = &(NeighDataP->PtsNeighsV[ Ipt ]) ;
    int         nbn = ptnP->NbNeigh ;

    if ( nbn > PtNeighsP->NbNeigh )
    {
        nbn = PtNeighsP->NbNeigh ;
    }

    memcpy( PtNeighsP->NeighsV, ptnP->NeighsV, nbn * sizeof( NeighT ) );

    return nbn ;        
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                                                                                                                                                                                                                                                                                   ./._nem_nei.h                                                                                       000755  000765  000765  00000000312 11541743334 012671  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_nei.h                                                                                           000755  000765  000024  00000000476 11541743334 012622  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\
Vers-mod  Date         Who Description
1.06-a    28-JUN-1998  MD  GetSpatialFunc instead of GetNeighImage, ...
\*/

#include "nem_typ.h"    /* NeighDataT */


int GetSpatialFunc  /* STS_OK or STS_E_FUNCARG */
    (
     TypeET          SpatialType,    /* I */
     GetNeighFT**    GetNeighFP      /* O */
    ) ;

                                                                                                                                                                                                  ./._nem_rnd.c                                                                                       000755  000765  000765  00000000312 11541743334 012674  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_rnd.c                                                                                           000755  000765  000024  00000004535 11541743334 012625  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\

    NEM_RND.C

    Programme NEM (Neighborhood EM) : routines de tirage aleatoire

    Van Mo DANG       Avril 96

Vers-mod  Date         Who  Description

1.04-a    10-OCT-1997  MD   add RandomPermutationAlgo()
1.04-b    05-NOV-1997  MD   use random/srandom instead of lrand48/srand48
1.04-c    12-JAN-1997  MD   add RandomReal()
\*/


#include <sys/types.h>   /* time_t */
#include <time.h>        /* time() */
#include <stdlib.h>      /* srand48(), lrand48() */

#include "nem_rnd.h"

#define MAXRAND   0x7fffffff   /* 2**31 - 1 = maximum value of random() */


void  RandomSeedByTime( void ) 
{
    time_t   t ;

    t = time( 0 ) ;

#ifdef __TURBOC__
    srand( (unsigned) t ) ;
#else
    srandom( t ) ; /*V1.04-b*/
#endif
}


int   RandomInteger( int Mini, int Maxi ) 
{
    int       span ;
    long int  nrandom ;
    int       result ;


    if ( Mini >= Maxi )
      {
        return Maxi ;
      }

    span = Maxi - Mini + 1 ;

#ifdef __TURBOC__
    nrandom = rand( ) ;
#else
    nrandom = random( ) ; /*V1.04-b*/
#endif

    result = ( (int) ( nrandom % span ) ) + Mini  ;

    return  result ;
}


/*V1.05-a*/
float   RandomFloat( float Mini, float Maxi )
{
    float     span ;
    long int  nrandom ;
    float     result ;


    if ( Mini >= Maxi )
      {
        return Maxi ;
      }

    span = Maxi - Mini ;

#ifdef __TURBOC__
    nrandom = rand( ) ;
#else
    nrandom = random( ) ; /*V1.04-b*/
#endif

    result = ( (float) nrandom / MAXRAND ) * span + Mini  ;

    return  result ;
}





/* =========================== */

void RandomPermutationAlgo( int* TabV , int Nb )      /*V1.04-a*/

/* =========================== */
{
  int icou ;
  int iech ;
  int valech ;

  for ( icou = 0 ; icou < Nb ; icou ++ )
    {
      iech = RandomInteger( 0 , Nb - 1 ) ;

      valech       = TabV[ iech ] ;

      TabV[ iech ] = TabV[ icou ] ;

      TabV[ icou ] = valech ;
    }
}



/* Tests 

#include <stdio.h>

main( int argc, char *argv[] )
{
  int       n ;
  int       mini, maxi ;
  int       i ;

  if ( argc < 4 )
    {
      printf( "Au moins 3 args\n" );
      return 2 ;
    }
  n = atoi( argv[ 1 ] ) ;

  mini = atoi( argv[ 2 ] ) ;
  maxi = atoi( argv[ 3 ] ) ;

  RandomSeedByTime( ) ;

  for ( i = 0 ; i < n ; i ++ )
    {
      int r = RandomInteger( mini , maxi  ) ;

      fprintf( stdout, "%8d  ", r ) ;
    }

  fprintf( stdout, "\n" ) ;
  return 0 ;
}
*/
                                                                                                                                                                   ./._nem_rnd.h                                                                                       000755  000765  000765  00000000312 11541743334 012701  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_rnd.h                                                                                           000755  000765  000024  00000000542 11541743334 012624  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         void  RandomSeedByTime( void ) ;
int   RandomInteger( int Mini, int Maxi ) ;
float   RandomFloat( float Mini, float Maxi ) ;      /*V1.05-a*/
void  RandomPermutationAlgo( int* TabV , int Nb ) ;  /*V1.04-a*/

/* external library random functions */
#if !defined(__GO32__) /* djgpp has its own random */
int srandom( unsigned seed );
long random();
#endif
                                                                                                                                                              ./._nem_rnd2.h                                                                                      000755  000765  000765  00000000312 11541743335 012764  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_rnd2.h                                                                                          000755  000765  000024  00000000571 11541743335 012711  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         void  RandomSeedByTime( void ) ;
int   RandomInteger( int Mini, int Maxi ) ;
float   RandomFloat( float Mini, float Maxi ) ;      /*V1.05-a*/
void  RandomPermutationAlgo( int* TabV , int Nb ) ;  /*V1.04-a*/

/* external library random functions */
#if !defined(__GO32__) /* djgpp has its own random */
int srandom( unsigned seed );
long random();
#endif /* RANDOM */
                                                                                                                                       ./._nem_typ.h                                                                                       000755  000765  000765  00000000312 11541743335 012733  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_typ.h                                                                                           000755  000765  000024  00000035030 11541743335 012656  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         #ifndef NEM_TYP_H
#define NEM_TYP_H

/*\
    NEM_TYP.H

    Programme NEM (Neighborhood EM) : tous les types echanges entre modules

    Van Mo DANG       Janvier 96


    Vers-mod  Date         Description

    V1.03-a   01-NOV-1996  Add log-likelihood criterion in CriterT
    V1.04-a   04-OCT-1997  Add BetaET for beta estimation mode
    V1.04-b   04-OCT-1997  Add BtaStep and BtaMode in NemParaT
    V1.04-c   10-OCT-1997  Add CvThres in NemParaT
    V1.04-d   10-OCT-1997  Add Seed in NemParaT
    V1.04-e   10-OCT-1997  Type OrderET and add VisitOrder in NemParaT
    V1.04-f   13-OCT-1997  Debug and RefName in NemParaT
    V1.04-g   13-OCT-1997  Error in CriterT
    V1.04-h   13-OCT-1997  CritET and Crit in NemParaT
    V1.05-a   12-JAN-1998  NbMissing in DataT, MisModeET, MissMode in NemParaT
    V1.05-b   12-JAN-1998  MissMode in parameters of EstimNoiseFT
    V1.05-c   16-JAN-1998  EstimNoiseFT : EmptyK_P as parameter, return status
    V1.05-d   05-FEB-1998  Add Markov fuzzy pseudo-like. (CritET and CriterT)
    V1.06-a   17-JUN-1998  NoiseModel -> StatModel, structure changed
    V1.06-b   23-JUN-1998  LEN_LINE <- nem_exe.c
    V1.06-c   30-JUN-1998  Add define EPSILON
    V1.06-d   03-AUG-1998  Add INIT_MIXFIX and INIT_MIXINI
    V1.06-e   10-SEP-1998  Add SiteUpdate in NemParaT and UpdET
    V1.06-f   10-SEP-1998  Add TieRule in NemParaT and TieET
    V1.06-g   15-SEP-1998  Add BtaPsGrad in NemParaT and BtaPsGradT struct
    V1.06-h   20-SEP-1998  Add ErrInfo in CriterT and ErrInfoT struct
    V1.06-i   21-SEP-1998  Change CvTest in CriterT and add CvemET
    V1.06-j   30-NOV-1998  Add EPSILON_INV
    V1.06-k   01-DEC-1998  FkP double* instead of *float in compudensft
    V1.07-a   26-FEB-1999  FAMILY_BERNOULLI added
\*/

/*
 *  Constant definitions
 */

#define    TRUE             1
#define    FALSE            0

#define    LEN_NOISEMODEL   20
#define    LEN_FILENAME     100
#define    LEN_LINE         500  /*V1.06-b*/


#define    MAX_PTS          1000000L
#define    MAX_VARS         100

#define    CONV_THRES       0.001

#define    EXT_OUTNAMEHARD  ".cf"
#define    EXT_OUTNAMEFUZZY ".uf"
#define    EXT_MFNAME       ".mf"
#define    EXT_LOGNAME      ".log"

#define    EPSILON          1e-20  /* to check for FP zero or equality */
#define    EPSILON_INV      1e20   /* multiply by this for small floats */

/*
 *  Enumerated types
 */

typedef enum { 
        STS_OK = 0,
        STS_I_DONE,
        STS_W_EMPTYCLASS,
        STS_E_ARG,
        STS_E_MEMORY, 
        STS_E_FILEIN,
        STS_E_FILEOUT,
        STS_E_FILE,
        STS_E_FUNCARG
        } 
        StatusET ;


typedef enum { 
        ALGO_NEM ,    
        ALGO_NCEM ,
	ALGO_GEM ,
	ALGO_NB
        } 
        AlgoET;


typedef enum { 
        BETA_FIX ,        /* Use given fixed beta */
        BETA_PSGRAD ,     /* Estimate beta using pseudo-likelihood gradient */
        BETA_HEUD,        /* Estimate beta using Hathaway heuristic */  
        BETA_HEUL,        /* Estimate beta using likelihood heuristic */  
	BETA_NB
        } 
        BetaET;           /* Beta estimation mode */ /*V1.04-a*/


typedef enum { 
        CRIT_U ,    /* Use NEM criterion */
	CRIT_M ,    /* Use markovian fuzzy log pseudo-likelihood */ /*V1.05-d*/
	CRIT_D ,    /* Use Hathaway criterion */
	CRIT_L ,    /* Use Mixture log-likelihood */
	CRIT_NB	  
        } 	  
        CritET;     /* Which criterion to choose local max */ /*V1.04-h*/



typedef enum { 
        TYPE_SPATIAL, 
        TYPE_IMAGE, 
        TYPE_NONSPATIAL,
        TYPE_NB
        } 
        TypeET;


/*V1.06-a*/
#if 0
typedef enum {
	   MODEL_P_VkI,
	   MODEL_PkVkI,
	   MODEL_P_V_I,
	   MODEL_PkV_I,
	   MODEL_NB
	   }
	   ModelET ;
#endif

typedef enum { 
  FAMILY_NORMAL,
  FAMILY_LAPLACE,
  FAMILY_BERNOULLI,
  FAMILY_NB
} 
FamilyET;

typedef enum { 
  DISPER___,  /* 1 dispersion for all classes and variables (V.I) */
  DISPER_K_,  /* 1 dispersion for each class, same in all variables (Vk.I) */
  DISPER__D,  /* 1 dispersion for each variable, same in all classes (V.B) */
  DISPER_KD,  /* 1 dispersion for each class and variable (VkBk) */
  DISPER_NB
} 
DisperET;

typedef enum { 
  PROPOR__,  /* equal proportion for all classes = 1/K (P) */
  PROPOR_K,  /* 1 proportion for each class (Pk) */
  PROPOR_NB
} 
ProporET;


typedef enum {
        FORMAT_HARD,
        FORMAT_FUZZY,
        FORMAT_NB
        }
        FormET ;

typedef enum {
        INIT_SORT,
        INIT_RANDOM,
	INIT_MIXINI,      /* EM mixture estimate at start */
        INIT_MIXFIX,      /* EM mixture estimate with fixed value */
        INIT_FILE,
	INIT_LABEL,
        INIT_NB
        }
        InitET ;

typedef enum
{
  MISSING_REPLACE,  /* Replace missing statistics with expectation as in EM */
  MISSING_IGNORE,   /* Ignore missing statistics as in CEM */
  MISSING_NB
}
MissET ;   /*V1.05-a*/


typedef enum {
        NEIGH_FOUR,
        NEIGH_FILE,
        NEIGH_NB
        }
        NeighET ;

typedef enum {
        ORDER_DIRECT,
        ORDER_RANDOM,
        ORDER_NB
        }
        OrderET ;       /*V1.04-e*/


typedef enum {
        UPDATE_SEQ,
        UPDATE_PARA,
        UPDATE_NB
        }
        UpdET ;         /*V1.06-e*/


typedef enum {
        TIE_RANDOM,
        TIE_FIRST,
        TIE_NB
        }
        TieET ;         /*V1.06-f*/


typedef enum {
        CVTEST_NONE,
        CVTEST_CLAS,
	CVTEST_CRIT,
        CVTEST_NB
        }
        CvemET ;        /*V1.06-g*/



/*
 *  Structured types
 */

typedef struct
{
    int         NbPts ;       /* number of observation vectors */
    int         NbVars ;      /* number of variables */
    int         NbMiss ;      /* number of missing data 0..Npts*NbVars */
                              /*V1.05-a*/
    float       *PointsM ;    /* observations (NbPts,NbVars) to allocate */
    int         *LabelV ;     /* fixed labels (NbPts) to allocate: 0..k */
    int         *SiteVisitV ; /* site to visit (NbPts) to allocate: 0..Npts-1*/
    int         *SortPos_ND ; /* PointsM[ SortPos_ND[i*D+d]*D+d ] : +++*/
}
DataT ;         /* Matrix of observed data (each line = 1 point) */


/*V1.06-g*/
typedef struct
{
  int     NbIter ;     /* Max number of iterations of gradient ascent */
  float   ConvThres ;  /* Convergence if gradient <  threshold * N */
  float   Step ;       /* >0 : bta += grad*(step/N), 0: bta += grad/dsec */
  int     RandInit ;   /* 1 = random initial beta, 0 = specified by -b */
}
BtaPsGradT ;   /* parameters of beta pseudo-likelihood gradient estimation */


typedef struct
{
    AlgoET  Algo ;      /* Type of algorithm */
  /*V1.06-a*/
#if 0
     float   Beta ;      /* context weight : >= 0 */
     BetaET  BtaMode ;   /* type of beta estimation */  /*V1.04-b*/
#endif
    float   BtaHeuStep ;/* step of beta for heuristic estimation */ /*V1.04-b*/
    float   BtaHeuMax ; /* maximum beta for heuristic */
    float   BtaHeuDDrop ;/* drop of Hathaway slope for beta heuristic */
    float   BtaHeuDLoss ;/* proportion of Hathaway loss for beta heuristic */
    float   BtaHeuLLoss ;/* proportion of likelihood loss for beta heuristic */
    BtaPsGradT BtaPsGrad ; /* parameters of beta gradient estimation */
    CritET  Crit ;      /* criterion to choose local max */ /*V1.04-h*/
    float   CvThres ;   /* convergence threshold */    /*V1.04-c*/
    CvemET  CvTest ;    /* which convergence test to use */
    int     DoLog ;     /* TRUE if log file requested */
    int     NbIters ;   /* nb of iterations for NEM */
    int     NbEIters ;  /* nb of iterations for E-step */
    int     NbRandomInits ;  /* nb of random initializations */
    long    Seed ;      /* random generator seed */   /*V1.04-d*/
    FormET  Format ;    /* output file format (hard or fuzzy) */
    InitET  InitMode ;  /* initialization mode (histogram, random, file) */
    MissET  MissMode ;  /* how to process missing statistics */ /*V1.05-a*/
    int     SortedVar ; /* variable to be sorted : 0..NbVars */
    NeighET NeighSpec ; /* neighborhood specification */
    OrderET VisitOrder ;/* order of visit at E-step */ /*V1.04-e*/
    UpdET   SiteUpdate ;/* site update scheme at E-step */
    TieET   TieRule ;   /* rule for equal probabilities when computing MAP */
    int     Debug ;     /* TRUE if in debug mode */    /*V1.04-f*/
    char    OutBaseName[ LEN_FILENAME + 1 ] ; /* base name of output file */
    char    OutName[ LEN_FILENAME + 1 ] ;   /* name of output file */
    char    LogName[ LEN_FILENAME + 1 ] ;   /* name of log file ("" = no log) */
    char    StartName[ LEN_FILENAME + 1 ] ; /* name of initial partition file */
    char    NeighName[ LEN_FILENAME + 1 ] ; /* name of neighborhood file */
    char    LabelName[ LEN_FILENAME + 1 ] ; /* name of fixed labels file */
    char    RefName[ LEN_FILENAME + 1 ] ;   /* name of reference labels file *//*V1.04-f*/
}
NemParaT ;      /* NEM running parameters */


typedef struct
{
    int     Dl ;    /* neighbour shift in line :     -2 .. 2 */
    int     Dc ;    /* neighbour shift in column :   -2 .. 2 */
    float   Weight ; /* neighbour weight >= 0 */
}
INeighT ;       /* one neighbour (in image configurations) */

typedef struct
{
    int     Nl ;        /* Image number of lines (height) : > 0 */
    int     Nc ;        /* Image number of columns (width) : > 0 */
    int     NbNeigh ;   /* nb of allocated neighbours : >= 0 */
    INeighT *NeighsV ;  /* to be allocated : pixel's neighbours */
}
ImageNeighT ;   /* neighbourhood system (in image configurations) */

typedef struct
{
    int         Index ; /* index of neighbour : 0 .. Npt-1 */
    float       Weight ;/* weight of neighbour (default : 1.0) */
}
NeighT ;        /* one neighbour (in non-image configurations) */

typedef struct
{
    int     NbNeigh ;   /* nb of allocated neighbours */
    NeighT  *NeighsV ;  /* to be allocated : point's neighbours */
}
PtNeighsT;      /* one point's neighbours (in non-image configuration) */

typedef union
{
    ImageNeighT Image ;
    PtNeighsT   *PtsNeighsV ;   /* to be allocated : all points' neighbours */
}
NeighDataT ;    /* generic neighbourhood system */

typedef struct
{
    NeighDataT  NeighData ; /* generic neighbourhood system */
    int         MaxNeighs ; /* maximum number of neighbours (>= 0) */
    TypeET      Type ;      /* type of spatial configuration */
}
SpatialT ;


/*V1.06-a*/

#if 0
typedef struct
{
    float   *Pk ;   /* proportions (K) */   /* to be allocated */
    float   *Vk ;   /* volumes (K) */       /* to be allocated */
    float   *Ck ;   /* shapes (d*d*K) */    /* to be allocated */
    float   *Mk ;   /* means (d*K) */       /* to be allocated */
}
NoiseParaT ;


typedef struct
{
    ModelET     ModelNum ;
    int         Nk ;
    NoiseParaT  NoisePara ;
}
NoiseModelT ;
#endif

typedef struct
{
  int       K ;                /* number of classes */
  FamilyET  ClassFamily ;
  DisperET  ClassDisper ;
  ProporET  ClassPropor ;
  BetaET    BetaModel ;
}
ModelSpecT ;  /* Model specification */

typedef struct
{
  float     Beta ;
  float*    Center_KD ;  /* Center in each class and variable (K*D) */
  float*    Disp_KD ;    /* Dispersion in each class and variable (K*D) */
  float*    Prop_K ;     /* Proportion of each class (in ]0,1[) */

  float*    NbObs_K ;    /* Nb of observations in each class (K) */
  float*    NbObs_KD ;   /* Nb of observations in each class/variable (K*D) */
  float*    Iner_KD ;    /* Inertia = sum_i cik * Dist(xid, mkd) (K*D) */
}
ModelParaT ;  /* Model parameters */

typedef struct
{
  float*     DispSam_D ; /* Dispersion of whole sample in each variable (D) */
  float*     MiniSam_D ; /* Minimum of whole sample in each variable (D) */
  float*     MaxiSam_D ; /* Maximum of whole sample in each variable (D) */
}
SampleDesT ;  /* Sample description */

typedef struct
{
  ModelSpecT   Spec ;
  ModelParaT   Para ;
  SampleDesT   Desc ;
}
StatModelT ;  /* Model description */


typedef struct
{
  int      Kc ;            /* # of classes in computed classification */
  int      Kr ;            /* # of classes in reference classification */
  int      Km ;            /* max( Kr, Kc ) */
  int      Kmfac ;         /* Km! */
  TieET    TieRule ;       /* same value as in NemPara */
  float*   Refclas_N_Kr ;  /* reference classification : 0/1 */
  int*     Perm_Kmfac_Km ; /* all permutation of Km classes */
}
ErrinfoT ;    /* Classification error information */


typedef struct
{
  float*   Agree_Km_Km ;   /* #common elements between found and ref classes */
  float*   Loclas_N_Kc ;   /* local copy of computed classification */
  int      Ibestpermut ;   /* index of best agreement permutation 0..Km!-1 */
  float    Errorrate ;     /* # misclassified objects of best permut / N */
}
ErrcurT ;    /* Currently computed classification error */


typedef struct
{
  float    D ; /* hathaway crit. D = sum[i]sum[k] cik (log pkfki-log cik) */
  float    G ; /* geog. cohesion G = sum[i]sum[j]sum[k] cik cjk wij */
  float    U ; /* NEM maximized criterion U = D + 0.5 * beta * G */
  float    M ; /* markovian fuzzy class. like. M = D + beta * G - Z */
  float    L ; /* mixture likelihood crit. L = sum[i] log sum[k] pkfki ) */
  float    Z ; /* log pseudo-l. Z =-sum[i]log(sum[k]e(bta*sum[j~i]wij cjk)) */
  ErrinfoT Errinfo ; /* information to compute error */   /*V1.06-h*/
  ErrcurT  Errcur ;  /* current error rate */   /*V1.06-h*/
} /*V1.05-d*/
CriterT ;       /*V1.03-a*/




typedef StatusET            /* ret : OK, W_EMPTYCLASS or E_MEMORY *//*V1.05-c*/
EstimNoiseFT          
 (
  const float   *Cih_nk,    /* I : classification matrix (N,K) */
  const DataT   *DataP,     /* I : observed points */
  int           K,          /* I : number of classes */
  MissET        MissMode,   /* I : how to treat missing data */ /*V1.05-b*/

  int           *EmptyK_P,  /* O : which empty class (1..K) or 0 *//*V1.05-c*/
  StatModelT    *StatModelP /* O : estimated parameters */
 ) ;

typedef int CompuDensFT         /* ret : 0 if OK, -1 if zero density */
        (
            int                Nd,      /* I : point dimension */
            int                Ik,      /* I : class number : 0..Nk-1 */
            const ModelParaT*  ParaP,   /* I : noise parameters *//*V1.06-a*/
            const float*       XV,      /* I : point (dim d) */
            double*            FkP,     /* O : density for class Ik */
            float*             LogFkP   /* O : log of density */
        ) ;


typedef int GetNeighFT         /* ret : nb of neighbours */
                ( 
                  int               Ipt ,           /* I : index of point */
                  const NeighDataT  *NeighDataP,    /* I : neighborhood data */
                  PtNeighsT         *PtNeighsP      /* O : indices/weights of neighbours */
                ) ;


#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ./._nem_ver.c                                                                                       000755  000765  000765  00000000312 11541743336 012707  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_ver.c                                                                                           000755  000765  000024  00000002152 11541743336 012631  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\
    NEM_VER.C

    Informations about software successive versions.
    To be updated at each modification.
\*/

#include "nem_ver.h"  /* Exported prototypes */

const char *NemVersionStrC = "1.07" ; /* Current version of NEM */

void PrintVersions( FILE* F )         /* Describes successive versions */
{
    fprintf( F , "\n" ) ;
    fprintf( F , " Vers  Date      Description\n" ) ;
    fprintf( F , " 0.00  01.02.96  First version released on WEB\n" ) ;
    fprintf( F , " 1.00  31.05.96  Random start, image defaults, long help\n" ) ;
    fprintf( F , " 1.01  27.06.96  NCEM, sequential E-step\n" ) ;
    fprintf( F , " 1.02  17.10.96  Partially known labels\n" ) ;
    fprintf( F , " 1.03  02.10.97  Init centers from known labels, log is optional, augmented .mf\n" ) ;
    fprintf( F , " 1.04  11.01.98  Estimation of beta, Gibbsian EM, longer help\n" ) ;
    fprintf( F , " 1.05  09.04.98  Missing data, initialization modified\n" ) ;
    fprintf( F , " 1.06  26.02.99  Pseudo-likelihood beta, Laplace distributions\n" ) ;
    fprintf( F , " 1.07  08.04.99  Bernoulli distributions\n" ) ;
    fprintf( F , "\n" ) ;
}

                                                                                                                                                                                                                                                                                                                                                                                                                      ./._nem_ver.h                                                                                       000755  000765  000765  00000000312 11541743336 012714  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_ver.h                                                                                           000755  000765  000024  00000000271 11541743336 012636  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         #include <stdio.h> /* FILE */


extern const char *NemVersionStrC ;    /* Current version of NEM software */

extern void PrintVersions( FILE* F ) ; /* Describes successive versions */
                                                                                                                                                                                                                                                                                                                                       ./._lib_io.c                                                                                        000755  000765  000765  00000000312 11541743315 012506  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      lib_io.c                                                                                            000755  000765  000024  00000046010 11541743315 012431  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\
    lib_io.c

    Routines for input/output and memory allocation

    June 96

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    Vers-mod  Date         Description

    1.05-a    26-JAN-1998  Remove temporarily added alloc_exit to lib_io.c
\*/

#include <stdio.h>   /* FILE */
#include <string.h>  /* strcpy */
#include <stdlib.h>  /* FILE */

#include "lib_io.h"   /* prototypes of exported functions */

#define TRUE       1
#define FALSE      0

#define LEN_LINE   500
#define LEN_FIELD  30



/* ------------------------------------------------------------------- */
int   ReadOpeningComments  /* 0 : OK, -1 : can't open, 1 : LenComment small */
      (
       const char*   FileName,    /* I : name of file to open */
       const char*   MarkerS,     /* I : comment marker of beginning */
       int           LenComment,  /* I : length of allocated CommentS */
       FILE**        FP ,         /* O : opened file, NULL if impossible */
       char*         CommentS     /* O : read comment */
      )
/* ------------------------------------------------------------------- */
{
    int     mlen = strlen( MarkerS ) ;
    int     iscomment ;
    int     nblines ;
    int     sts ;
    char    line[ LEN_LINE + 1 ] ;
    int     iline ;

    if ( ( (*FP) = fopen( FileName , "r" ) ) == NULL )
      return -1 ;

    /* Read comment lines */
    strcpy( CommentS, "" ) ;
    for ( iscomment = TRUE , nblines = 0 , sts = 0 ; 
          iscomment && (! feof( (*FP) )) ; 
          nblines ++ )
    {
        if ( fgets( line, LEN_LINE, (*FP) ) != NULL )
        {
            iscomment = ( strstr( line , MarkerS ) == line ) ;
            if ( iscomment )
            {
              if ( sts == 0 )
                {
                  if ( (int)( strlen( CommentS ) + strlen( &line[ mlen ] ) ) > 
                      LenComment )
                    sts = 1 ;
                
                  strncat( CommentS, &(line[ mlen ]) , LenComment ) ;
                }
            }
        }
    }
    nblines -- ;

    fclose( (*FP) ) ;

    /* Skip comment lines */
    (*FP) = fopen( FileName, "r" ) ;
    for ( iline = 0 ; iline < nblines ; iline ++ )
      {
        fgets( line , LEN_LINE , (*FP) ) ;
      }

    return sts ;

} /* end of ReadOpeningComments() */

/* ------------------------------------------------------------------- */
int AskFileToRead           /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of file to read */
        char*       NameF   /* O : name of (readable) file */
    )
/* ------------------------------------------------------------------- */
{
    int     exists ;    /* TRUE if last entered file name exists */
    int     nbask ;     /* nb of repeated asks */

    for ( exists = FALSE , nbask = 1 ; 
          ( ! exists ) && ( nbask <= MAX_ASK ) ; 
          nbask ++ )
    {
        if ( nbask == 1 )
            printf( "Name of  %s  file  (RETURN to quit) : ", Desc ) ;

        gets( NameF ) ;

        if ( strlen(NameF) != 0 )
        {
            FILE*   f ;

            if ( ( f = fopen( NameF , "r" ) ) != NULL )
            {
                fclose( f ) ;
                exists = TRUE ;
            }
            else
            {
                printf( " '%s' does not exist. " , NameF ) ;
                if ( nbask < MAX_ASK )
                    printf( "Please type again : " ) ;
                else
                    printf( "\n" ) ;
                exists = FALSE ;
            }
        }
        else nbask = MAX_ASK ;
    }

    if ( exists )
       return 0 ;   /* last file typed in exists */
    else
        return -1 ; /* MAX_ASK or more unsuccessful tries */

} /* end of AskFileToRead() */

/* ------------------------------------------------------------------- */
int AskFileToWrite          /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of file to write */
        int         Conf ,  /* I : TRUE if ask confirmation for overwrite */
        char*       NameF   /* O : name of (readable) file */
    )
/* ------------------------------------------------------------------- */
{
    int     writeok ; /* TRUE if entered file name creation ok */
    int     nbask ;   /* nb of repeated asks */

    for ( writeok = FALSE , nbask = 1 ; 
          ( ! writeok ) && ( nbask <= MAX_ASK ) ; 
          nbask ++ )
    {
        FILE*   f ;
        int     accept ;

        printf( "Name of  %s  file to create : " , Desc ) ;
        gets( NameF ) ;

        if ( strlen(NameF) != 0 )
        {
            /* By default, accept file overwriting ; this acceptation
               will be unvalidated if confirmation asked, and
               file already exists, and user refuses to overwrite file
            */
            accept = TRUE ;
            if ( ( Conf ) && ( ( f = fopen( NameF , "r" ) ) != NULL ) )
            {
                char c ;

                fclose( f ) ;   /* after successful open "r" , close file */

                printf( "File %s already exists. Overwrite it ? (y/n/q) " ,
                        NameF ) ;

                c = getchar() ;
                getchar() ; /* to empty buffer */
                switch( c )
                {
                    case 'y' : accept = TRUE ; break ;
                    case 'q' : accept = FALSE ; nbask = MAX_ASK ; break ;
                    default  : accept = FALSE ;
                }
            }

            if ( accept )
            {
                /* Try to create file, if fails : directory may be invalid */
                if ( ( f = fopen( NameF , "w" ) ) != NULL )
                {
                    fclose( f ) ;   /* after open "w" , close file */
                    remove( NameF ) ; /* remove file */
                    writeok = TRUE ;
                }
                else
                {
                    printf( " Cannot create '%s' (check name/permission)\n" ,
                            NameF ) ;
                }
            }
        }
        else    nbask = MAX_ASK ;
    }

    if ( writeok )
       return 0 ;   /* last file typed in could be created */
    else
        return -1 ; /* MAX_ASK or more unsuccessful tries */

} /* end of AskFileToWrite() */


/* ------------------------------------------------------------------- */
int AskInteger              /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of number to type in */
        int         Def ,   /* I : default value */
        int         Min ,   /* I : minimum value */
        int         Max ,   /* I : maximum value */
        int*        NbReadP /* O : number read */
    )
/* ------------------------------------------------------------------- */
{
    int     numberok ; /* TRUE if entered file name creation ok */
    int     nbask ;    /* nb of repeated asks */

    for ( numberok = FALSE , nbask = 1 ; 
          ( ! numberok ) && ( nbask <= MAX_ASK ) ; 
          nbask ++ )
    {
        char    stringread[ 132 + 1 ] ;

        printf( "Enter  %s  ( %d <= n <= %d )  [%d]  : " , 
                Desc , Min , Max , Def ) ;
        gets( stringread ) ;

        if ( strlen( stringread ) != 0 )
        {
            if ( ( sscanf( stringread , "%d" , NbReadP ) == 1 ) &&
                 ( Min <= (*NbReadP) ) && ( (*NbReadP) <= Max )  )
            {
                numberok = TRUE ;
            }
            else printf( " Invalid number\n" ) ;
        }
        else
        {
            (*NbReadP) = Def ;
            numberok = TRUE ;
        }
    }

    if ( numberok )
       return 0 ;   /* last file typed in could be created */
    else
        return -1 ; /* MAX_ASK or more unsuccessful tries */

} /* end of AskInteger() */


/* ------------------------------------------------------------------- */
int AskFloat                /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of number to type in */
        float       Def ,   /* I : default value */
        float       Min ,   /* I : minimum value */
        float       Max ,   /* I : maximum value */
        float*      NbReadP /* O : number read */
    )
/* ------------------------------------------------------------------- */
{
    int     numberok ; /* TRUE if entered file name creation ok */
    int     nbask ;    /* nb of repeated asks */

    for ( numberok = FALSE , nbask = 1 ; 
          ( ! numberok ) && ( nbask <= MAX_ASK ) ; 
          nbask ++ )
    {
        char    stringread[ 132 + 1 ] ;

        printf( "Enter  %s  ( %g <= n <= %g )  [%g]  : " , 
                Desc , Min , Max , Def ) ;
        gets( stringread ) ;

        if ( strlen( stringread ) != 0 )
        {
            if ( ( sscanf( stringread , "%f" , NbReadP ) == 1 ) &&
                 ( Min <= (*NbReadP) ) && ( (*NbReadP) <= Max )  )
            {
                numberok = TRUE ;
            }
            else printf( " Invalid number\n" ) ;
        }
        else
        {
            (*NbReadP) = Def ;
            numberok = TRUE ;
        }
    }

    if ( numberok )
       return 0 ;   /* last file typed in could be created */
    else
        return -1 ; /* MAX_ASK or more unsuccessful tries */

} /* end of AskFloat() */


/* ------------------------------------------------------------------- */
int CountTokens               /* ret : nb of tokens in Line */
    (
        const char* Line ,    /* I : line to analyze */
        const char* SeparS    /* I : separator between tokens */
    ) 
/* ------------------------------------------------------------------- */
{
    int            NbTokens ;                /* to be returned */
    static char    myline[ LEN_LINE + 1 ] ;  /* static to avoid stacking */
    int            len ;
    char*          p;

    /* Copy to local string, and eventually strip off newline char */
    strncpy( myline , Line , LEN_LINE ) ;
    len = strlen( myline ) ;
    if ( myline[ len - 1 ] == '\n' )
       myline[ len - 1 ] = '\0' ;

    for ( NbTokens = 0 , p = strtok( myline , SeparS ) ;
          p != NULL ;
          p = strtok( NULL , SeparS ) )
    {
        NbTokens ++ ;
    }

    return NbTokens ;

} /* end of CountTokens() */


/* ------------------------------------------------------------------- */
int CountLinesColumns          /* 0/1 = same/dif. nb of col, -1 = problem */
    (
        const char* NameF ,    /* I : file name to analyze */
        const char* SeparS ,   /* I : separator between columns */
        int*        MinColP ,  /* O : minimum number of columns*/
        int*        MaxColP ,  /* O : maximum number of columns*/
        int*        NbLinesP   /* O : number of lines */
    )
/* ------------------------------------------------------------------- */
{
    FILE*   finp ;                  /* input file handler */
    char    line[ LEN_LINE + 1 ] ;  /* last line read */
    int     nblines ;
    int     mincols=0 , maxcols=0 ;
    int     noteq ;


    if ( ( finp = fopen( NameF , "r" ) ) == NULL )
    {
        printf( "Error : can't open file %s\n" , NameF ) ;
        return -1 ;
    }

    /* Read first line */
    nblines = 0 ;
    if ( fgets( line , LEN_LINE , finp ) != NULL )
    {
        maxcols = CountTokens( line , SeparS ) ;
        mincols = maxcols ;
        if ( maxcols > 0 )
           nblines ++ ;
    }

    /* Read following lines */
    noteq = FALSE ;
    while ( ! feof( finp ) )
    {

        if ( fgets( line , LEN_LINE , finp ) != NULL )
        {
            int nbcols = CountTokens( line , SeparS ) ;

            if ( nbcols > 0 )
            {
                nblines ++ ;
                if ( nbcols != maxcols )
                {
                    noteq = TRUE ;
                    if ( nbcols > maxcols )
                       maxcols = nbcols ;
                    else
                        mincols = nbcols ;
                }
            }
        }
    }

    (*NbLinesP) = nblines ;
    (*MinColP) = mincols ;
    (*MaxColP) = maxcols ;

    /* Close file */
    fclose( finp ) ;

    if ( noteq )
       return 1 ;
    else
        return 0 ;

} /* end of CountLinesColumns() */


/* ------------------------------------------------------------------- */
int ReadSelectedColumns        /* 0 = OK , -1 = problem */
    (
        const char* NameF ,    /* I : name of file to read */
        int         Npt ,      /* I : number of lines to read */
        int         Ntot ,     /* I : total number of columns per line */
        int         Nsel ,     /* I : number of selected columns */
        const int*  SelCol ,   /* I : selected columns [ at least Nsel ] */
        float*      PtsM       /* O : points read [ Npt * Nsel ] */
    )
/* ------------------------------------------------------------------- */
{
    FILE*   finp ;  /* input file handler */
    int     ok ;    /* TRUE if file format is correct */
    int     i ;     /* current line :      0..Npt-1 */
    int     c ;     /* current column :    0..Ntot-1 */
    int     sel ;   /* current selection : 0..Nsel-1 */
    char    field[ LEN_FIELD + 1 ] ;    /* last token read */


    /* Open file */
    if ( ( finp = fopen( NameF , "r" ) ) == NULL )
    {
        printf( " Error : can't open file %s\n" , NameF ) ;
        return -1 ;
    }

    /* Read field by field */
    for ( i = 0 , ok = TRUE ; ( i < Npt ) && ok ; i ++ )
    {
        for ( c = 0 ; ( c < Ntot ) && ok ; c ++ )
        {
            if ( fscanf( finp , "%s" , field ) == 1 )
            {
                float x ;
                int   isfloat = ( sscanf( field , "%f" , &x ) == 1 ) ;

                for ( sel = 0 ; ( sel < Nsel ) && ok ; sel ++ )
                {
                    if ( SelCol[ sel ] == c )
                    {
                        if ( isfloat )
                        {
                            PtsM[ ( i * Nsel ) + sel ] = x ;
                        }
                        else
                        {
                          printf( " In '%s', [%d,%d] = '%s' not a number\n" ,
                                    NameF , i + 1 , c + 1 , field ) ;
                          ok = FALSE ;
                        }
                    }
                }
            } /* end - if read successful for element [i,c] */
            else 
            {
                printf( " File '%s', cannot read line %d, column %d\n" ,
                        NameF , i + 1 , c + 1 ) ;
                ok = FALSE ;
            }
        } /* end - for each column c */
    } /* end - for each line i */

    /* Close file */
    fclose( finp ) ;

    if ( ok )
        return 0 ;
    else
        return -1 ;

} /* end of ReadSelectedColumns() */


/* ------------------------------------------------------------------- */
/* Test of routine ReadOpeningComments 

#define FNAME         "t.dat"
#define LEN_COMMENT   1000

main()
{
  int   sts ;
  FILE* f ;
  char  comS[ LEN_COMMENT + 1 ] ;
  char  line[ LEN_LINE + 1 ] ;

  sts = ReadOpeningComments( FNAME , "//" , LEN_COMMENT , &f , comS ) ;

  if ( sts != -1 )
    {
      printf( "Comments %s of file %s :\n%s\n" , 
             (sts == 1) ? "(shortened)" : "",  FNAME , comS ) ;
      printf( "Remaining lines : \n" ) ;

      while( !feof( f ) )
        {
          if ( fgets( line , LEN_LINE , f ) != NULL ) 
            printf( line ) ;
        }
      fclose( f ) ;
    }
  return 0 ;
}

*/

/* ------------------------------------------------------------------- */
/* Test of routine AskFileToRead 

main()
{
    int     sts ;
    char    fname[ LEN_FILE + 1 ] ;

    sts = AskFileToRead( "test" , fname ) ;

    printf( "*** AskFileToRead returned file name '%s' (status %d)\n" ,
            fname , sts ) ;
    return sts ;
}

*/


/* ------------------------------------------------------------------- */
/* Test of routine AskFileToWrite

main()
{
    int     sts ;
    char    fname[ LEN_FILE + 1 ] ;

    sts = AskFileToWrite( "test" , TRUE, fname ) ;

    printf( "*** AskFileToWrite returned file name '%s' (status %d)\n" ,
            fname , sts ) ;
    return sts ;
}

*/

/* ------------------------------------------------------------------- */
/* Test of routine AskInteger 
#include <stdlib.h>
main( int argc , char *argv[] )
{
    int     sts ;
    int     n ;

    if ( argc < 4 ) return 1 ;

    sts = AskInteger( "test number" , atoi( argv[1] ) , atoi( argv[2] ) ,
                      atoi( argv[3] ) , & n ) ;

    printf( "*** AskInteger returned '%d' (status %d)\n" ,
            n , sts ) ;
    return sts ;
}

*/

/* ------------------------------------------------------------------- */
/* Test of routine AskFloat 
#include <stdlib.h>
main( int argc , char *argv[] )
{
    int     sts ;
    float   x ;

    if ( argc < 4 ) return 1 ;

    sts = AskFloat( "test number" , atof( argv[1] ) , atof( argv[2] ) ,
                      atoi( argv[3] ) , & x ) ;

    printf( "*** AskFloat returned '%f' (status %d)\n" ,
            x , sts ) ;
    return sts ;
}

*/

/* ------------------------------------------------------------------- */
/* Test of routine CountTokens : 2 command line args = args of function 
main( int argc , char *argv[] )
{
    int     sts ;

    if ( argc < 3 ) return 1 ;

    sts = CountTokens( argv[1] , argv[2] ) ;

    printf( "*** CountTokens( '%s' , '%s' )  returned '%d'\n" , 
            argv[1] , argv[2] , sts ) ;
    return sts ;
}
*/

/* ------------------------------------------------------------------- */
/* Test of routine CountLinesColumns : 
   2 command line args = args of function 

main( int argc , char *argv[] )
{
    int     sts ;
    int     minc , maxc , nbl ;

    if ( argc < 3 ) return 1 ;

    sts = CountLinesColumns( argv[1] , argv[2] , &minc , &maxc , &nbl ) ;

    printf( "*** CountLinesColumns( '%s' , '%s' ) returned '%d'\n" , 
            argv[1] , argv[2] , sts ) ;
    printf( "*** minc = %d   maxc = %d    nbl = %d\n" , minc,maxc,nbl ) ;

    return sts ;
}

*/


/* ------------------------------------------------------------------- */
/* Test of routine ReadSelectedColumns : 
   command line args : 
   1 = filename , 2 = nb sel. col , 3, 4, ... = sel. col

#include <stdlib.h>
main( int argc , char *argv[] )
{
    int     sts ;
    int     minc , nbl , nbc , nbsel , s ;
    int     selcV[ 5 ] ;
    float*  ptM ;

    if ( argc < 4 ) return 1 ;
    nbsel = atoi( argv[2] ) ;
    if ( nbsel > 5 ) return 1 ;
    for ( s = 0 ; s < nbsel ; s ++ )
        selcV[ s ] = atoi( argv[ 3 + s ] ) - 1 ;

    sts = CountLinesColumns( argv[1] , " " , &minc , &nbc , &nbl ) ;
    if ( sts != 0 )
    { printf( "Error : CountLinesColumns returned %d\n" , sts ) ;
      return 2 ; }

    if ( ( ptM = malloc( nbl * nbsel * sizeof( float ) ) ) == NULL )
       return 3 ;

    sts = ReadSelectedColumns( argv[1] , nbl , nbc , nbsel , selcV , ptM ) ;

    printf( "*** ReadSelectedColumns( '%s' , %d , %d , %d , [%d,%d,%d] ) returned '%d'\n" , 
            argv[1] , nbl , nbc , nbsel , selcV[ 0 ] , selcV[ 1 ] , 
            selcV[ 2 ] , sts ) ;

    {
        int i , s ;

        for ( i = 0 ; i < nbl ; i ++ )
        {
            printf( "*** Pt %d : " , i + 1 ) ;
            for ( s = 0 ; s < nbsel ; s ++ )
                printf( "  %4.2f" , ptM[ i * nbsel + s ] ) ;
            printf( "\n" ) ;
        }
    }

    free( ptM ) ;
    return sts ;
}

*/


/* ------------------------------------------------------------------- */
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ./._lib_io.h                                                                                        000755  000765  000765  00000000312 11541743317 012515  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      lib_io.h                                                                                            000755  000765  000024  00000010235 11541743317 012440  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\
    lib_io.h

    Prototypes of routines for input/output

    June 96

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    Vers-mod  Date         Description

    1.05-a    26-JAN-1998  Add ExitET
\*/

#include <stdio.h>   /* FILE */

#define LEN_FILE     132     /* maximum length of file name */
#define MAX_ASK      5       /* maximum number of repeated asked inputs */

typedef enum
{

  EXIT_OK,           /* Good:    program achieved processing normally */
  EXIT_W_RESULT,     /* Warning: program achieved with unusable result */
  EXIT_E_ARGS,       /* Error:   user gave wrong argument syntax */ 
  EXIT_E_FILE,       /* Error:   files not found or wrong format */
  EXIT_E_MEMORY,     /* Error:   program ran out of memory */
  EXIT_E_SYSTEM,     /* Error:   a system call failed */
  EXIT_E_BUG,        /* Error:   internal program inconsistency */
  EXIT_NB

} ExitET ;

/* ------------------------------------------------------------------- */
int   ReadOpeningComments  /* 0 : OK, -1 : can't open, 1 : LenComment small */
      (
       const char*   FileName,    /* I : name of file to open */
       const char*   MarkerS,     /* I : comment marker of beginning */
       int           LenComment,  /* I : length of allocated CommentS */
       FILE**        FP ,         /* O : opened file, NULL if impossible */
       char*         CommentS     /* O : read comment */
      ) ;

/* ------------------------------------------------------------------- */
int AskFileToRead           /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of file to read */
        char*       NameF   /* O : name of (readable) file */
    ) ;

/* ------------------------------------------------------------------- */
int AskFileToWrite          /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of file to write */
        int         Conf ,  /* I : TRUE if ask confirmation for overwrite */
        char*       NameF   /* O : name of (readable) file */
    ) ;

/* ------------------------------------------------------------------- */
int AskInteger              /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of number to type in */
        int         Def ,   /* I : default value */
        int         Min ,   /* I : minimum value */
        int         Max ,   /* I : maximum value */
        int*        NbReadP /* O : number read */
    ) ;

/* ------------------------------------------------------------------- */
int AskFloat                /* ret : 0 if OK, -1 if problem */
    (
        const char* Desc ,  /* I : description of number to type in */
        float       Def ,   /* I : default value */
        float       Min ,   /* I : minimum value */
        float       Max ,   /* I : maximum value */
        float*      NbReadP /* O : number read */
    ) ;

/* ------------------------------------------------------------------- */
int CountLinesColumns          /* 0/1 = same/dif. nb of col, -1 = problem */
    (
        const char* NameF ,    /* I : file name to analyze */
        const char* SeparS ,   /* I : separator between columns */
        int*        MinColP ,  /* O : minimum number of columns*/
        int*        MaxColP ,  /* O : maximum number of columns*/
        int*        NbLinesP   /* O : number of lines */
    ) ;

/* ------------------------------------------------------------------- */
int CountTokens               /* ret : nb of tokens in Line */
    (
        const char* Line ,    /* I : line to analyze */
        const char* SeparS    /* I : separator between tokens */
    ) ;

/* ------------------------------------------------------------------- */
int ReadSelectedColumns        /* 0 = OK , -1 = problem */
    (
        const char* NameF ,    /* I : name of file to read */
        int         Npt ,      /* I : number of lines to read */
        int         Ntot ,     /* I : total number of columns per line */
        int         Nsel ,     /* I : number of selected columns */
        const int*  SelCol ,   /* I : selected columns [ at least Nsel ] */
        float*      PtsM       /* O : points read [ Npt * Nsel ] */
    ) ;
                                                                                                                                                                                                                                                                                                                                                                   ./._exemain.c                                                                                       000755  000765  000765  00000000312 11541743312 012674  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      exemain.c                                                                                           000755  000765  000024  00000001443 11541743312 012620  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\

    main.c

    A main() interface to compile a standalone program which is
    executable from the operating system.

    Jan 1998

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

Vers-mod  Date         Who  Description

1.05-a    25-JAN-1998  MD   Create to make an interface to operating system
1.05-b    30-JAN-1998  MD   Prototype of called mainfunc() in mainfunc.h
\*/


#include "mainfunc.h"   /* mainfunc() */

/* ------------------------------------------------------------------- */
int main( int argc, const char *argv[] )
/*\
    Interface to a main routine.
\*/
/* ------------------------------------------------------------------- */
{
  return mainfunc( argc, argv ) ;
}
/* ~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                                                                                                                                                                                                                             ./._exememo.c                                                                                       000755  000765  000765  00000000312 11541743313 012706  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      exememo.c                                                                                           000755  000765  000024  00000002547 11541743313 012640  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\
    exememo.c

    Specific routines for operating system memory allocation

    Jan 1998

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    Vers-mod  Date         Description

    1.05-a    26-JAN-1998  Creation
\*/


#include <stdlib.h>  /* calloc, free, size_t */
#include <stdio.h>   /* stderr */

/* ------------------------------------------------------------------- */
void* GenAlloc
(
 size_t       nelem,        /* I : number of elements to allocate */ 
 size_t       elsize,       /* I : size in bytes of each element */
 int          doexit,       /* I : 1 if failure exits, 0 if only return NULL */
 const char*  where,        /* I : name of calling function */
 const char*  what          /* I : name of allocated object */
)
{
  void *result = calloc (nelem , elsize);
  if ( result != NULL )
    {
      return result;
    }
  else
  {
    fprintf(stderr, "Fatal: in %s, no memory for %s (%ld elements size %ld)\n",
	    where, what, nelem, elsize);
    if ( doexit )
      exit( EXIT_FAILURE );
    else
      return result ;
  }
}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
void GenFree( void* ptr )
{
  if ( ptr != NULL )
    free( ptr ) ;
}
/* ------------------------------------------------------------------- */
                                                                                                                                                         ./._mexmain.c                                                                                       000755  000765  000765  00000000312 11541743320 012703  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      mexmain.c                                                                                           000755  000765  000024  00000013470 11541743320 012632  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\

    mexmain.c

    A mexFunction() interface to compile a matlab routine which is
    executable from matlab.

    Jan 1998

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

Vers-mod  Date         Who  Description

1.05-a    25-JAN-1998  MD   Create to make an interface to matlab
1.05-b    30-JAN-1998  MD   Prototype of called mainfunc() in mainfunc.h
\*/


#include "mex.h"
#include "mainfunc.h"   /* mainfunc() */

#include <stdio.h>      /* printf, ... */
#include <string.h>     /* strlen, ... */


#define  SEPAR_S    " \t"
#define  CMDNAME_S  "mexprog"


/* Local functions */

static void make_args
(
 const char* cmd,      /* I : command line to analyse */
 int*        argc_p,   /* O : number of arguments */
 char***     argv_p    /* O : array of argc arguments (allocated) */
) ;

static int fetch_args  /* returns number of counted strings */
( 
 const char* cmd,  /* I : command line to analyze */
 int         nmax, /* I : number of allocated args in argv */
 char**      argv  /* O : table [0..nmax-1] strings (each string alloc here)*/
) ;



/* 
 *   Function called from matlab.
 *   A string parameter is expected = the command line string
 *   to invoke the executable program.
 */

void mexFunction(
                 int nlhs, Matrix *plhs[],
                 int nrhs, Matrix *prhs[]
		 )
{
  int    len ;   /* length of given argument string */
  char*  arg_s ; /* given argument string - allocated locally */
  char*  cmd ;   /* cmd + given argument string - allocated locally */

  int    argc ;  /* number of arguments + 1 in string */
  char** argv ;  /* array of the arguments (+ "command") - allocated locally */
  int    iarg ;  /* counter from 0 to argc-1 to free elts of argv */

  int    sts ;   /* status returned from mainfunc */

  /* If args given in a string */
  if ( nrhs >= 1 )
    {
      /* Allocate and get string */
      len = mxGetN( prhs[0] ) + strlen( CMDNAME_S ) + 1 ;
      arg_s = mxCalloc( len, sizeof( char ) ) ;
      if ( arg_s == NULL )
	mexErrMsgTxt( "Could not allocate string of args" ) ;
      mxGetString( prhs[0], arg_s, len ) ;

      cmd = mxCalloc( strlen( arg_s ) + strlen( CMDNAME_S ) + 1 + 1, 
		      sizeof( char ) ) ;
      if ( cmd == NULL )
	mexErrMsgTxt( "Could not allocate 2nd string of args" ) ;
      strcpy( cmd, CMDNAME_S ) ;
      strcat( cmd, " " ) ;
      strcat( cmd, arg_s ) ;

      /* Translate string into mainfunc args */
      make_args( cmd, & argc, & argv ) ;

      /* Free strings */
      mxFree( cmd ) ;
      mxFree( arg_s ) ;
    }
  else  /* no arg */
    {
      /* Translate no arg into mainfunc args */
      make_args( CMDNAME_S, & argc, & argv ) ;
    }

  /* Call function */
  printf( "=== Call mainfunc with %d args :\n===  " , argc ) ;
  for ( iarg = 0 ; iarg < argc ; iarg ++ )
    {
      printf( "[%s] ", argv[ iarg ] ) ;
    }
  printf( "\n" ) ;

  sts = mainfunc( argc, (const char* *) argv ) ;

  printf( "=== Returned from mainfunc (status = %d) \n", sts ) ;

  /* Free each element of array of args and the array of pointers itself */
/*  for ( iarg = 0 ; iarg < argc ; iarg ++ )
/*    {
/*	mxFree( argv[ iarg ] ) ;
/*    }
/*  mxFree( argv ) ;
 */

}



/* ------------------------------------------------------------------------ */
static void make_args
(
 const char* cmd,      /* I : command line to analyse */
 int*        argc_p,   /* O : number of arguments */
 char***     argv_p    /* O : array of argc arguments (allocated) */
)
/*\

    A command line CMD is taken as input.  Each distinct string within
    CMD is copied into a slot of (*ARGV_P).  After the call, the table
    (*ARGV_P) contain (*ARGC_P) allocated strings.

\*/
/* ------------------------------------------------------------------------ */
{
  char** tabstr ; /* table of n strings - allocate and copy in (*argv_p) */  

  /*
   *  1 : Count number of strings 
   */
  (*argc_p) = fetch_args( cmd, 0, NULL ) ;
  
  /* 
   *  2 : Allocate and fill array of these strings 
   */
  (*argv_p) = mxCalloc( (*argc_p), sizeof( char* ) ) ;
  if ( tabstr == NULL )
    mexErrMsgTxt( "Could not allocate array of args" ) ;

  fetch_args( cmd, (*argc_p), (*argv_p) ) ;
}


/* ------------------------------------------------------------------------ */
static int fetch_args  /* returns number of counted strings */
( 
 const char* cmd,  /* I : command line to analyze */
 int         nmax, /* I : number of allocated args in argv */
 char**      argv  /* O : table [0..nmax-1] strings (each string alloc here)*/
)
/*\ 

    A command line string CMD is taken as input.  If ARGV is not NULL,
    a string is allocated in each of the NMAX slots of ARGV, and the
    distinct strings of CMD are copied in those slots.  The function
    returns the number of distinct strings in CMD.

    Typically, first call fetch_args() to count in N the number of
    strings of CMD, then allocate a table of N pointers, then call
    again fetch_args() to copy the strings of CMD into the table.

\*/
/* ------------------------------------------------------------------------ */
{
  char*  copy_s ; /* copy of input string - allocated locally */
  char*  this_s ; /* pointer to current arg in copy_s */

  int    n ;      /* number of given strings - to copy in (*argc_p)*/

  copy_s = mxCalloc( strlen( cmd ) + 1 , sizeof( char ) ) ;
  if ( copy_s == NULL )
    mexErrMsgTxt( "Could not allocate copy of string of args" ) ;

  strncpy( copy_s, cmd, strlen( cmd ) + 1 ) ;
  n = 0 ;
  for ( this_s = strtok( copy_s , SEPAR_S ) ;
	this_s != NULL ;
	this_s = strtok( NULL , SEPAR_S ) )
    {
      if ( ( argv != NULL ) && ( n < nmax ) )
	{
	  argv[ n ] = mxCalloc( strlen( this_s ) + 1 , sizeof( char ) ) ;
	  if ( argv[ n ] == NULL )
	    mexErrMsgTxt( "Could not allocate an elt of array of args" ) ;

	  strncpy( argv[ n ] , this_s , strlen( this_s ) + 1 ) ;
	}
      
      n ++ ;
    }

  mxFree( copy_s ) ;

  return n ;
}

                                                                                                                                                                                                        ./._mexmemo.c                                                                                       000755  000765  000765  00000000312 11541743320 012714  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      mexmemo.c                                                                                           000755  000765  000024  00000002534 11541743320 012642  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\
    mexmemo.c

    Specific routines for MATLAB memory allocation

    Jan 1998

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    Vers-mod  Date         Description

    1.05-a    26-JAN-1998  Creation
\*/


#include "mex.h"    /* mxCalloc, mxFree */
#include <stdio.h>  /* stderr */



/* ------------------------------------------------------------------- */
void* GenAlloc
(
 size_t       nelem,        /* I : number of elements to allocate */ 
 size_t       elsize,       /* I : size in bytes of each element */
 int          doexit,       /* I : 1 if failure exits, 0 if only return NULL */
 const char*  where,        /* I : name of calling function */
 const char*  what          /* I : name of allocated object */
)
{
  void *result = mxCalloc( nelem , elsize );
  if ( result != NULL )
    {
      return result ;
    }
  else
  {
    fprintf(stderr, "Fatal: in %s, no memory for %s (%d elements size %d)\n",
	    where, what, nelem, elsize);
    if ( doexit )
      mexErrMsgTxt( " " );
    else
      return result ;
  }
}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
void GenFree( void* ptr )
{
  if ( ptr != NULL )
    mxFree( ptr ) ;
}
/* ------------------------------------------------------------------- */
                                                                                                                                                                    ./._genmemo.h                                                                                       000755  000765  000765  00000000312 11541743313 012703  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      genmemo.h                                                                                           000755  000765  000024  00000002103 11541743313 012621  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         #ifndef GENMEMO_H
#define GENMEMO_H

/*\
    genmemo.h

    Prototypes of generic routines for memory allocation

    Jan 1998

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    Vers-mod  Date         Description

    1.05-a    26-JAN-1998  Creation
    1.06-a    20-SEP-1998  Add macro freenull
\*/


#include <stdlib.h>  /* size_t */


/* ------------------------------------------------------------------- */
/* #define freenull( p ) do { \ 
   GenFree( p ); \
   p = NULL; \
   } while(0) 
*/


/* ------------------------------------------------------------------- */
void* GenAlloc
(
 size_t       nelem,        /* I : number of elements to allocate */ 
 size_t       elsize,       /* I : size in bytes of each element */
 int          doexit,       /* I : 1 if failure exits, 0 if only return NULL */
 const char*  where,        /* I : name of calling function */
 const char*  what          /* I : name of allocated object */
) ;

/* ------------------------------------------------------------------- */
void GenFree( void* ptr ) ;




#endif
                                                                                                                                                                                                                                                                                                                                                                                                                                                             ./._err2.c                                                                                          000755  000765  000765  00000000312 11541743312 012120  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      err2.c                                                                                              000755  000765  000024  00000003710 11541743312 012043  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*
 * program to compute permuted classification error in 2 cluster case.
 *
 * Mo Dang, oct 1997.
 */

#include <stdio.h>
#include <stdlib.h>

static int ReadIntVector( const char *fname, int* VP[], int nbmax ) ;


main(int argc, char *argv[])
{
  int    nbmax ;
  int*   clas1V ;
  int*   clas2V ;
  int    nbelts ;
  int    nbelts2 ;
  int    ielt ;
  int    nbdif ;
  float  percent ;


  if ( argc < 4 )
    {
      fprintf( stderr, "Usage : %s file1 file2 nbmax\n" , argv[0] ) ;
      fprintf( stderr, "  computes the percentage of misclassification\n" ) ;
      return 1 ;
    }

  nbmax = atoi( argv[ 3 ] ) ;
  nbelts = ReadIntVector( argv[1] , &clas1V , nbmax ) ;
  if ( nbelts <= 0 )
    return 2 ;

  nbelts2 = ReadIntVector( argv[2], &clas2V, nbmax ) ;
  if ( nbelts2 != nbelts )
    {
      fprintf( stderr, "Error : %d elts in %s and %d elts in %s\n", 
	       nbelts, argv[1], nbelts2, argv[2] ) ;
      return 3 ;
    }

  for ( ielt = 1, nbdif = 0 ; ielt <= nbelts ; ielt ++ )
    {
      if ( clas1V[ ielt ] != clas2V[ ielt ] )
	nbdif ++ ;
    }

  if ( nbdif > nbelts / 2 ) 
    nbdif = nbelts - nbdif ;

  percent = ( 100.0 * nbdif ) / nbelts ;

  fprintf( stdout , "Error = %3.1f %%  (%d on %d)\n", percent ,
	   nbdif , nbelts ) ;

  return 0 ;
}





static int ReadIntVector( const char *fname, int* VP[], int nbmax )
{
  FILE* fid ;
  int   nbelts ;

  if ( ( fid = fopen( fname , "r" ) ) == NULL )
    {
      fprintf( stderr , "Cannot open file %s\n" ,  fname ) ;
      *VP = NULL ;
      return 2 ;
    }

  if( ( *VP = calloc( nbmax + 1 , sizeof( int ) ) ) == NULL )
    {
      fprintf( stderr , "Not enough memory\n" ) ;
      fclose( fid ) ;
      return 2 ;
    }

  for ( nbelts = 0 ; ( nbelts <= nbmax ) && (! feof( fid ) ) ; )
    {
      int n ;

      if ( fscanf( fid , "%d", &n ) >= 1 )
	{
	  nbelts ++ ;
	  (*VP)[ nbelts ] = n ;
	}
    }

  if ( nbelts == 0 )
    {
      free( *VP ) ;
      *VP = NULL ;
    }

  fclose( fid ) ;

  return nbelts ;
}
                                                        ./._geo2nei.c                                                                                       000755  000765  000765  00000000312 11541743314 012600  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      geo2nei.c                                                                                           000755  000765  000024  00000055543 11541743314 012536  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\

    geo2nei.c

    This utility takes as input a set of geographic coordinates
    and computes the neighborhood system induced by the distances
    between the objects.

    November 1997

    Mo Van Dang
    Universite de Technologie de Compiegne
    URA CNRS 817

\*/

#include <stdio.h>      /* FILE , printf ... */
#include <stdlib.h>     /* malloc, ... */
#include <values.h>     /* MAXFLOAT */
#include <math.h>       /* sqrt */
#include <string.h>     /* memcpy, ...       */
#include "lib_io.h"     /* AskFileToRead ... */

#define VERSION   "1.00"


#define MAT( A , nl , nc , l , c ) ( ( A[ ( l * nc ) + c ] ) )


#define TRUE        1
#define FALSE       0

#define MAX_DIM     4
#define LEN_DESC    500

#define EXPFAC      3.0

#define SEPAR_STR   " \t"


typedef enum {
  WEI_NONE,
  WEI_CONST,
  WEI_EXP,
  WEI_NB
}
WeiModT ;


typedef struct
{
    int     Nelts;                /* number of elements */
    float   AveSumWei ;           /* ave sum of weights */
    float   MinSumWei ;           /* min sum of weights */
    float   MaxSumWei ;           /* max sum of weights */
    float   AveDisNearest ;       /* ave distance to 4th nearest neighbor */
    float   MinDisNearest ;       /* min distance to 4th nearest neighbor */
    float   MaxDisNearest ;       /* max distance to 4th nearest neighbor */
    float   AveNbNei ;            /* ave number of neighbors */
    int     MinNbNei ;            /* minimum number of neighbors */
    int     MaxNbNei ;            /* maximum number of neighbors */
}
StatsT ;



main( int argc , char *argv[] )
{
    int     sts = 0 ;                   /* status returned by the program */
    char    namedat[ LEN_FILE + 1 ] ;   /* input data file name */
    char    namenei[ LEN_FILE + 1 ] ;   /* output neighbors file name */
    char    descnei[ LEN_DESC + 1 ] ;   /* comment of neigbors file */
    int     nd ;                        /* number of spatial variables */
    int     d ;                         /* current spatial variable */
    int     spacol[ MAX_DIM ] ;         /* columns of spatial variables */
    int     mincol , ncol ;             /* input file : nb of columns */
    int     npt ;                       /* input file : nb of lines */
    float*  locM ;                      /* spatial locations : npt * nd */

    int     ranknei ;                   /* rank of nearest neighbor */
    float   weifac ;                    /* factor -> average sum weights = 4 */
    float   disthres ;                  /* pts neighbors if dis<=threshold */
    WeiModT weimod ;                    /* weight computation mode */

    StatsT  stats ;                     /* neighbor stats */

    void PrintHelp( const char* CmdName ) ;
    void PrintVersion( const char* CmdName ) ;

    int ComputeNeighbors
        (
         const float*   PtsM , 
         int            Npt , 
         int            Nc , 
	 int            RankNei ,
         float          Dth, 
         WeiModT        WeiMod , 
	 float          WeiFac ,
         const char*    NameNei , 
         const char*    Desc ,
         int            DoSave , 
	 StatsT*        NeiStatsP
        ) ;

    void DisplayStats
      (
       int              RankNei, /* I : rank of nearest neighbor */
       const StatsT*    StatsP   /* I : distances/neighbors statistics */
      ) ;


    if ( argc > 1 )
      {
	if ( strcmp( argv[ 1 ] , "-h" ) == 0 )
	  {
	    PrintHelp( argv[ 0 ] ) ;
	    return 1 ;
	  }
	else if ( strcmp( argv[ 1 ] , "-v" ) == 0 )
	  {
	    PrintVersion( argv[ 0 ] ) ;
	    return 1 ;
	  }
	else
	  {
	    printf( "Type %s -h to get help or -v to know current version\n" ,
		    argv[ 0 ] ) ;
	    return 2 ;
	  }
      }

    printf( "* * *  Welcome to GEO2NEI program V%s  * * *\n\n" ,
	    VERSION ) ;

    /* Ask for all parameters : 
       input file name, number and columns of spatial variables ;
    */
    if ( AskFileToRead( "input objects-variables", namedat ) != 0 )
       return 2 ;
    if ( AskFileToWrite( "output neighborhood" , FALSE , namenei ) != 0 )
       return 2 ;
    printf( "Enter comment in neighbors file : " ) ;
    gets( descnei ) ;

    if ( AskInteger( "number of spatial coordinates" , 2 , 1 , MAX_DIM , 
                     & nd ) != 0 )
       return 2 ;

    for ( d = 0 ; d < nd ; d ++ )
    {
        char    msg[ 120 ] ;

        sprintf( msg , "column of spatial coordinate %d" , d + 1 ) ;
        if ( AskInteger( msg , d + 1 , 1 , 100 , & spacol[ d ] ) != 0 )
           return 2 ;
        spacol[ d ] -- ; /* C-like indices start from 0 */
    }

    if ( AskInteger( "Average number of neighbors" , 4 , 1 , 20 , 
                     & ranknei ) != 0 )
      return 2 ;

    if ( AskInteger( "Weight computation {0=>none, 1=>const 2=>exp[(-d/dmax)^2]}" , 
		     WEI_NONE , WEI_NONE , WEI_NB-1 , (int*) & weimod ) != 0 )
      return 2 ;
    
/*    if ( AskFloat( "Distance threshold" , 1.00 , 0.01 , 100.00 , 
		     & disthres ) != 0 )
	 return 2 ;
 */

    printf( "Analyzing file %s ...\n" , namedat ) ;
    switch( CountLinesColumns( namedat , SEPAR_STR , 
                               &mincol , &ncol , &npt ) )
    {
        case 0 : printf( "File %s : %d points, %d columns\n", 
                         namedat, npt, ncol ) ;
                 break ;
        case 1 : printf( "Error : '%s' unequal columns (%d to %d)\n" ,
                         namedat , mincol , ncol ) ;
                 return 2 ;
        default :
                printf( "Error while analyzing file %s\n" , namedat ) ;
                return 2 ;
    }

    for ( d = 0 ; d < nd ; d ++ )
    {
        if ( spacol[ d ] >= ncol )
        {
            printf( "Error : spatial var. %d's column = %d > nb columns\n" ,
                    d + 1 , spacol[ d ] + 1 ) ;
            return 2 ;
        }
    }


    /* Allocate and read spatial locations from file to memory 
    */
    if ( ( locM = malloc( npt * nd * sizeof( float ) ) ) ==  NULL )
    {
        printf( "Error : out of memory for %d site locations\n" , npt ) ; 
        return 2 ;
    }

    if ( ReadSelectedColumns( namedat , npt , ncol , nd , spacol , 
                              locM ) != 0 )
    {
        printf( "Error while reading spatial locations from file\n" ) ;
        free( locM ) ;
        return 2 ;
    }

    /* Compute and save neighborhood system
    */
    printf( "Computing neighborhood graph ...\n" ) ;

    /* Step 1 : compute distance stats to set distance threshold */
    disthres = 1.0 ;
    weifac = 1.0 ;
    if ( ComputeNeighbors( locM , npt , nd , ranknei , disthres , 
			   WEI_NONE , weifac , namenei , descnei , FALSE , 
			   &stats ) != -1 )
      {
	/* Step 2 : compute neighborhood graph with non-normalized weights */
	disthres = stats.AveDisNearest * 1.1 ;
	if ( ComputeNeighbors( locM , npt , nd , ranknei , disthres , 
			       weimod , weifac , namenei , descnei , FALSE , 
			       &stats ) != -1 )	
	  {
	    /* Step 3 : save neighborhood graph with normalized weights */
	    weifac = 4.0 / stats.MaxSumWei ;
	    printf( "weifac = %g\n", weifac ) ;
	    if ( ComputeNeighbors( locM , npt , nd , ranknei , disthres , 
			       weimod , weifac , namenei , descnei , TRUE , 
			       &stats ) != -1 )	
	      DisplayStats( ranknei, &stats ) ;
	    else
	      printf( "Error while computing neighborhood graph 2\n" ) ;
	  }
	else
	  printf( "Error while computing neighborhood graph 1\n" ) ;
	    
      }
    else
      printf( "Error while computing distance statistics\n" ) ;


    /* Free allocated memory */
    free( locM ) ;

    return sts ;

} /* main() */




/* ------------------------------------------------------------------- */
int ComputeNeighbors
        (
         const float*   PtsM , 
         int            Npt , 
         int            Nd , 
	 int            RankNei ,
         float          Dth , 
         WeiModT        WeiMod , 
	 float          WeiFac ,
         const char*    NameNei , 
         const char*    Desc ,
         int            DoSave , 
	 StatsT*        StatsP
        ) 
/* ------------------------------------------------------------------- */
{
    FILE*  fnei ;       /* File to write into */
    float  sqrDth ;     /* Squared distance threshold */

    float* sqrDisV ;    /* [Npt] Current point's distance to other points */
    float* weiV ;       /* [Npt] Current point's neighbors weights > 0.0 */
    int*   neiV ;       /* [Npt] Current point's neighbors indices 0..Npt-1 */

    int    pt ;         /* Index of current point 0..Npt-1 */
    int    ok ;         /* FALSE if problem in writing file */

    void StartStats
      (
       int      Npt ,          /* I : number of points */
       StatsT*  StatsP         /* O : distances/neighbors statistics */
      ) ;

    void CompuSqrDisV
      (
       int           Pt ,      /* I : index of current point 0..Npt-1 */
       int           RankNei,  /* I : neighbor rank to compute distance */
       const float*  PtsM ,    /* I : [Npt*Nd] matrix of point coordinates */
       int           Npt ,     /* I : number of points */
       int           Nd ,      /* I : number of coordinates */
       float*        SqrDisV , /* O : [Npt] distances to other points */
       float*        SqrDnearP /* O : distance to nearest point */
      ) ;

    void SqrDisVToNeiV
      ( 
       int           Pt ,      /* I : index of current point 0..Npt-1 */
       float         SqrDth ,  /* I : squared distance threshold */ 
       const float*  SqrDisV , /* I : [Npt] distances to other points */
       int           Npt ,     /* I : number of points */
       WeiModT       WeiMod ,  /* I : weight computation mode */
       float         WeiFac ,  /* I : weight normalizing factor */
       int*          NbNeiP ,  /* O : number of neighbors 0..Npt-1 */
       int*          NeiV ,    /* O : [Npt] neighbors indices 0..Npt-1 */
       float*        WeiV ,    /* O : [Npt] neighbors weights > 0.0 */
       float*        SumWeiP   /* O : sum of weights > 0.0 */
      ) ;

    void UpdateStats
      (
       int           Pt ,      /* I : index of current point 0..Npt-1 */
       float         SqrDnear, /* I : distance to nearest point */
       int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
       float         SumWei ,  /* I : sum of neighbor weights */
       StatsT*       StatsP    /* I/O : distances/neighbors statistics */
      ) ;

    int SaveNeiV
      ( 
       FILE*         Fnei ,    /* I/O : file to write into */
       int           Pt ,      /* I : index of current point 0..Npt-1 */
       int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
       const int*    NeiV      /* I : [Npt] neighbors indices 0..Npt-1 */
      ) ;

    int SaveWeiV
      ( 
       FILE*         Fnei ,    /* I/O : file to write into */
       int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
       const float*  WeiV      /* I : [Npt] neighbors weights > 0.0 */
      ) ;

    void CloseStats
      (
       StatsT*  StatsP         /* I/O : distances/neighbors statistics */
      ) ;



    /* Convert parameters to squared distances */
    sqrDth = Dth * Dth ;

    /* If save to file requested, create file and write header */
    if ( DoSave ) 
      {
	if ( ( fnei = fopen( NameNei , "w" ) ) == NULL )
	  {
	    printf( "Error : cannot create file %s\n", NameNei ) ;
	    return -1 ;
	  }

	fprintf( fnei , "# %s\n" , Desc ) ;
	if ( WeiMod == WEI_NONE )
	  fprintf( fnei , "0\n" ) ;
	else
	  fprintf( fnei , "1\n" ) ;
      }

    /* Allocate memory for local computations */
    sqrDisV = malloc( Npt * sizeof( float ) ) ;
    weiV    = malloc( Npt * sizeof( float ) ) ;
    neiV    = malloc( Npt * sizeof( int ) ) ;
    if ( ( sqrDisV == NULL ) || ( weiV == NULL ) || ( neiV == NULL ) )
      {
	printf( "Could not allocate distance/weight/neigbor vectors [%d]\n" ,
		 Npt ) ;
	if ( sqrDisV != NULL ) free( sqrDisV ) ;
	if ( weiV    != NULL ) free( weiV ) ;

	return -1 ;
      }


    /* For each point */
    for ( pt = 0 , 
	    ok = TRUE , 
	    StartStats( Npt , StatsP ) ; 
	  ( pt < Npt ) && ok ; 
	  pt ++ )
      {
	float   sqrdisNearest ;   /* Distance to nearest point */
	int     nbnei ;           /* Number of neighbors */
	float   sumwei ;          /* Sum of neighbor weights */

	/* Compute its distance to other points */
	CompuSqrDisV( pt , RankNei, PtsM , Npt , Nd , 
		      sqrDisV , &sqrdisNearest ) ;

	/* Threshold the computed distances to get neighbors and weights */
	SqrDisVToNeiV( pt , sqrDth , sqrDisV , Npt , WeiMod , WeiFac ,
		       &nbnei , neiV , weiV , &sumwei ) ;

	/* printf( "sumwei %d = %g\n" , pt , sumwei ) ; */

	/* Update distance/neighbor statistics */
	UpdateStats( pt , sqrdisNearest , nbnei , sumwei , StatsP ) ;

	/* Eventually save neighbors and weights to file */
	if ( DoSave )
	  {
	    if ( SaveNeiV( fnei , pt , nbnei , neiV ) != 0 )
	      ok = FALSE ;
	    else
	      if ( WeiMod != WEI_NONE )
		if ( SaveWeiV( fnei , nbnei , weiV ) != 0 )
		  ok = FALSE ;
	    fprintf( fnei , "\n" ) ;
	  }
      }
    CloseStats( StatsP ) ;


    /* Free memory used for local computations */
    free( sqrDisV ) ;
    free( neiV ) ;
    free( weiV ) ;

    /* If save to file requested, close file */
    if ( DoSave ) 
      {
	fclose( fnei ) ;
      }

    if ( ok ) 
      return 0 ;
    else
      return -1 ;

} /* end of ComputeNeighbors() */



/* ------------------------------------------------------------------- */
void StartStats
 (
  int      Npt ,          /* I : number of points */
  StatsT*  StatsP         /* O : distances/neighbors statistics */
 ) 
/* ------------------------------------------------------------------- */
{
    StatsP->Nelts = 0 ;

    StatsP->AveSumWei = 0.0 ;
    StatsP->MinSumWei = MAXFLOAT ;
    StatsP->MaxSumWei = 0.0 ;

    StatsP->AveDisNearest = 0.0 ;
    StatsP->MinDisNearest = MAXFLOAT ;
    StatsP->MaxDisNearest = MINFLOAT ;

    StatsP->AveNbNei     = 0.0 ;
    StatsP->MinNbNei      = Npt ;
    StatsP->MaxNbNei      = 0 ;
}


/* ------------------------------------------------------------------- */
void UpdateStats
 (
  int           Pt ,      /* I : index of current point 0..Npt-1 */
  float         SqrDnear, /* I : distance to nearest point */
  int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
  float         SumWei ,  /* I : sum of neighbor weights */
  StatsT*       StatsP    /* I/O : distances/neighbors statistics */
 ) 
/* ------------------------------------------------------------------- */
{

    StatsP->Nelts ++ ;


    StatsP->AveSumWei += SumWei ;

    if ( SumWei < StatsP->MinSumWei )
      StatsP->MinSumWei = SumWei ;

    if ( SumWei > StatsP->MaxSumWei )
      StatsP->MaxSumWei = SumWei ;


    StatsP->AveDisNearest += sqrt( SqrDnear ) ;

    if ( SqrDnear < StatsP->MinDisNearest )
      StatsP->MinDisNearest = SqrDnear ;

    if ( SqrDnear > StatsP->MaxDisNearest )
      StatsP->MaxDisNearest = SqrDnear ;


    StatsP->AveNbNei += NbNei ;

    if ( NbNei < StatsP->MinNbNei ) 
      StatsP->MinNbNei = NbNei ;

    if ( NbNei > StatsP->MaxNbNei ) 
      StatsP->MaxNbNei = NbNei ;

}


/* ------------------------------------------------------------------- */
void CloseStats
 (
  StatsT*  StatsP         /* I/O : distances/neighbors statistics */
 ) 
/* ------------------------------------------------------------------- */
{
  if ( StatsP->MinDisNearest >= 0 )
    StatsP->MinDisNearest = sqrt( StatsP->MinDisNearest ) ;

  if ( StatsP->MaxDisNearest >= 0 )
    StatsP->MaxDisNearest = sqrt( StatsP->MaxDisNearest ) ;

  if ( StatsP->Nelts > 0 )
    {
      StatsP->AveSumWei     /= StatsP->Nelts ;
      StatsP->AveDisNearest /= StatsP->Nelts ;
      StatsP->AveNbNei      /= StatsP->Nelts ;
    }
}


/* ------------------------------------------------------------------- */
void DisplayStats
 (
  int              RankNei, /* I : rank of nearest neighbor */
  const StatsT*    StatsP   /* I : distances/neighbors statistics */
 ) 
/* ------------------------------------------------------------------- */
{
    printf( "Statistics (ave, min, max) : \n" ) ;
    printf( "  Distance to %dth nearest site :  %8.2f  (%8.2f to %8.2f)\n" ,
	    RankNei,
	    StatsP->AveDisNearest, 
	    StatsP->MinDisNearest , StatsP->MaxDisNearest );
    printf( "  Number of neighbors          :  %8.2f  (%8d to %8d)\n" , 
	    StatsP->AveNbNei ,
	    StatsP->MinNbNei , StatsP->MaxNbNei );
    printf( "  Sum of neighbor weights      :  %8.2f  (%8.2f to %8.2f)\n" ,
	    StatsP->AveSumWei, 
	    StatsP->MinSumWei , StatsP->MaxSumWei );
}

/* ------------------------------------------------------------------- */
static int floatcompare(const void *x, const void *y)
/* ------------------------------------------------------------------- */
{
  const float* xx = x;
  const float* yy = y;

  if (*xx > *yy)
    return (1);
  if (*xx < *yy)
    return (-1);
  return (0);
}

/* ------------------------------------------------------------------- */
/* Computes the distances from a given point to other points 
 */
void CompuSqrDisV
 (
  int           Pt ,      /* I : index of current point 0..Npt-1 */
  int           RankNei,  /* I : neighbor rank to compute distance */
  const float*  PtsM ,    /* I : [Npt*Nd] matrix of point coordinates */
  int           Npt ,     /* I : number of points */
  int           Nd ,      /* I : number of coordinates */
  float*        SqrDisV , /* O : [Npt] distances to other points */
  float*        SqrDnearP /* O : distance to RankNei nearest point */
 )
/* ------------------------------------------------------------------- */
{
    int   j ;
    float dnear ;

    for ( j = 0 , dnear = MAXFLOAT ; j < Npt ; j ++ )
      {
	if ( j != Pt )
	  {
	    int   d ;
	    float sqrdis ;

	    for ( d = 0 , sqrdis = 0.0 ; d < Nd ; d ++ )
	      {
	        float cpt = MAT( PtsM , Npt , Nd , Pt , d ) ;
		float cj  = MAT( PtsM , Npt , Nd , j , d ) ;
		float dif = cpt - cj ;

		sqrdis = sqrdis + dif * dif ;
	      }

	    /* +++ */ if ( sqrdis == 0.0 ) 
	      printf( "*** Warning : Points %d and %d have same location\n" ,
		      Pt+1 , j+1 ) ;

	    SqrDisV[ j ] = sqrdis ;
	    if ( sqrdis < dnear ) 
	      dnear = sqrdis ;
	  }
	else
	  SqrDisV[ j ] = 0.0 ;
      }

    /* Compute distance to RankNei'th nearest location */
    {
      float* sortdis_1n = calloc( Npt , sizeof( float ) ) ;

      if ( sortdis_1n == NULL )
	{
	  (*SqrDnearP) = dnear ;
	  return ;
	}

      /* Sort the distances to other locations 
       */
      memcpy( sortdis_1n , SqrDisV , Npt * sizeof( float ) ) ;

      qsort( sortdis_1n , Npt, sizeof( float ), floatcompare ) ;

      (*SqrDnearP) = sortdis_1n[ RankNei - 1 ] ;

      free( sortdis_1n ) ;
    }
}


/* ------------------------------------------------------------------- */
/* Thresholds the computed distances to get neighbors and weights */
void SqrDisVToNeiV
 ( 
  int           Pt ,      /* I : index of current point 0..Npt-1 */
  float         SqrDth ,  /* I : squared distance threshold */ 
  const float*  SqrDisV , /* I : [Npt] distances to other points */
  int           Npt ,     /* I : number of points */
  WeiModT       WeiMod ,  /* I : weight computation mode */
  float         WeiFac ,  /* I : weight normalizing factor */
  int*          NbNeiP ,  /* O : number of neighbors 0..Npt-1 */
  int*          NeiV ,    /* O : [Npt] neighbors indices 0..Npt-1 */
  float*        WeiV ,    /* O : [Npt] neighbors weights > 0.0 */
  float*        SumWeiP   /* O : sum of weights > 0.0 */
 ) 
/* ------------------------------------------------------------------- */
{
    int   j ;

    for ( j = 0 , 
	    (*NbNeiP) = 0,
	    (*SumWeiP) = 0.0 ; 
	  j < Npt ; 
	  j ++ )
      {
	if ( j != Pt )
	  {
	    if ( SqrDisV[ j ] <= SqrDth )
	      {
		NeiV[ (*NbNeiP) ] = j ;

		switch( WeiMod )
		  {
		  case WEI_CONST:
		    WeiV[ (*NbNeiP) ] = 1.0 ;
		    break ;

		  case WEI_EXP:
		    WeiV[ (*NbNeiP) ] = exp( - EXPFAC * 
					     SqrDisV[ j ] / SqrDth ) ;
		    break ;

		  default:
		    WeiV[ (*NbNeiP) ] = 1.0 ;
		  }
		WeiV[ (*NbNeiP) ] *= WeiFac ;

		(*SumWeiP) += WeiV[ (*NbNeiP) ] ;

		(*NbNeiP) ++ ;
	      }
	  }
      }
}



/* ------------------------------------------------------------------- */
int SaveNeiV
 ( 
  FILE*         Fnei ,    /* I/O : file to write into */
  int           Pt ,      /* I : index of current point 0..Npt-1 */
  int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
  const int*    NeiV      /* I : [Npt] neighbors indices 0..Npt-1 */
 ) 
/* ------------------------------------------------------------------- */
{
    int nei ;
    int ok ;

    fprintf( Fnei , "%4d  %3d  " , Pt + 1 , NbNei ) ;
    for ( nei = 0 , ok = TRUE ; ( nei < NbNei ) && ok ; nei ++ )
      {
	ok = ( fprintf( Fnei , "%4d " , NeiV[ nei ] + 1 ) > 0 ) ;
      }

    return ( ok ? 0 : -1 ) ;
}



int SaveWeiV
 ( 
  FILE*         Fnei ,    /* I/O : file to write into */
  int           NbNei ,   /* I : number of neighbors 0..Npt-1 */
  const float*  WeiV      /* I : [Npt] neighbors weights > 0.0 */
 ) 
{
    int nei ;
    int ok ;

    for ( nei = 0 , ok = TRUE ; ( nei < NbNei ) && ok ; nei ++ )
      {
	ok = ( fprintf( Fnei , " %4.2f" , WeiV[ nei ] ) > 0 ) ;
      }

    return ( ok ? 0 : -1 ) ;
}



void PrintHelp( const char* CmdName )
{
  printf( "\n" ) ;
  printf( "This program computes a neighborhood graph given a set of \n" ) ;
  printf( "spatial coordinates. It saves the resulting graph in a \n" ) ;
  printf( "neighborhood file which may be used as input to the \n" ) ;
  printf( "program nem_exe. It uses a simple thresholding of the \n" ) ;
  printf( "euclidean distances between the objects.\n" ) ;
  printf( "\n" ) ;
  printf( "Before running %s, the spatial coordinates should be given in \n" ,
	  CmdName ) ;
  printf( "an ASCII file, 1 line/object, 1 column/spatial coordinate.\n" ) ;
  printf( "Example (spatial coordinates in columns 2 and 3 of file) : \n" ) ;
  printf( "   x 15 50 x x   => 1st object is at position x = 15  y = 50\n" ) ;
  printf( "   x 30 20 x x   => 2nd object is at position x = 30  y = 20\n" ) ;
  printf( "The elements on a line should be separated by spaces or tabs\n" ) ;
  printf( "\n" ) ;
  printf( "--- Press ENTER for more ---\n" ) ;
  getchar( ) ;
  printf( "The program will prompt you for :\n" ) ;
  printf( "1 - The name of the input spatial coordinates file\n" ) ;
  printf( "2 - The name of the output neighborhood file\n" ) ;
  printf( "3 - How many spatial coordinates and their column numbers\n" ) ;
  printf( "4 - The desired average number of neighbors of an object ;\n" ) ;
  printf( "    this parameter is used to compute the distance threshold\n" ) ;
  printf( "5 - How to calculate the weights : no weigths or exponential\n" ) ;
  printf( "    exponential means  w_ij = A * exp - 3 (d_ij / d_th)^2\n" ) ;
  printf( "    i.e. weights decreasing as  exp (- squared distance).\n" ) ;
  printf( "    Constant A is computed to have max_i( sum_j w_ij ) = 4,\n" ) ;
  printf( "    as in 4 nearest neighbor unweighted graphs\n" ) ;
  printf( "\n" ) ;
}


void PrintVersion( const char* CmdName )
{
  printf( "\n" ) ;
  printf( "Version 1.00 (14-NOV-1997)\n" ) ;
  printf( "==========================\n" ) ;
  printf( "14-NOV-1997\n" ) ;
  printf( "First complete version. \n" ) ;
  printf( "Added help, computation of distance threshold and weights.\n" ) ;
  printf( "\n" ) ;
}
                                                                                                                                                             ./._randord.c                                                                                       000755  000765  000765  00000000312 11541743336 012705  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      randord.c                                                                                           000755  000765  000024  00000005215 11541743336 012632  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\

    randord.c

    This utility outputs a random permutation of set {1,...,N},
    where N is given on the command line. The result is printed
    on standard output.

    05-NOV-1997

    Mo Van Dang
    UMR CNRS 6599
    Universite de Technologie de Compiegne, France

\*/


#include <stdio.h>     /* printf() */
#include <stdlib.h>    /* srandom() */

#include "nem_rnd.h"   /* RandomSeedByTime(), RandomInteger() */


static void PrintUsage( const char* CmdS ) ;
static void Error( const char *MsgS , int ExitCode ) ;
static void GiveRandomSeed( int  MySeed ) ;


/* =========================== */

int main( int argc, char *argv[] )

/* =========================== */
{
  int  nbElts ;
  int* tabEltsV ;
  int  i ;


  /* Check arguments */
  if ( argc < 2 )
    {
      PrintUsage( argv[0] ) ;
      return 1 ;
    }

  nbElts = atoi( argv[ 1 ] ) ;
  if ( nbElts <= 0 )
    Error( "Number of elements must be greater than zero" , 1 ) ;


  if ( argc < 3 )
    RandomSeedByTime( ) ;
  else
    GiveRandomSeed( atoi( argv[ 2 ] ) ) ;


  /* Allocate and intialize the vector of integers */
  if ( ( tabEltsV = calloc( nbElts , sizeof( int ) ) ) == NULL )
    Error( "Not enough memory for that many elements" , 2 ) ;

  for ( i = 0 ; i < nbElts ; i ++ )
    tabEltsV[ i ] = i + 1 ;

  /* Run random permutation algorithm */
  RandomPermutationAlgo( tabEltsV , nbElts ) ;


  /* Print result to standard output */
  for ( i = 0 ; i < nbElts ; i ++ )
    fprintf( stdout , "%3d " , tabEltsV[ i ] ) ;
  fprintf( stdout , "\n" ) ;

  free( tabEltsV ) ;
  return 0 ;
}


/* =========================== */

static void PrintUsage( const char* CmdS )

/* =========================== */
{
  fprintf( stderr , "\nSyntax :\n\n" ) ;
  fprintf( stderr , "  %s  N  [ seed ]\n\n" , CmdS ) ;
  fprintf( stderr , "  This utility prints to standard output a random permutation\n" ) ;
  fprintf( stderr , "  of set {1,...,N}.  A seed of the random generator may be specified\n" ) ;
  fprintf( stderr , "  (by default, seed = system time).\n\n" ) ;
}


/* =========================== */

static void Error( const char *MsgS , int ExitCode )

/* =========================== */
{

  switch( ExitCode )
    {
    case 1 : 
      fprintf( stderr , "\n*** Input Error : %s\n\n" , MsgS ) ;
      break ;

    case 2 :
      fprintf( stderr , "\n*** Runtime Error : %s\n\n" , MsgS ) ;
      break ;

    default :
      fprintf( stderr , "\n*** Fatal Error : %s\n\n" , MsgS ) ;
    }

  exit( ExitCode ) ;

}



/* =========================== */

static void GiveRandomSeed( int  MySeed )

/* =========================== */
{
#ifdef __TURBOC__
    srand( (unsigned) MySeed ) ;
#else
    srandom( MySeed ) ;
#endif
}



                                                                                                                                                                                                                                                                                                                                                                                   ./._txt2hlp.c                                                                                       000755  000765  000765  00000000312 11541743342 012656  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      txt2hlp.c                                                                                           000755  000765  000024  00000013064 11541743342 012604  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*\
    txt2hlp.c

    This utility converts a text file into a C file consisting of a
    'help' function. "Batch" version (parameters given by command line
    arguments, no keyboard input requested).

    New version 2.0 11-JAN-1998 : allows several help functions.
\*/

#include <stdio.h>      /* printf, FILE */
#include <string.h>     /* strcpy */
#include <time.h>       /* time */


#define TRUE        1
#define FALSE       0

#define LEN_FILE    100
#define LEN_FUNC    32 
#define LEN_LINE    100
#define LEN_COMMENT 500

#define NB_PARA     4

#define MARKER_CHAR '%'

main( int Argc , char* Argv[] )
{
    char    ftextName[ LEN_FILE + 1 ] ;
    char    fhelpName[ LEN_FILE + 1 ] ;
    char    funcName[ LEN_FUNC + 1 ] ;
    char    CommentS[ LEN_COMMENT + 1 ] ;
    char    rep ;
    FILE*   fh ;

    int     hFileYes ;

    FILE*   ftext ;

    void PrintUsage( const char* Cmd ) ;


    printf( "\nWelcome to program TXT2HLP\n\n" ) ;


    /* Fetch command line arguments */
    if ( ( Argc - 1 ) == 0 )
      {
	/* No args -> print help */
	PrintUsage( Argv[ 0 ] ) ;
	return 1 ;
      }

    if ( ( Argc - 1 ) < NB_PARA )
      {
	printf( "Error : only %d parameters (%d requested)\n" ,
	         ( Argc - 1 ) ,  NB_PARA ) ;
	PrintUsage( Argv[ 0 ] ) ;
	return 2 ;
      }


    strncpy( ftextName , Argv[ 1 ] , LEN_FILE ) ;
    strncpy( fhelpName , Argv[ 2 ] , LEN_FILE ) ;
    rep                = Argv[ 3 ][ 0 ] ;
    strncpy( CommentS  , Argv[ 4 ] , LEN_COMMENT ) ;

    printf( "Input text file name  : %s\n" , ftextName ) ;
    printf( "Output help file name : %s.c\n" , fhelpName ) ;
    printf( "___.h file ? (y/n)    : %c\n" , rep ) ;
    hFileYes = ( rep == 'y' ) ;
    printf( "Comment               : %s\n" , CommentS ) ;
    printf( "\n" );


    /* Open input file */
    if ( ( ftext = fopen( ftextName , "r" ) ) == NULL )
    {
        printf( "Error : cannot read file %s\n" , ftextName ) ;
        return 2 ;
    }


    /* If ___.h file requested, create it */
    if ( hFileYes )
    {
        char    fhName[ LEN_FILE + 1 ] ;

        strcpy( fhName , fhelpName ) ;
        strncat( fhName ,  ".h" , LEN_FILE ) ;

        /* Open ___.h file */
        if ( ( fh = fopen( fhName , "w" ) ) == NULL )
        {
            printf( "Error : cannot write file %s\n" , fhName ) ;
            fclose( ftext ) ;
            return 2 ;
        }

        printf( "Writing file %s ...\n" , fhName ) ;

	/* Start conditional inclusion of .h file */
        fprintf( fh , "#ifndef %s_H\n", fhelpName ) ;
	fprintf( fh , "#define %s_H\n\n", fhelpName ) ;

        /* Include necessary ___.h files */
        fprintf( fh , "#include <stdio.h>  /* FILE */\n\n" ) ;

    }


    /* Treat C file */
    {
        char    fcName[ LEN_FILE + 1 ] ;
        FILE*   fc ;
	time_t  timer = time( NULL ) ;
	int     in_function ;  /* 1 if in a function, 0 else */

        strcpy( fcName , fhelpName ) ;
        strncat( fcName ,  ".c" , LEN_FILE ) ;

        /* Open ___.c file */
        if ( ( fc = fopen( fcName , "w" ) ) == NULL )
        {
            printf( "Error : cannot write file %s\n" , fcName ) ;
            fclose( ftext ) ;
            return 2 ;
        }

        printf( "Writing file %s ...\n" , fcName ) ;

        /* Start ___.c file */
        fprintf( fc , "/*\\\n\n" ) ;

        fprintf( fc , "    %s.c\n\n" , fhelpName ) ;

        fprintf( fc , "    %s\n\n" , CommentS ) ;

	fprintf( fc , "    %s\n", asctime( localtime( &timer ) ) ) ;

        fprintf( fc , "\\*/\n\n" ) ; 

        fprintf( fc , "#include \"%s.h\"  /* Exported prototypes */\n\n" , 
                 fhelpName ) ;

        /* Scan lines from text file to ___.c file */
	in_function = 0 ;
	funcName[ 0 ] = '\0' ;
        while ( ! feof( ftext ) )
        {
            char    line[ LEN_LINE + 1 ] ;

            if ( fgets( line , LEN_LINE , ftext ) != NULL )
	      {
		/* Suppress newline character */
                int     len     = strlen( line ) ;
                line[ len - 1 ] = '\0' ;

		/* If this line is a structure marker */
		if ( line[ 0 ] == MARKER_CHAR )
		  {
		    /* If not currently in function, start a new function 
		       and include it in header file */
		    if ( ! in_function )
		      {
			in_function = 1 ;
			strncpy( funcName, &line[ 1 ], LEN_FUNC ) ;
			fprintf( fc , "void %s( FILE* F )\n{\n" , funcName ) ;
			fprintf( fh , "extern void %s( FILE* F ) ;\n\n" , 
				 funcName ) ;
		      }
		    else /* currently in a function, end it */
		      {
			in_function = 0 ;
			fprintf( fc , "} /* end of %s() */\n\n\n", funcName ) ;
		      }
		  }
		else /* not a structure marker, if in function, 
			make print out instruction */
		  {
		    if ( in_function )
		      fprintf( fc , "    fprintf( F , \"%s\\n\" ) ;\n" , 
			       line ) ;
		  }

	      } /* end if fgets() != NULL */
        }

	/* If current function not ended, end it */
	if ( in_function )
	  {
	    printf( "Warning : function %s not ended, I will do it\n", 
		    funcName ) ;
	    fprintf( fc , "} /* end of %s() */\n\n\n", funcName ) ;
	  }

        /* Close ___.c file */
        fclose( fc ) ;
    }


    /* If ___.h file requested, close it */
    if ( hFileYes )
      {
	/* Close conditional inclusion of .h file */
        fprintf( fh , "#endif\n" ) ;

        /* Close ___.h file */
        fclose( fh ) ;
      }

    /* Close input file */
    fclose( ftext ) ;

    printf( "\nBye\n" ) ;

    return 0 ;
} /* end of main() */



/* -------------------------------------------------------------- */
void PrintUsage( const char* Cmd )
{
    printf( "\nSyntax :\n" ) ;
    printf( "   %s   in_file  out_file  hfile_yn  comment\n\n" , 
	    Cmd ) ;
}



                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ./._tcpu.c                                                                                          000755  000765  000765  00000000312 11541743342 012224  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      tcpu.c                                                                                              000755  000765  000024  00000001533 11541743342 012150  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         /*
 *  tcpu.c
 *
 *  This utility computes the cpu time used by a given command
 *  and writes this time into a specified file.
 */

#include <stdio.h>      /* FILE */
#include <stdlib.h>     /* malloc() */
#include <sys/times.h>  /* times() */
#include <time.h>  /* CLK_TCK() */

int main(int argc, char *argv[] ) 
{
   clock_t begintime, endtime;
   struct tms *a_tms;

   FILE* fid ;
   int   sts ;


   if ( argc < 3 )
     {
       fprintf( stderr , "Usage : %s  cmd  timefile\n" , argv[ 0 ] ) ;
       return 1 ;
     }

   a_tms = ( struct tms *) malloc( sizeof (struct tms));

   times(a_tms); begintime = a_tms->tms_cutime;

   sts = system( argv[1] );

   times(a_tms); endtime = a_tms->tms_cutime;

   fid = fopen( argv[2], "w" ) ;
   fprintf( fid , " %5.3f \n", ((double)(endtime-begintime)/CLOCKS_PER_SEC));
   fclose( fid ) ;

   return sts ;
}
                                                                                                                                                                     ./._mainfunc.h                                                                                      000755  000765  000765  00000000312 11541743317 013060  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      mainfunc.h                                                                                          000755  000765  000024  00000000310 11541743317 012774  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         #ifndef MAINFUNC_H
#define MAINFUNC_H

/*
 * mainfunc.h
 *
 * Prototype of main function called either by mexFunction() or main()
 *
 */

extern int mainfunc( int argc, const char** argv ) ; 

#endif
                                                                                                                                                                                                                                                                                                                        ./._simimg.str                                                                                      000755  000765  000765  00000000312 11541743342 013124  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      simimg.str                                                                                          000755  000765  000024  00000000150 11541743342 013042  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         # File created by nemcdon.m on 05/11/1997 at 02h23
# Simulated image beta=1.5, k=2, sigma=1
I 100 100 1
                                                                                                                                                                                                                                                                                                                                                                                                                        ./._simimg.dat                                                                                      000755  000765  000765  00000000312 11541743341 013063  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      simimg.dat                                                                                          000755  000765  000024  00000257620 11541743341 013021  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                           2.1839
  1.5378
  0.5381
  1.1025
  2.1614
  0.8342
  2.0234
 -0.1198
 -0.9461
  0.7311
 -0.4083
 -0.4864
 -1.9502
 -0.8633
 -1.0229
 -2.3825
  0.9896
  0.8265
 -0.3908
 -2.0716
 -0.7092
 -0.0032
  0.7295
  1.8306
  1.2264
  1.0246
  2.0928
  2.1985
  1.8558
  1.0050
  0.3243
  0.0160
  1.6184
 -0.1057
  1.3768
  1.1969
  0.8186
  1.1499
  1.8195
  1.2787
 -0.1573
  0.1681
  0.4787
 -0.3984
 -1.1890
 -0.9077
  1.7228
  0.9916
 -1.4672
 -0.8577
  0.8030
  0.6713
  0.0686
  0.2751
  0.0830
  0.4669
  0.6556
  3.9432
  0.1200
 -0.8769
  0.9851
  0.0232
  1.1139
  0.5892
  0.4235
 -0.0790
 -0.1927
  0.2345
  0.1185
  1.2854
 -1.6853
 -1.4921
 -0.6032
 -0.7846
  0.4844
  0.5923
 -0.4456
 -0.6237
  1.5588
 -0.6800
 -0.0530
  2.6234
  0.5046
  1.4754
  0.6137
  1.5201
 -0.2498
 -1.3894
 -0.6972
  0.7407
 -0.9530
  0.2355
  1.2474
 -1.6984
  0.4735
 -1.0464
  0.9288
  0.3524
  0.9146
  3.3722
 -0.7651
  2.9263
  0.4368
  2.3436
  1.8016
  0.0296
  1.7078
  2.9591
  1.8263
 -0.4643
  0.1687
 -0.3494
 -0.3652
  0.2029
  0.9573
 -0.2913
 -2.0186
  0.1027
  0.2976
 -0.5089
  2.1188
 -0.3065
  1.7018
 -0.1679
  1.3241
  0.1781
  1.4997
  1.5597
  2.1743
  0.9212
  0.4238
  1.2578
  1.1360
  1.0541
  1.3083
  1.8185
  0.4440
  2.4781
  0.9332
 -0.1413
 -0.5645
 -1.0847
 -1.2270
  0.2200
 -0.7476
  1.0070
 -0.4286
 -0.8195
 -1.0035
  0.8287
 -0.1998
  0.6981
  0.7189
  0.8392
  0.8504
  0.6812
  2.7543
  0.8017
  1.7689
  1.5121
  2.0977
 -0.1712
  0.0589
 -2.7085
 -0.0803
 -1.1455
 -0.9357
  1.7273
 -0.5912
  0.1391
 -1.2718
 -1.0035
  0.0510
 -1.1730
 -1.3775
  0.5556
  0.6134
 -0.1314
 -0.1251
  1.3412
 -1.9578
  1.1511
  0.1956
  2.1118
  0.5133
 -0.0286
 -0.2703
 -1.1616
 -1.4996
 -1.3609
  0.3732
  0.0102
  0.8138
 -1.0343
  0.9848
 -0.0361
  0.6139
 -0.4186
  0.2533
  2.2973
  2.9785
  2.3169
  0.9154
  2.3529
  0.0091
  1.2513
  2.1344
  2.9136
 -0.1282
  0.5228
  0.0933
 -2.4468
 -1.3882
 -0.2392
 -0.6292
  0.0665
 -0.4852
  0.2701
 -0.3113
 -0.6043
 -0.1496
 -0.6563
  1.2846
  0.4712
  2.3781
  0.7010
  0.5717
  1.3295
 -0.3253
  1.1914
  1.8079
  2.1126
 -0.9467
  1.4257
  2.0556
  1.6713
  0.7388
  2.5621
  0.4654
  0.8377
  0.1361
 -0.4877
 -2.7888
 -0.7203
 -0.9741
  1.7905
 -1.0570
  0.5480
 -0.6165
 -1.9486
  1.7378
 -1.2411
  1.4831
  0.8626
  0.6000
  1.5341
  2.8737
  4.6152
  0.3254
 -1.4916
  1.0906
  0.3442
  0.3180
  2.0353
  0.2979
  0.0027
 -1.1497
 -0.0821
 -0.1545
  0.8119
 -0.0519
  0.3644
 -0.6389
 -1.5735
 -1.0769
  0.7453
  1.2239
  1.4715
 -0.3536
  0.9820
 -0.2631
  0.8607
  1.5146
  0.7321
 -0.2036
  1.3065
  1.0757
 -1.1396
  0.3052
  0.6581
  1.5454
 -0.3322
  0.6804
  0.2558
  1.2499
  0.6250
  2.0858
 -0.1183
  0.4314
  2.0459
  2.4377
  1.4895
 -0.5289
  0.8554
  0.0720
  0.9532
  2.2229
 -1.8137
  1.6323
  1.5767
 -0.4691
  1.3429
  1.2253
  0.7270
  0.0443
 -0.4250
 -0.1670
  0.9833
 -0.9707
 -1.2279
  2.2464
  0.1439
  0.9818
  0.8588
  1.1286
  1.7404
  0.1096
  0.3598
  0.9800
 -0.4234
  1.6063
  1.2876
  1.0393
 -0.2272
  3.8226
  1.6100
  1.8736
  0.4415
  1.0612
 -1.9151
 -1.7879
  0.5562
  0.9032
 -1.7943
  0.5729
 -0.0578
 -0.5933
 -0.0268
 -0.6973
  0.8970
 -2.0243
 -0.8477
  1.9503
  1.4627
  1.6024
  0.8224
  0.0238
  2.6034
  1.7462
  1.2054
  0.4421
 -0.5109
 -0.8829
  1.2515
 -0.4351
  1.3799
  0.0560
  0.2796
  1.6042
 -0.2836
 -1.0133
  0.6435
 -0.0145
 -0.2098
 -1.1149
 -1.2215
  0.5436
 -0.4474
 -0.1633
  1.2437
  0.2008
  1.2397
  0.3678
  2.5320
  1.3002
 -0.3156
  1.7484
 -0.9156
  0.0905
 -0.1708
  1.2439
  0.3553
 -0.6932
  0.5487
 -0.9060
  0.8068
  2.6622
  0.2416
 -1.3735
  1.7908
  1.1926
  2.3957
  0.8427
  1.2618
 -1.7457
  3.0368
  1.3515
 -1.0273
  0.5567
  1.2380
  0.5456
 -1.5122
 -0.7164
  0.6443
  1.6402
  1.0674
  1.6788
  0.5386
 -0.9347
  2.5432
 -0.0623
 -0.2610
  0.0327
  1.5514
  2.2661
 -0.1366
 -0.1564
  0.8201
 -0.6340
 -0.0504
  1.8875
  1.6582
  2.4943
  1.9274
  0.1133
  2.5256
  0.8874
  0.1902
  1.8898
 -1.5148
 -0.6160
  0.3779
 -1.0688
  0.2780
  0.1257
  2.2453
  0.3698
 -0.5521
  1.3167
 -2.3077
 -0.0619
 -1.0548
  1.1282
 -1.1032
  2.7209
 -0.2402
  0.8032
 -0.0229
  0.5128
 -1.7063
 -1.0667
 -0.1037
 -1.7370
 -0.3099
  1.2911
  1.0410
 -1.1283
 -0.8052
 -0.4843
  0.3185
 -0.2639
  0.0574
 -0.6518
  1.7856
  0.7017
  1.3763
 -1.5694
 -1.0706
 -0.7860
  0.2914
  2.1443
  2.1745
  1.2730
  0.9147
  1.9978
 -0.8705
 -0.4841
  0.0601
  0.0837
  0.5276
 -0.6012
  0.3948
  0.5910
  0.9413
  0.3592
  1.6441
  2.4323
  2.3895
  0.0496
  1.5942
  1.0956
  1.6953
  0.3414
  0.1255
  1.2596
 -0.1436
 -0.1825
 -0.3101
 -0.0281
 -1.5738
  0.9110
 -2.5487
  0.6364
 -1.2549
 -1.8689
  1.4455
 -0.1778
  0.4297
 -0.1757
 -0.0856
 -0.5593
  0.9786
 -0.9297
 -0.1657
  0.9447
  0.3395
  1.5249
  0.8183
 -0.3407
  0.6033
  0.0396
  2.6663
  1.6725
  1.7660
  1.3116
  1.5840
  1.5787
  0.8338
  1.8547
  1.5272
  0.8194
 -0.4439
  0.5249
 -0.2512
  0.5774
 -0.2089
  1.3404
  1.2040
 -0.6175
 -0.0470
  0.1914
  0.9309
 -1.0461
  1.1732
  0.7399
  1.2381
  0.9567
  1.1008
 -0.9272
  0.1334
 -0.7251
  0.5104
 -0.3228
 -0.3267
 -1.0595
  0.7107
 -1.8202
 -0.3294
  0.1661
  1.4795
 -0.7052
  0.0579
 -0.8105
  2.5049
  0.4658
 -0.0009
 -0.2345
  0.3605
 -0.2602
  2.0818
  0.3850
  0.1891
 -0.6505
  0.4319
  1.6072
  2.0738
  0.1250
  0.2955
  0.7298
  0.2652
 -1.0223
  0.2508
 -0.8253
 -0.8780
  0.4652
 -0.3219
 -2.5657
  2.3486
  0.5613
  0.5636
  1.2556
  1.2846
  0.9249
 -0.7081
  1.0824
  0.8992
  2.2809
 -0.1173
  1.1040
  1.3661
  0.1882
 -0.1510
  2.1794
 -0.2334
 -0.0688
 -2.0333
  0.7443
 -1.2790
  0.7519
 -0.1113
  0.0668
  0.0422
  0.8842
  0.7450
  0.5678
  0.8790
  0.3691
  1.1833
  1.9144
  0.1327
  1.7759
  1.2955
  0.5189
 -0.2198
 -0.3445
  1.3622
  2.9614
  0.6122
  1.2408
  2.3589
 -0.6859
  0.8452
  1.6536
  0.6679
 -0.1675
  0.4667
 -0.9586
  0.0531
 -1.1757
 -1.7995
 -0.9009
  0.5971
  1.1891
  0.8808
 -0.4756
  0.8220
  0.6714
 -0.1877
  1.3360
 -0.7121
 -0.8585
 -0.1262
 -0.5406
 -1.4576
  1.2475
  2.7312
 -0.4487
 -0.1724
  0.0376
  0.1288
  1.0589
  0.5793
  0.3341
  1.3285
 -0.7647
 -1.4775
  0.4156
 -0.9059
  0.4384
  1.9932
  1.1936
 -0.6631
 -0.9759
  2.3656
  2.0118
  1.4101
  2.3368
 -1.0131
 -1.6445
 -0.1663
 -0.4308
  0.5876
  0.8486
 -2.0294
 -0.0474
  1.1955
  1.5045
 -2.7642
 -0.8307
  0.3372
  0.0741
  0.0673
  0.6248
  1.6686
  2.5396
  0.6571
  1.0617
 -1.5758
 -0.6510
  0.8586
  1.0049
 -0.9827
  0.9982
  1.0526
  0.6707
 -0.1445
  2.0595
  1.0661
  1.1689
  0.7007
  2.6415
  1.2469
 -0.5563
  1.6905
  1.8098
  1.4136
  1.4584
  0.3213
  0.9871
 -0.1319
  0.1478
 -0.2129
  2.0778
  1.1290
  0.4284
  1.3101
  3.2174
 -0.4056
  2.2867
  0.5563
  0.2382
 -0.8927
 -0.0332
 -0.4839
  0.4766
  0.1318
 -1.2477
 -0.9183
  1.1804
 -0.3540
 -0.3856
 -0.2943
 -0.3316
 -0.9172
 -2.1999
  1.0969
 -0.1623
 -0.1277
  0.0829
  0.5298
 -1.3750
  1.3467
 -1.5676
 -0.7089
  0.2999
  1.7410
  0.0684
 -1.0917
 -0.4124
 -0.6716
 -0.0780
  0.0877
  1.2390
  0.6776
 -1.1447
 -0.6246
  0.3101
  2.7798
  1.7712
  1.9917
  1.0760
  2.4581
  1.3687
 -0.0313
  0.4188
  1.2145
  0.0110
 -1.4707
 -2.0623
  0.0074
  1.0343
  0.0881
 -2.1362
  0.5692
 -1.2488
 -1.1724
  0.5727
  2.2452
 -0.6225
 -2.0080
 -1.0259
 -1.2227
  1.5336
  0.9941
  1.1884
 -0.3185
 -0.5285
  0.3122
  0.8740
  1.2969
 -0.5322
 -0.4697
 -0.3784
  0.5856
  0.2398
 -2.5857
  0.1646
 -0.3323
  0.7554
  1.3710
 -0.2287
  1.9004
  0.2505
  0.6170
  2.2324
  0.6648
  2.6378
 -1.5283
  0.8812
 -0.0016
  2.3452
  0.5148
 -0.0713
  0.5941
  0.2488
 -1.2759
 -0.8708
  1.1026
  1.3038
  0.3686
  2.0125
 -1.2109
 -0.2435
 -0.8551
  0.2237
 -2.1551
 -1.2617
 -0.1814
 -0.6887
  1.6733
  0.1248
 -0.1435
  0.8005
  1.4497
 -2.4224
 -0.1854
 -1.8093
  1.0381
  1.5290
  0.0296
 -0.4424
  0.9532
 -0.5332
  0.1711
 -1.0541
 -0.8792
 -1.2939
  0.2048
 -0.2752
 -0.9949
 -1.2027
  0.7505
 -0.5388
 -1.2520
  1.7399
  1.2564
 -1.0774
  1.6137
  0.8674
  2.1428
  0.0621
  0.3525
  1.5634
  0.7181
  4.1093
  1.8299
  1.4035
 -0.6012
 -0.0461
  0.6944
  0.1501
  0.7773
  2.7251
  0.8419
 -0.0445
  0.3844
  0.1257
  0.2509
 -0.1927
  0.2624
  2.3998
 -0.3449
 -0.1012
  0.6477
  0.6141
  0.5105
  0.4304
  1.1784
  0.1750
 -0.1418
  0.5238
 -0.5045
  0.4895
  0.6449
  0.4762
 -1.5874
  0.8908
  0.1367
  0.1401
  0.6874
  0.9721
 -0.6621
  2.9853
  1.6893
 -0.7317
  1.8278
 -0.2402
  1.2504
  1.0622
  0.9736
 -1.3359
  2.1605
  0.0733
 -0.0655
  1.2681
  0.6620
  0.3035
  0.9244
  1.0622
  0.9148
  2.3002
 -0.9350
  1.1791
 -0.9129
  0.3229
  0.0173
  0.6217
  0.2690
 -1.1902
  0.3733
 -0.2963
 -0.6762
 -0.1446
  0.8816
 -0.7125
 -0.3082
  0.2778
  1.1339
 -0.8267
 -0.8350
 -0.9495
 -0.5669
  1.5095
  0.6008
  0.7836
  1.4677
  0.4184
 -0.5140
 -0.0520
 -0.9885
  0.8411
  0.2726
  0.0937
  0.9551
 -0.4714
 -0.8487
  0.6775
  1.5624
  0.7864
  0.5453
  0.6111
  1.3695
  0.1839
  1.4480
 -0.5624
  1.3850
  2.2048
  0.2889
  1.4752
  0.2218
  0.7690
  0.5425
  2.0857
 -0.6133
 -1.2655
 -0.0907
  2.5236
 -1.4175
  1.6649
  0.3408
 -1.3326
 -0.6585
  0.6337
 -0.4866
  0.1367
 -0.0094
 -2.5851
  0.6081
  0.2245
 -0.8942
 -0.7335
  0.3503
  0.7264
  1.2610
 -0.2299
 -1.1822
  0.0287
  0.1843
  0.4817
  1.6043
  0.4034
 -0.5484
 -0.1049
  0.3281
  1.8100
  1.5436
  1.2632
  1.2823
  0.8096
  0.6684
  0.1149
  0.1797
  1.5124
  2.5602
  0.4390
  3.1097
  0.9639
  0.4047
  0.5273
  0.0539
  1.0981
 -0.3709
  1.3677
 -1.0253
  1.1299
  0.0452
 -0.5894
 -1.4102
 -1.8521
 -0.8963
  0.7087
 -0.9731
 -0.9171
  0.2703
  1.2716
 -0.3857
 -0.4808
 -0.9022
 -2.3373
 -0.8615
  2.5831
  1.2623
 -0.7800
 -0.2734
 -0.5008
 -0.7603
 -1.0637
  0.0408
 -0.2525
 -1.8259
 -1.9040
  2.8300
 -0.5828
 -0.7817
  1.1466
  1.5891
 -0.6634
 -0.0376
  0.6098
  2.3274
 -1.3072
  1.4453
 -0.3539
  1.0387
  2.3210
  0.2164
 -0.0512
  1.0416
 -1.1322
  1.3994
  0.1618
 -0.4685
 -0.4260
  0.1864
  1.0405
  0.7070
 -1.3748
 -0.2850
  0.0301
  0.5480
  0.8614
 -0.6602
  0.0586
  0.4648
  0.9143
  0.6278
  0.8945
 -0.2835
 -0.1280
  0.5032
 -0.3561
 -0.9844
 -0.9123
 -0.3271
 -0.4946
 -1.9316
 -1.0650
  0.3433
  0.8560
 -2.3705
  1.2524
  1.1073
 -0.4768
  1.1053
 -0.5795
  1.4758
  1.2778
  2.6811
  0.1032
  1.5878
 -1.9016
  0.8425
  1.7998
  1.8279
 -0.5679
  1.9866
  0.3298
  0.7666
  0.9611
  0.0546
  1.2502
  3.0614
 -0.1077
 -0.4125
  0.0060
  1.2278
 -1.3299
  1.9753
  0.8207
  0.7763
 -0.2398
 -0.7794
 -0.0951
 -1.0404
 -0.0515
 -1.1715
 -0.6800
  0.5696
 -0.8748
 -0.4915
 -0.8024
  1.9198
  0.5329
  1.7671
  0.3995
  0.2773
  0.2352
  0.7763
 -0.5251
  0.4493
  1.0400
 -0.1670
  0.2457
 -1.5845
 -1.2015
  0.2106
  0.2462
  2.1305
  0.7627
 -0.4335
  2.2914
 -0.1562
  0.7089
  0.4896
  0.6132
  0.4557
  0.1191
  0.3385
  1.1331
 -0.0897
  2.0673
 -1.5316
  2.4980
  0.4469
 -1.2111
  0.2284
  1.3750
 -0.8339
 -1.8254
 -0.6836
  0.2627
  0.6162
 -0.1887
 -0.6549
 -0.9693
 -0.1304
 -2.8799
  1.8671
 -1.2661
  1.6267
  0.7188
 -0.8474
  0.8119
  0.8572
  2.4356
  1.2091
 -0.3201
  0.1602
  1.9857
  0.1638
  0.4691
  0.0682
  2.5037
 -2.3667
 -1.1606
  1.0262
  1.3389
  2.6447
  0.9111
  2.6538
 -0.9493
  2.1353
  1.0190
  2.4914
 -1.0946
  0.4496
  1.5938
  2.2132
  1.7921
  0.8164
  1.3060
  1.4510
  0.5544
  0.8593
 -0.7098
  0.2214
  1.4187
  0.2556
 -0.7800
 -0.2200
 -2.1658
  1.1346
 -2.4427
 -1.2083
  0.2007
 -1.4389
 -0.8721
  1.5944
  0.2817
  0.6818
 -0.7948
  1.0348
 -0.5761
  1.5437
 -1.1753
 -0.7432
 -0.0944
 -1.8104
 -1.7786
 -1.3511
  0.9714
 -0.5264
  0.4351
  1.5908
 -1.3202
  1.8366
  0.7764
  0.1355
  1.9246
 -0.9667
  0.6933
  2.4176
  1.0743
  0.1786
  0.5257
  2.3881
  2.1307
  2.2134
  1.6828
  2.8384
  0.1415
  0.0483
 -0.4966
  0.4814
  0.1879
  0.7531
 -0.5682
  1.4103
  0.1395
  0.6136
 -1.0830
 -1.0370
  1.1686
  0.5193
 -1.1824
 -1.6364
  1.0166
  0.2089
  0.1066
  0.7999
  0.7277
 -1.8428
 -0.3036
 -0.6235
  0.1978
  0.0351
  1.2308
  0.5980
  0.7500
  0.3467
 -0.7713
  0.1032
  0.2806
  0.1474
  0.5430
  1.1373
 -0.8821
  2.3855
  1.8416
  1.1405
  0.8331
  1.1367
  1.0591
  0.2989
  0.6806
 -0.8798
  0.3368
 -0.0218
  0.9330
  0.7753
  2.3596
  2.5869
  0.8880
  0.6289
  0.3519
  0.5241
 -0.1630
  1.6432
 -0.7393
 -0.5157
 -0.6993
  0.0247
 -1.2275
  1.0135
 -0.0392
  0.9447
 -1.2222
  0.0622
 -1.5319
 -1.6490
 -1.0132
 -1.5223
  1.3779
 -0.2031
 -1.1374
 -0.1991
  1.4991
  1.3708
  0.2287
  0.8177
 -1.6370
 -0.8386
 -0.6617
 -1.0469
 -1.7094
 -0.3518
  2.0082
  2.1987
  0.9022
  1.6956
  0.9024
  0.3363
  1.9988
  3.1719
  2.0578
  1.0924
 -1.5800
 -0.3581
  2.1690
  0.3578
 -0.4256
  0.7790
  1.5966
 -0.2801
 -1.4040
 -0.3393
  2.4542
  0.2410
 -0.8616
  0.5047
 -1.1660
  0.2752
 -1.9678
  2.1176
  0.3063
  2.0649
  0.5118
  0.3944
  0.3766
 -0.4974
  0.7983
 -1.1206
  0.0232
 -1.6519
 -2.2597
 -0.6138
  1.0126
  0.6073
 -1.5402
 -0.1180
 -2.0069
  1.0011
 -0.4675
  0.1056
  0.1086
 -0.0997
  0.2817
  0.7510
 -0.2338
  0.3239
  0.5122
  1.0781
  1.7004
  0.9181
  0.2085
  1.1008
  1.8142
  0.7279
  0.4890
  1.0432
  3.9534
  1.7740
  0.4193
  1.4822
  1.1313
  1.3979
 -0.3517
 -0.2146
  0.7295
  2.5247
  0.3424
 -3.0290
  1.7803
 -0.0068
  2.2554
  2.1833
 -1.1047
  0.6960
 -1.5001
 -2.8691
 -0.2202
 -0.2046
  0.4870
  1.2064
 -0.2718
  1.3491
  0.9400
  0.4245
 -0.2373
  0.0024
  1.5740
  0.5525
  1.9015
  1.2831
 -0.3048
 -1.2667
  0.9887
  0.7050
  1.4009
  1.4053
  1.1251
  0.0884
  2.4145
  1.6975
  1.7612
  0.0777
  1.8167
 -1.8996
  1.7641
  0.4360
 -0.1409
  0.4694
 -0.1546
  0.3068
 -1.6933
 -0.6569
 -1.4656
  1.4286
  0.3268
  2.3784
  0.8269
  2.2460
  2.6619
  1.1758
  0.1670
 -0.1014
  1.3257
 -0.4508
 -1.1566
  0.6560
 -0.1305
  0.5109
  0.1156
  0.1831
  1.9790
 -0.5713
 -0.8565
 -0.1032
 -0.3607
 -0.3633
 -0.5166
 -0.7855
 -2.2175
 -1.2923
 -1.9774
 -0.1146
 -1.0943
 -0.2382
  0.6909
  2.0736
 -0.0429
  2.4465
  0.1638
  1.3605
  2.4565
  0.7868
  0.6654
  1.4390
  0.4217
 -0.8263
  1.5309
  2.5850
  0.6970
  1.4313
  1.8700
 -0.1607
  0.4506
  0.6760
 -0.6519
 -1.2972
 -0.7290
 -0.5808
  0.0557
  2.0060
  1.3089
  3.2220
  0.2811
  0.2545
  0.3882
  2.1302
 -0.7054
  0.9565
  1.0202
 -0.4532
  0.0388
 -0.7614
  0.0840
  0.8318
  0.7807
 -1.4262
 -0.1679
 -0.3652
 -0.1566
 -1.5369
 -1.9274
  0.5408
  1.3152
  1.4461
  1.4341
  0.1836
  1.1794
  1.7963
  2.5927
  1.1058
  0.6760
  0.2778
  0.6143
  1.0753
  1.5708
  2.1131
  1.2524
  0.0630
  1.8935
  2.4876
 -0.3445
  1.1399
  2.0593
  2.4471
  1.8894
  0.4167
  0.9558
  0.9547
  0.0918
  2.5482
  1.7848
  1.5278
 -0.3203
 -0.3153
  0.1149
  1.7150
  0.7447
  0.8318
  0.4365
  0.6837
  1.3691
 -0.2677
 -0.1203
  0.9249
  1.3273
 -0.3527
 -0.5955
 -0.8742
 -1.1895
  1.4514
  1.5785
 -1.0223
 -0.1037
  0.6512
 -0.0213
  2.1965
  0.4879
  0.9785
  1.5357
  0.3807
  0.8011
  1.6015
 -0.0382
  2.0178
  0.8292
  1.5535
  0.1925
 -0.3195
  0.8637
  3.1900
  1.0442
  0.4858
 -2.2337
  0.4900
 -0.4568
  0.5395
  3.6075
  1.2147
  0.2363
  2.2812
 -0.3705
 -0.1266
  2.6913
  1.6422
  1.5983
  2.2033
  1.4559
 -0.1368
  0.5373
  2.0413
  1.2436
 -1.1261
  1.6445
  0.2372
 -1.4504
  0.6318
 -0.0852
  1.4011
 -0.2761
 -1.9728
 -0.1600
  1.1865
  0.6565
  1.3035
 -1.0884
  0.9418
  2.4133
  2.5185
 -0.8409
  0.4272
 -0.1550
  1.7252
  0.9220
  1.6063
  2.4249
  2.7513
  0.3321
  0.5372
  0.7599
  1.3010
  0.0728
 -0.7918
  0.6728
  0.9108
  0.7444
 -0.3882
  0.1623
  0.7309
  1.5239
  0.5749
  2.1484
  0.8095
  0.6177
  1.7905
  1.5199
  1.5589
  1.4160
  0.8412
  0.6336
  0.7715
 -1.3866
  0.2721
  0.8178
  0.4745
  0.9146
 -1.5404
 -0.4793
  0.7899
  0.4836
 -0.9432
  0.5178
 -1.1976
 -0.4551
  1.1899
  1.1341
  1.3444
 -0.8741
  2.0227
  1.1109
  1.8743
 -1.1995
  0.9965
  0.7236
  2.5713
  1.3223
  0.4222
  1.0220
  2.0492
  2.2589
  0.3491
  3.5567
  1.3969
 -0.2123
  1.0011
 -0.5204
 -0.6746
  1.6589
  0.6162
 -0.0779
  1.0096
 -0.2417
  0.7588
  1.4566
 -0.0960
  0.3156
  0.9688
  0.4074
 -1.4244
 -0.3846
 -1.1487
 -0.7138
 -0.3906
 -1.6159
 -0.0647
  0.2410
 -2.6450
 -0.0001
  0.1121
  1.5122
  0.7052
  0.2318
  0.1255
  2.2568
 -0.0251
  0.7438
  0.4189
  0.0443
  0.9016
  1.3798
  1.3389
  0.2954
 -0.1183
  1.0126
  0.9725
  2.3483
  0.2727
  0.7251
  2.3638
  0.6261
 -0.3197
  0.5951
 -1.1476
 -0.1176
  0.1985
 -0.7886
  2.3662
 -0.2293
 -0.7119
  1.6512
  0.4926
  0.7664
  2.0739
 -0.6830
  0.9393
  0.9156
  1.4836
  0.5280
  1.2822
  0.0847
 -2.2530
  0.3455
 -0.3340
 -1.9006
 -1.1503
 -0.6720
 -1.7613
 -0.4014
 -1.4630
  0.2531
  0.0810
  1.5167
  0.0283
  0.1947
 -0.4759
 -0.0634
 -0.5899
  1.1610
  0.5366
  0.5616
  2.0512
  1.6797
  0.2079
  1.0716
  2.0875
  0.6232
  1.4753
  1.1020
  0.6163
  2.1125
  2.2050
  0.2453
 -0.9628
  2.1006
 -1.0271
  2.0363
  0.7518
 -0.1651
  0.5574
  1.7438
  1.7951
 -0.2367
  1.8190
 -0.2708
  0.8128
  0.5486
 -0.6843
  0.0736
  0.0080
 -2.0500
  1.1782
  2.5377
  0.8436
  2.4607
  0.9662
 -1.7296
  1.4476
 -1.0469
 -0.2788
 -0.0977
  2.3824
  1.2263
  0.5102
  2.0045
 -0.1776
  1.0458
  0.9992
  2.5406
  1.8700
  0.3401
  0.5620
  1.6040
  1.7137
 -0.7524
  2.2109
  1.8613
  1.7520
  2.1536
  1.1970
  0.6149
  0.1917
  2.7364
  1.1074
  0.1805
  1.0615
  0.6933
  0.1397
 -1.1712
  1.0117
  1.1835
  1.3806
  1.0346
  3.3119
 -0.7054
  1.2018
  1.5925
  1.2365
  0.2009
 -1.2284
 -0.2980
  2.2123
  0.3065
 -1.7554
 -1.2285
  0.9219
  1.6115
 -0.6240
 -0.2451
  0.1927
 -0.6966
  0.5250
  0.8050
 -0.5890
  0.8243
  0.9880
  0.2970
  2.7369
  3.0771
  1.4392
  0.3136
  0.0564
  2.5255
  1.9139
 -0.1564
  1.7842
 -0.0859
  1.0788
  3.6195
  1.1751
  1.2512
  0.7965
 -0.6105
 -0.7631
  1.0641
  1.7112
  2.6934
  0.9615
  2.3921
  1.3466
  0.9167
  0.4597
  0.6313
  1.3860
  1.6330
  1.8394
  1.4232
  0.1067
  0.6246
  3.3063
  0.3769
  0.3559
  0.6326
  1.9637
  2.5329
  0.9933
 -0.0603
 -0.9094
  0.2024
  0.8025
  0.5797
  2.8087
  2.0391
  1.5771
  0.6897
  2.4220
  1.4883
  2.1350
  1.2680
  1.6852
 -1.1634
  1.9184
  3.0105
  0.9252
  0.3793
  1.6952
 -0.3842
  2.4135
  1.2687
 -0.3265
  1.2123
  3.6667
  1.1422
  4.1583
  1.1804
  0.0697
  1.5317
  0.6330
  1.4636
 -0.5870
  0.2295
 -2.7810
  1.4904
  1.3241
  0.7068
  0.1404
 -0.3683
 -0.2811
  1.6879
  1.3976
  1.1107
 -2.1012
  0.1788
  0.5867
 -1.1251
 -0.9719
 -0.2507
 -0.1976
  0.7961
 -0.1011
 -0.1873
  1.4933
  0.6717
 -1.7111
  0.6295
 -0.1345
 -0.2950
  0.1317
 -0.8205
  0.2008
  0.4435
 -0.4402
  1.2043
 -1.7951
  1.0133
 -0.1132
 -0.4772
 -0.7441
  1.3180
  1.1072
  0.1914
  1.5379
  0.7883
  1.8863
  1.4608
 -0.3922
  0.7022
  1.8992
  0.6223
 -0.4056
  0.7187
  1.1001
  0.2947
  0.9188
  2.5613
  1.6233
  0.6230
  0.5370
  0.7731
 -0.3713
  0.6104
  2.1583
  2.9739
  2.3722
  0.6846
  2.1690
 -1.3948
  0.8317
  0.1133
  0.4284
  1.3253
 -0.1530
 -1.1646
  0.6160
 -0.3476
  1.3459
  0.3412
  1.4531
  0.0546
  0.2853
  1.5933
  1.4144
  1.7549
  3.2831
  0.8955
  1.0700
  1.1593
  2.0613
  1.5851
  0.6017
 -0.3709
 -1.1217
  0.4777
  0.5616
  2.2163
 -1.5366
  0.6238
  2.2403
 -0.0896
  1.7744
  2.4752
  0.7637
  1.1835
  4.2782
  0.3456
  1.2765
  1.1891
  0.7528
  0.4379
 -0.2483
 -1.2538
 -0.3212
  0.4684
  1.3331
 -1.3317
  1.0057
 -0.0605
 -0.0565
  0.2528
 -0.7023
 -1.2175
  0.2933
  0.0635
 -1.4472
  0.5725
  0.3583
 -0.0692
  0.0645
 -0.8393
  2.1520
  0.2587
 -1.4848
  0.2393
 -0.6360
 -0.0048
 -0.0822
 -0.0704
  0.7535
  3.3300
  1.4549
  1.2534
  3.9956
  0.2126
  1.6557
  1.6958
  3.3921
  1.5410
  1.0852
 -0.3171
  0.5628
  0.7051
 -0.4023
  2.5380
 -0.7316
  0.6016
 -0.0737
 -0.0652
  3.4249
  1.0322
  0.4383
  2.3610
  0.7103
  1.8130
  1.0490
  1.3878
  0.4098
  1.4541
  0.7335
 -1.1620
 -0.5545
 -0.7590
 -0.4185
  1.8087
  1.4589
  2.2320
  0.7842
 -0.7336
  0.2417
  1.8938
  1.4083
  1.8575
  0.0187
  1.6619
  1.0644
  1.9125
  1.2717
  1.4322
 -1.0448
  1.7091
  1.7339
  1.8596
  0.1826
 -0.2109
  0.9217
 -1.4333
  0.7518
 -0.6061
  2.8039
  1.8160
  0.4236
 -0.1114
  2.0131
  0.8417
  1.2008
 -0.3165
  1.2968
  1.2328
  0.3287
 -0.4223
  0.1065
  0.7485
  1.0702
 -0.6836
  0.4013
  1.7988
 -0.7294
 -0.4796
  0.8269
 -0.8313
  2.3283
 -0.8119
  0.5050
  0.4490
  0.1240
  0.6587
  2.5473
 -0.2709
  1.6508
 -0.6118
 -1.2557
 -1.0051
 -0.5238
  0.6407
  0.3043
 -0.7999
 -0.5542
 -0.8667
 -0.7202
 -1.0474
 -1.7744
 -0.4311
  1.0509
 -0.0726
  2.1924
 -0.4181
  2.3655
  1.0434
  1.7671
  2.8424
  0.9055
  1.4651
  1.0980
  2.3925
  1.9645
  1.8430
  1.8590
  0.8992
 -0.1666
 -0.0267
  3.3332
  0.4712
  0.0574
 -0.9967
  1.5285
  1.7444
  1.3496
  0.1589
  2.2447
  2.6841
 -0.0563
  0.4331
  1.6704
  0.7083
  0.5933
  0.6343
  0.7352
  1.5578
  1.2067
  1.8849
 -0.1624
  1.9547
  1.3409
  2.0659
  3.2664
  0.3838
 -0.1160
  0.8889
  1.7345
  0.1049
  2.0109
  1.9141
  0.7172
  2.2972
  0.5879
  2.1940
  1.2315
 -0.6789
  0.8192
  1.0306
 -0.0652
 -0.1538
  0.7237
 -0.0736
  2.1162
  0.4735
 -0.4227
  0.0134
  0.4756
  1.4824
  0.4141
 -0.3119
 -1.4840
 -0.9822
  1.3585
  1.0329
  0.9947
  0.9809
 -3.2994
 -1.4000
 -0.8955
 -1.7504
  0.8951
  0.9800
 -1.0399
  0.2310
  0.6805
  1.3201
  1.4810
  1.3659
  0.5322
 -1.4374
 -0.6232
  0.2636
  0.0958
 -1.4644
 -0.4962
 -0.4264
  1.5996
 -1.0735
  0.4903
  2.0696
  1.1214
  0.5764
  2.8603
  0.6775
  0.9182
  1.5399
  2.6299
  1.7732
  0.9607
  0.5593
  0.7897
  1.7862
  0.3890
 -0.1437
  1.7401
  3.0737
 -0.6339
  0.9027
  0.8250
 -0.0839
  0.1897
  0.1838
  2.4934
 -0.8418
  1.2945
  2.1630
  0.4954
  1.7948
  1.5425
 -1.2011
  0.8847
  1.4028
  0.8762
 -0.6764
  2.0860
  1.5080
  1.4310
  1.0459
  3.7817
  1.0393
 -0.2990
  2.5333
  1.2869
  1.0667
  0.1407
  0.6645
  3.0090
  0.6816
  1.3167
 -0.1113
 -0.1123
  1.3168
  2.2562
 -0.1659
 -0.9364
  1.1800
  0.9826
  1.0573
  1.1225
  2.0914
  0.0566
  1.9268
 -0.8699
  0.7616
  1.3120
  0.1762
 -0.3650
 -0.2003
  1.3197
  0.3012
  0.2965
  0.5135
  1.6566
  0.8306
  0.9250
 -1.9910
  0.7904
 -0.6681
 -0.0319
 -2.2776
 -0.9250
  1.2485
 -1.1579
 -0.0357
  1.5430
  1.6152
  0.2934
  0.7041
  1.0493
  1.4052
 -0.8564
  0.7181
 -0.3793
 -0.1715
  0.5695
  1.5627
 -2.5697
  1.8632
 -0.7123
 -0.6778
  0.0388
  1.2461
  0.4994
 -0.6541
  0.9579
  1.6682
  1.2712
  0.4003
  3.5584
  3.0911
 -0.0663
 -0.5450
 -1.1148
  1.1167
  2.4734
  2.4576
  0.7286
  2.5076
  1.0076
  1.3406
  2.2050
  2.6386
  0.6022
  0.9318
  1.9430
  0.6036
  1.2656
  0.2072
  1.1448
  1.1832
  0.7399
 -0.1388
  2.1472
  1.9090
  1.0234
  0.4493
  2.8071
  1.1662
  2.1906
  1.4955
 -0.2005
  2.0628
  1.2291
  2.4841
  2.1764
  1.0835
  0.6043
  0.9086
  0.3901
  0.0116
  2.0253
  0.4664
 -1.0598
  1.9451
  1.5973
  0.7927
 -0.3480
  1.0534
 -0.4604
  1.5064
  2.3487
  1.5655
  1.6661
  0.5802
  0.1833
 -0.0675
  0.3544
 -0.4058
 -3.9587
  0.9690
 -2.4282
  0.9073
  1.3234
 -1.0914
  2.7827
  2.2407
  2.6364
  0.4946
 -0.8382
 -0.5109
  0.4852
 -1.2194
  0.0272
  1.0566
  1.1138
 -0.8091
 -0.6829
  0.2370
  0.3551
  0.5114
 -0.8012
  0.0377
  1.2928
  1.8864
  0.6097
  1.4643
 -0.0400
 -0.2237
  0.7431
  1.3458
  0.9795
  2.2115
  2.6931
  0.1479
  2.2792
  1.8277
  1.4342
  1.0089
  1.3346
  1.3828
  0.4594
  1.1650
  0.3702
  0.4777
  0.7657
  0.1199
  0.4569
  0.7162
  0.7324
 -0.2625
  1.0516
  0.6681
  2.2326
  3.1997
  1.2317
  0.3601
  1.7360
 -0.4292
  1.9749
 -0.2575
  1.3095
  0.6501
 -0.1747
  0.7794
  1.5701
  2.1137
 -0.1354
  1.2706
  1.4610
  0.9718
  1.3251
  0.8518
  0.8916
  2.2863
  1.4585
  2.5207
 -0.6135
  0.7048
  1.7473
 -1.5392
  2.0714
  3.1895
  1.4546
 -1.6425
 -0.0358
  0.3510
 -0.1152
 -0.0801
 -0.6517
  1.7559
  0.2076
  1.9721
 -0.7852
  0.2545
 -0.5781
  0.9092
 -0.1836
 -1.3960
 -1.6079
  1.1336
 -0.4133
  0.9497
 -1.0794
 -0.8729
  2.4394
 -0.1424
 -0.4506
  0.8522
 -0.7404
 -0.4343
 -0.6164
 -0.6049
  0.2333
 -0.0906
 -0.1192
 -0.4990
 -2.0988
 -1.5317
  0.4924
  0.2711
  0.2512
 -0.0676
 -0.2603
  0.5041
 -0.5468
 -0.8620
 -1.2205
  1.6303
 -1.0751
  0.6994
  1.4156
  2.6262
  0.4046
  0.6797
  1.1289
  1.7058
 -1.6870
  0.3260
 -0.1954
  1.5122
  1.2911
  2.5485
 -0.1609
  0.2294
  0.0711
  0.1760
  2.1563
  0.0831
  0.7694
  1.7303
 -0.0946
 -0.4560
  1.2327
  0.3626
  1.3488
  0.8783
  1.4720
  1.0835
  1.8756
  0.7260
  1.2183
  0.8841
  3.0295
  1.2811
  0.8900
  0.8123
  1.7208
  2.6783
  0.2265
 -0.3874
  0.3815
  0.6881
  0.1342
 -0.3367
  0.5071
  0.4824
  1.0559
 -0.1809
  1.5086
 -0.2763
 -0.8619
  2.0249
  2.8669
 -0.7949
 -0.3136
  0.5457
 -0.4372
  0.1863
  0.1102
 -0.2743
  0.1246
 -0.3775
  1.2753
  1.2338
  0.6563
  2.0406
 -0.9188
  0.9321
 -0.2415
 -0.8236
 -1.0052
  0.4128
 -0.0342
 -1.5529
  1.2304
 -0.4838
 -0.5926
 -1.0883
 -0.1918
 -2.1796
  1.6299
 -0.0851
 -0.2333
 -0.0168
  0.5454
 -0.6362
 -0.0037
  0.7189
 -0.2035
 -1.4344
  0.3145
 -0.5030
 -0.9449
 -0.4547
 -0.5524
  0.8572
 -1.9878
 -0.5206
  0.6362
  1.5577
  0.5473
  0.4768
 -0.0125
  1.9459
 -0.4700
  1.2436
  1.3240
 -0.3528
  1.0018
  1.5064
  0.1063
  1.0337
  1.4913
 -0.0014
  1.5170
  1.9135
  2.3476
  1.2571
  0.4212
  0.8699
 -0.9907
  0.5879
  1.7803
  1.5403
  2.8077
  0.0300
  0.3416
 -1.9730
  0.5154
  0.8292
  0.2321
  0.9989
  0.9013
  1.9267
  2.6401
  0.1975
  1.1342
  1.4586
  1.2044
  1.1510
  1.5620
  1.3161
  1.5379
 -0.1902
 -0.4972
  1.3033
 -1.0475
  1.3869
  2.4442
 -0.8549
 -0.7446
  1.0346
 -1.0097
  0.3420
 -0.0944
 -0.7528
  0.7792
  2.1175
 -1.2090
 -1.4245
  1.0257
  0.3612
 -1.3139
  0.7999
 -2.4699
 -2.2220
 -0.3134
  0.2041
 -0.8017
  1.0886
 -1.1885
  0.2460
  0.6301
  0.4935
  0.1647
 -0.7054
 -0.0322
  0.7423
  0.7523
  0.4327
 -2.0284
 -0.9726
  2.0705
 -1.2213
 -0.9801
 -0.4131
 -1.0024
  1.3926
 -0.9747
  0.8055
 -0.7511
 -1.1662
 -0.5662
 -1.9703
 -1.1751
  0.0733
  1.1310
 -0.2048
  1.2979
  0.0527
 -0.6483
  0.8074
  1.0942
  1.9163
  0.7792
  0.8911
  0.6923
  3.6270
  2.9035
  0.6256
  1.1786
  0.6935
  1.7276
  0.7753
  0.7733
 -0.5144
  1.0054
  1.4909
  0.9977
  0.8699
  1.8386
  2.4775
 -0.5326
  2.0487
  0.3980
  0.1573
  2.9633
  1.2661
 -0.1227
  1.3462
  0.8161
  1.0863
  2.4155
  1.7833
  0.5597
  1.7542
  0.6897
  0.5562
  1.0769
  1.2199
  0.0115
  2.9759
  0.7658
  0.8374
  4.7354
  0.8598
  1.0798
  1.0775
  0.6777
 -0.3711
 -0.4221
  1.0529
 -2.4737
 -1.5995
 -1.1972
  0.0227
 -1.1283
  0.2395
 -0.4790
 -0.1069
  0.2791
  0.9870
 -1.3258
  0.3809
 -0.2091
 -0.5935
  0.0944
  0.3462
  1.8525
  0.7827
 -0.2520
  0.6225
 -0.1012
  1.2569
  0.5964
  0.6867
 -0.3962
  1.2350
  2.7589
  0.6149
  0.0789
  2.1217
  0.1682
 -0.1135
  0.8293
  2.5437
  1.7251
  0.8352
 -0.3320
  1.0130
  2.6719
  0.7401
  2.6300
 -0.6545
 -0.2737
  1.1420
 -0.2153
  0.1532
 -2.1660
 -0.4186
 -1.9372
  1.6734
  0.7798
  1.0421
  1.9898
  0.3765
  2.0719
  0.6522
  1.1837
 -0.0605
  1.6067
  1.1597
 -1.7363
 -0.9675
  1.2701
  0.5462
  1.5442
  0.4528
  0.2765
  0.7727
 -0.1464
  1.0219
 -0.5018
  2.0440
  2.1696
  1.2428
  0.4057
  0.2503
 -0.8549
  1.3046
  0.4229
  0.8232
  1.4432
  1.0300
  0.6112
  1.1696
  1.2268
  1.0786
  0.5347
  1.2781
  2.7335
 -1.4156
  0.4681
  2.0696
  1.9464
  0.9222
  0.5054
  1.5725
  0.1327
 -0.0485
  1.1718
  0.6084
 -0.0634
 -0.0417
 -0.2004
  1.1772
  0.9247
  1.1867
 -1.5809
  0.5420
 -0.5732
 -1.7712
 -0.3219
  0.6234
  1.0630
  1.5294
  1.4950
 -0.2192
 -1.9674
 -1.2938
  1.5877
  1.7793
 -0.4378
  0.8267
  0.0854
  2.0508
 -0.5518
  1.9337
  1.9473
  2.0667
  1.3495
  0.4199
 -0.2186
  1.2866
  0.7835
  1.7853
 -0.0881
  0.9649
 -0.5412
 -1.0983
 -0.0871
  0.3615
  0.4577
  2.1631
  0.3619
  0.0472
 -1.2075
  0.8302
  0.8312
  0.6685
  0.3927
  0.9277
 -0.1565
 -2.4385
  1.2796
  1.5091
  1.6049
  1.0885
  1.2390
 -1.0220
  0.6685
  0.8669
 -1.2896
  1.9791
  1.7490
 -0.4088
  2.1078
  1.8071
 -0.4917
  3.3501
  1.5288
  2.1338
  0.8997
 -0.9526
  1.0645
  0.6565
  0.9335
  0.3434
 -0.2650
  0.0228
  0.4986
  1.3865
  0.7464
  1.6891
  0.3796
  2.3242
  0.5708
  0.3123
 -1.1807
  2.8702
  1.4806
  1.4555
  1.1035
 -0.0128
  2.4542
  0.9277
  2.5999
  0.0461
 -0.1301
  0.5678
  0.0326
  0.5299
  0.6816
  1.5924
  1.1954
 -0.2070
  1.2878
 -0.8350
 -1.3921
 -1.0843
 -1.4345
 -0.7199
 -0.2935
 -1.4171
  0.8109
  1.0708
 -1.5613
 -0.2578
  0.1479
  2.2286
  0.5114
  1.6444
  0.7622
  3.5287
  2.1546
  1.5346
  1.8513
  0.1163
  1.5295
 -0.1783
  1.1505
  1.6015
  0.7379
  2.1928
  1.9946
 -0.0598
  1.1766
  1.8685
  1.1123
 -0.0566
 -0.7061
  1.4084
  1.5155
 -0.9701
 -0.4922
 -1.4238
  1.0880
  2.4125
  0.6051
  0.5565
  0.7106
  0.5504
 -1.5626
  0.5626
  0.6671
  1.9725
  0.1933
  1.0999
 -0.6818
  2.1827
  0.1012
  0.1390
 -0.3091
  0.7187
  1.2948
  0.6132
  1.9241
  0.6837
  1.7421
  1.1868
  2.2481
  0.2149
  1.4740
  1.4402
  0.8230
 -1.2135
  0.2403
 -1.6208
 -0.3362
  0.9430
  0.6296
  0.8887
 -0.0175
  1.6591
 -0.4424
  0.9466
  0.2178
  0.7640
 -0.1874
  1.4368
  2.2441
  1.3668
  0.4773
  1.8339
  1.2443
  0.6740
 -0.3492
 -1.5895
 -1.8984
  0.8897
 -0.9720
 -0.2891
 -0.4612
  0.0030
  1.4528
  1.0625
  1.8746
 -1.1742
 -0.9876
  0.1371
 -0.1950
 -1.8867
  1.1538
  1.1323
  1.6759
  1.0617
  1.2721
  0.8928
  0.9407
  0.4111
  0.6902
  1.8910
 -1.3600
 -0.1351
  0.4249
  0.0382
 -1.5134
  1.0548
  0.7454
  3.2618
  1.0221
  2.7513
 -0.5169
  1.0427
  1.5677
  1.6678
  1.4944
  3.3600
  1.5300
  0.4453
 -0.5281
  0.6174
 -0.4043
 -0.1076
 -0.5568
 -0.1345
  0.9000
  1.9978
  1.8877
  0.5488
  0.3272
  0.4613
 -0.3478
  0.1570
  1.6937
  1.1044
 -0.0002
  0.9235
 -0.7490
  0.3771
 -1.3146
  1.2455
 -1.5046
  2.2909
  0.9513
  2.3667
  1.0984
 -0.4941
  0.9190
  1.7306
 -0.0516
  4.4401
  2.1946
  1.5131
 -0.8715
  1.6062
 -0.0268
 -1.2575
 -0.3329
  1.5630
 -0.0064
  1.1957
  0.1911
  1.2024
  1.6306
  0.6439
  0.2512
 -1.0523
  0.5242
 -1.6222
  0.7041
  0.6532
 -0.5390
  1.9918
  2.9666
  2.5118
  0.7174
 -0.6668
 -0.1014
 -1.6874
  0.0480
  0.9152
  0.1406
  1.4346
 -0.2149
  0.3575
  0.0890
 -0.7094
  0.4188
 -0.1304
  0.0076
 -0.9469
 -0.1530
 -0.7715
  0.5115
 -0.1325
  2.1636
  2.9926
  0.1230
  0.4221
  2.8876
  1.0511
  3.3164
  1.7048
 -0.4461
  0.0316
  0.3836
  0.8401
  2.1960
  1.0234
  1.8063
 -0.6292
 -0.7415
  0.8629
  1.7753
  0.6617
  1.5996
  0.4269
  0.2433
  0.2047
  1.4365
  2.0191
  0.6755
 -0.0425
  1.2027
 -0.1273
  1.0353
 -0.3781
  0.4261
 -0.3089
  1.1967
  0.7824
  1.9836
 -0.7753
 -0.4604
 -0.1510
  1.5760
 -0.1711
  0.9370
 -1.3320
  2.0204
  0.3022
  1.1019
  0.1153
  0.8180
  1.0990
  2.7568
 -0.8006
  0.1864
  0.5517
  0.6338
 -0.0585
  0.7583
 -0.0017
 -1.0536
 -0.5205
 -0.4189
 -0.2209
  1.2208
  0.1627
 -1.3688
  1.3897
 -0.6428
  1.0932
  0.6254
  1.6580
  1.3262
  1.3772
  0.9282
  0.1942
  0.6064
  0.7242
  2.1473
  1.3895
  0.3678
  1.0525
 -0.6881
 -2.7894
 -0.5750
 -0.1479
 -0.9499
  0.6247
 -1.2335
 -0.2724
 -1.3357
  1.1229
 -0.3712
 -0.4684
 -0.8002
 -1.4520
 -0.2542
  0.0488
 -0.4601
 -0.4184
  0.8004
  1.9043
 -0.5092
  2.6751
  1.0363
  0.2274
 -1.4556
  1.2790
  0.4551
  1.0885
  0.6514
  0.2688
  1.1995
  1.7829
 -0.2235
  1.5522
  1.2108
  3.0243
  0.2709
  0.7301
  1.1797
 -0.0954
  0.2970
  1.7896
 -0.1913
  0.5841
 -2.7657
  0.4908
  1.3581
  1.0281
 -0.0563
  2.2539
  0.9174
  0.6573
  3.0801
  0.1031
  0.4884
  0.5445
 -1.9115
 -0.1850
  1.8725
  0.6537
  0.5286
 -0.1594
 -1.0139
 -0.0061
 -0.8178
  0.9669
  2.0827
 -1.6668
 -1.9062
 -1.3588
  0.0260
  0.2962
 -0.3056
  0.8052
  1.2971
  0.3639
  0.4090
 -0.2982
 -0.0468
  1.7057
 -1.4868
  1.9505
 -0.3577
 -0.6409
 -0.4907
  0.2626
  0.3710
 -0.4382
  1.0190
 -0.0506
 -0.4877
 -1.0701
  0.5642
  0.2682
  0.1855
  1.3613
  1.1011
  1.7357
  2.3210
 -0.3684
  0.7600
  1.6722
 -0.3752
  0.8757
 -0.5094
  2.3266
  0.6030
 -0.8515
 -0.7793
 -3.2211
 -0.1824
 -0.5900
 -0.1338
  0.5582
 -1.1565
 -0.2567
  1.0833
 -0.6177
 -1.3749
 -0.5015
  1.1690
  0.8140
  2.2991
  1.0250
  0.0279
  0.7455
  0.7689
  0.4341
  3.2889
  1.3524
  1.1523
  2.5005
  1.2114
  1.6244
  1.5009
  0.0588
 -0.8149
 -0.4248
  0.6902
  0.8444
  0.8380
  2.0111
  2.8138
  1.0537
  0.8827
  2.0718
  1.7220
 -0.5200
  0.5189
  1.3410
 -1.1602
  1.1801
  1.4105
  0.5505
  1.3644
  1.8974
  1.4404
 -0.0970
  1.2724
 -0.7598
  0.1647
 -2.2841
 -0.6610
 -0.6911
  0.0749
  0.9112
 -1.1339
  1.3349
 -0.1183
  0.6901
  1.6237
 -0.1664
  2.1333
  0.7976
  0.8708
  0.2161
  0.0571
  1.3589
  0.5111
  1.3357
  0.4281
  0.1816
 -0.4653
 -0.6443
 -1.2091
  0.1528
 -0.0934
 -0.5689
  0.7998
  1.8559
  0.1220
  1.3348
  1.5850
  1.2195
  0.7485
  1.1645
  1.6043
 -0.7033
  1.9716
 -0.0424
 -1.2832
 -0.0574
  0.0365
 -0.3931
  0.8176
  2.1282
  0.7085
  0.9225
 -0.3443
  0.0332
  0.4301
  0.4805
  0.0425
 -2.2242
 -0.9604
  1.7421
 -0.2352
  0.2960
  0.3184
  1.0558
  1.6457
  0.2337
  0.7517
  1.2750
  1.6169
  0.8621
  0.5262
  1.3503
  1.0897
  1.0817
  2.1785
  0.2321
  2.8772
  0.0192
 -0.4153
  1.0901
  0.5495
  0.3799
  1.5853
  0.9260
  0.4108
  0.2835
  0.1896
  1.7734
  2.3366
  1.8186
  1.1572
  1.8940
  1.0422
  0.9210
 -0.1497
  1.3522
  0.4815
  2.2986
 -0.9936
  1.6641
  0.5914
  2.3888
  1.7281
  0.3885
 -1.2902
 -0.8782
  0.7943
 -0.2487
  0.4734
  0.3742
 -0.0824
  0.2778
  1.7815
  1.6461
  0.3267
 -0.1879
 -1.8832
  3.0882
  1.1096
 -0.3946
  2.3963
  0.9830
 -2.3146
  1.3793
  1.2974
  0.7935
 -0.4983
 -0.4137
  0.6885
  0.9305
  0.3721
  1.5782
  1.5220
 -1.4086
  1.2922
  0.9810
  1.1608
 -2.3885
  1.9169
 -1.0357
 -0.1611
  1.0502
  1.7068
  0.5709
 -0.4674
  0.2800
 -0.0752
  1.2325
 -1.0839
 -0.3671
  0.0724
  1.9515
  1.0028
 -2.1918
 -1.2788
  1.0865
  0.0401
 -0.3354
  1.8516
 -0.5235
 -1.4863
 -1.6584
 -0.0364
  0.7520
  0.3510
 -0.5676
 -0.2747
  1.3455
  1.0816
  2.8795
  0.8835
 -0.0429
  1.4551
  0.6868
  0.7515
  0.2867
  0.9822
  1.8629
  1.1910
  0.7695
  1.9186
  0.8060
  1.6458
  2.0261
  0.7631
  1.4771
  1.5884
  1.5212
  0.5812
  0.6958
  1.3410
  1.7492
  0.2721
  1.7281
  1.3066
  0.5833
  0.3349
  1.1156
  1.1886
  0.6670
 -0.0820
  0.5762
  0.8935
 -1.2258
 -0.3190
  3.7863
 -0.4198
 -0.4292
  0.1183
 -1.0840
  0.1262
 -0.7497
  0.6513
  0.8102
  0.1754
 -0.2277
 -0.8459
  2.0766
  2.2870
 -1.2361
  2.1226
  1.1913
  1.9824
 -0.0004
 -1.2296
  0.1025
  1.7716
 -0.3268
  0.1990
 -0.4014
  0.4381
  0.5549
  0.6081
  1.7061
  2.2267
  1.0508
  0.2338
  1.6972
  0.5219
  2.7996
 -0.0680
  2.9513
  2.0004
 -0.4316
  0.7135
 -1.7758
 -0.4888
 -0.6172
  1.3103
  0.7927
 -1.2837
  0.1729
  0.9543
  1.1855
 -1.2093
 -0.1172
  0.4642
  0.2236
  1.1088
 -0.5920
 -0.1908
 -0.8079
  0.4668
  0.5567
  0.4994
  1.5421
 -0.2166
  0.3372
  2.6025
  1.9031
  2.0614
 -0.6419
 -1.5278
  0.6383
  0.7621
  2.0710
  0.9336
 -0.7005
  1.1980
  2.8659
  1.3022
  1.4441
  1.3452
 -0.0216
  1.2511
  1.7528
  0.8956
  0.8512
  0.9531
  1.0993
  1.2333
  1.1010
  1.3712
 -1.9598
  0.7726
  1.2447
  0.8393
  2.4706
  0.9866
  0.2441
  0.5251
  2.6234
  1.3203
  3.3491
  1.5447
 -0.4841
  1.5821
  0.0192
 -0.8047
  1.3314
  0.0474
  0.8125
 -1.2863
  0.5215
  0.0320
 -0.6533
 -0.9074
  1.5158
  0.3815
  2.5332
  1.4535
  0.7420
  0.7800
  1.6806
 -0.2519
  0.6382
  2.0394
  1.7616
  2.1692
 -0.5027
  0.9990
  0.5710
  1.0651
  0.3851
 -1.0046
  2.0635
  2.6854
  1.4879
 -0.1385
  2.1050
  1.1287
 -0.6616
 -1.2135
 -0.3417
 -1.0172
 -0.8402
  0.4368
  1.0346
  0.9500
 -1.0121
 -0.1565
 -0.1323
 -1.4275
 -1.1326
 -0.1842
  0.7839
  0.4708
 -1.3712
  0.3523
  0.4341
 -0.2065
 -0.3859
 -0.5952
  0.6778
  3.2339
  2.3512
  1.0842
  2.0308
  3.6442
  0.4282
  1.7301
 -0.0553
  1.1710
  0.1795
  1.4604
  0.6592
  3.4478
  0.2570
  0.4643
  0.0890
  0.8684
 -0.1901
  1.3339
  1.0496
  0.1713
  1.5267
  1.8457
  1.9216
  1.5291
  0.6634
  2.0147
 -0.9941
  1.6584
  2.0491
  0.6474
  0.6402
  3.0794
  0.2901
  0.4620
  0.1393
 -0.3191
  1.3027
  3.4875
  1.7753
  1.2702
  1.6232
  0.9136
 -0.4113
  0.2431
  0.6765
  0.5759
  0.5407
  0.8429
  0.9094
 -2.1886
 -1.2068
  0.7136
 -0.5027
  0.3724
  1.2284
  2.6804
  1.3702
  1.5687
  0.9135
  1.6889
  0.9292
  1.0500
  0.1328
  0.4409
  2.0144
  0.3396
  1.6945
  0.9443
  1.4440
  1.6531
  1.8528
  1.1572
  0.3751
  0.3860
  0.8055
  2.1444
  0.1400
 -0.2289
 -0.7430
 -0.2682
  0.0525
 -1.4218
  0.0066
  0.0571
  0.6942
  0.2853
  1.1703
 -0.0207
  1.3627
  0.5452
  0.1747
 -0.1152
  1.3812
 -1.8736
 -1.0383
  0.8498
  1.0340
  0.7720
 -0.4168
  1.8950
  0.4522
  0.8967
 -1.6075
  0.5991
  3.2610
  0.6050
 -0.5432
  2.2666
 -0.2481
  0.8778
  1.1720
  1.9950
  1.4999
  1.9570
  1.7372
  1.1379
  1.4841
  1.5370
  2.2771
  1.3651
  2.7149
  3.3138
  0.4595
 -0.5261
  0.8023
  0.6659
  1.2732
  1.0267
  2.0685
  1.3361
  0.4811
  2.6551
 -0.1192
  1.0592
  0.7049
 -0.2187
  0.3276
  0.5516
 -0.1413
 -0.7494
 -0.2264
 -2.0652
 -1.3686
  0.1429
  1.4793
  0.1897
 -0.5325
 -0.6685
  0.0050
 -0.7011
 -0.8830
 -0.2675
 -1.8324
  0.4171
  0.9417
  0.0845
  0.6039
  0.1434
  1.0511
 -0.0752
  0.7905
  0.4915
  1.0144
  0.9657
  1.1018
  0.0269
  2.5989
 -0.6277
  1.7265
  1.1487
  0.0517
 -0.0703
  1.8860
  1.0951
  0.4332
 -1.0435
  0.7140
  0.4565
 -0.2244
  2.0075
 -0.2359
 -0.4162
 -0.7105
  0.1378
  0.3204
 -0.9549
 -0.0858
  1.3750
  0.9824
 -0.3351
  0.6658
  0.4390
 -0.0863
  0.7827
  0.2537
  0.4740
  0.2182
  1.0182
  0.4427
  1.1221
  1.4740
  2.5008
  1.7944
 -0.2333
  0.6866
  2.0270
  0.3089
  1.3608
  2.3866
 -0.2438
  0.1971
 -0.5058
  0.3899
  0.5433
  0.8134
 -0.0525
  1.7833
  1.7376
 -0.8452
  1.4962
  0.9644
  1.5718
  1.6167
  0.2944
  0.0109
  0.5512
 -0.3950
  0.3031
  0.6335
  0.2334
  1.8005
  0.8763
  0.8361
  1.8298
  2.0040
  1.0994
  1.0309
 -0.7542
  0.6563
  0.4631
 -0.1414
 -0.4196
  0.2835
  0.6250
 -1.6866
  0.1072
 -0.5463
  0.3700
 -0.0185
  0.9919
 -2.5904
 -0.0009
 -0.2877
 -0.5503
  0.9250
  3.2276
  0.9877
  0.0361
 -0.0056
  0.1108
  1.3339
  2.1200
  1.3322
  0.6807
  2.3149
 -0.7764
  0.9585
  0.8585
 -0.6734
  0.7745
  0.8272
  0.4576
  0.2478
  0.0518
  1.8158
  1.0024
  1.9464
  0.3796
  2.9197
  0.0576
 -0.4900
  1.6623
 -0.2566
  0.4709
 -0.3322
  0.1240
  1.0207
  0.5053
 -1.0574
  0.7583
  1.7073
 -0.1035
 -0.8813
  2.2749
 -1.0549
 -1.4977
 -0.8816
  0.4849
  3.1706
  0.5961
  1.3673
  1.4904
  1.3146
 -0.5035
  1.4709
  1.2629
  1.9915
  1.8432
  0.4720
  0.2353
  1.6675
  2.5839
  0.8535
  0.1829
  1.1226
  0.0063
  1.9017
  1.3747
  0.1701
  1.6809
  0.8399
  2.0449
  2.2485
  1.2806
  0.3391
  1.3753
 -0.3536
  2.6588
  1.0936
  1.2452
 -0.5196
  1.0827
 -0.9309
  0.2306
  0.8946
  2.3972
 -0.1158
  0.9329
 -0.6203
 -0.3148
  1.6808
 -1.4143
 -1.4222
  0.3249
  0.2631
 -0.6222
 -0.4702
 -3.5494
  0.3412
 -2.3122
 -1.3876
 -0.0688
  0.5056
  0.3530
  1.5961
  2.1177
  0.7779
  1.0348
  3.0116
  1.1190
  0.5213
  0.0213
  2.5338
  2.0465
 -0.5106
  0.1144
  0.8052
  0.4784
  0.0383
  0.7065
  2.4325
  1.2193
  2.7077
 -0.0487
  2.2022
  2.5109
  2.5820
 -0.8702
 -0.9622
 -0.8462
 -1.0475
 -0.4589
 -1.9730
  0.0298
  1.5630
  1.1815
  2.0939
 -0.5682
  1.0880
  1.3542
 -0.2310
 -0.5148
  0.1639
  1.2777
 -1.6706
  0.3181
  0.4550
  0.4789
  0.1088
  0.5215
 -0.1674
  0.9396
 -0.8673
  1.6550
  0.7012
  1.3521
 -0.4930
 -0.2019
  0.2257
  1.1720
  0.2476
  0.6666
  0.0710
  2.0408
  1.3418
  1.3916
  0.7812
  1.2343
  2.8982
 -0.2333
  1.0513
  2.7059
  1.6243
  1.3176
  2.4060
  0.1900
  0.1159
  1.9428
  1.9728
  1.5610
  2.0003
  0.1294
 -0.3393
  1.9090
  0.2867
  0.6402
  1.5730
  0.9015
  2.0122
  0.9615
  0.8990
  1.2529
 -0.4117
  0.6138
  1.5923
 -0.1314
  0.6463
 -0.0023
 -1.5953
  0.9388
  1.3362
 -1.2509
 -1.2966
  0.9987
 -1.4236
  1.5025
  0.2285
  1.7060
  0.7397
  0.5487
  2.5865
  0.8985
  0.1689
  1.2064
 -0.2008
 -1.3840
  0.5774
 -0.6047
  0.6939
  0.8584
 -0.1160
  1.0666
 -0.6450
  2.3351
  1.4052
  0.7294
  0.4218
 -1.3074
  0.0418
 -0.1151
 -0.0770
 -0.4368
 -0.6830
 -0.1251
  0.1572
  0.2400
 -0.5977
 -1.2113
  1.3461
  0.0661
  0.0230
  1.4549
 -0.0066
 -0.4669
 -1.4050
  1.1523
  0.8618
  2.0182
  0.1422
  1.1152
  2.0396
  4.6407
  1.2633
  1.5076
  0.6764
 -1.2458
  2.6644
  0.7665
  1.2845
  1.6450
  1.8146
 -0.2091
  0.9499
 -0.2116
 -0.5028
  0.8232
  1.7439
  1.8129
 -0.1883
  2.3021
  0.9942
  1.7011
  1.0275
  0.2516
  1.4165
  1.4362
  0.1290
  2.7314
  1.3898
 -0.9139
  1.2403
  0.3527
  0.9721
  1.5108
  0.6222
 -1.3403
 -0.0692
 -0.5027
  0.2293
  0.3919
  2.2291
  0.1828
  1.2516
  1.6865
 -0.6915
  0.9904
  0.8705
 -0.7929
  0.6878
  0.7900
  0.2593
 -0.5046
  0.7971
  1.2168
 -0.5397
  1.5547
  1.3212
  1.6769
  0.4280
  2.5956
  2.0242
  2.7758
  0.6524
  0.8592
  0.4494
  0.9375
 -0.3721
 -0.2845
 -1.6431
  0.9320
  0.6636
  1.4701
 -0.1657
  1.0496
  0.6044
  1.0227
  1.9521
 -0.3599
 -0.4871
  0.1007
 -1.3416
 -0.7875
 -0.3512
  0.4130
 -0.6369
  0.4448
  0.3808
  2.1075
  0.6384
 -0.9145
  1.0990
  0.6650
 -0.3651
  0.8883
  1.0667
 -0.1023
 -1.5665
 -1.2307
 -2.0404
  1.8122
  1.8724
  1.6846
  3.1811
 -0.0000
  0.6721
  1.0686
 -0.1130
  0.6590
  0.0611
  2.1051
  0.5119
  1.5820
 -0.0040
  1.8034
  1.8937
 -0.2167
  0.6249
  0.4276
  0.6524
  0.7446
  0.2367
  1.5882
  2.0563
  1.9290
  2.4743
  0.2387
  0.2831
  1.1262
  0.5904
  0.9077
  0.3826
  0.4572
  0.8393
  1.7457
  1.0864
  2.1395
  0.7033
  2.4183
 -1.8902
 -0.7918
  0.4299
  0.5661
 -1.2100
  0.0911
  0.0230
 -1.0184
 -1.5205
 -0.4599
 -2.0747
  0.4634
 -0.8784
  0.4403
  0.8951
  1.0019
  0.3251
  2.2364
  1.8741
  0.6378
  1.1925
 -0.2415
  2.4081
  0.6704
  2.3906
  0.1905
  0.4172
 -0.0434
 -0.7327
  0.2484
 -1.6006
  0.3018
 -0.5085
 -0.4849
  0.1822
 -0.2493
  1.7922
  1.3316
  1.8078
 -1.4712
  0.5782
  1.5760
  0.9076
  0.2369
 -1.7386
 -0.4680
 -1.0827
  0.9846
 -0.9702
  1.0352
 -0.7119
 -1.3816
 -1.1210
 -1.0350
  0.4667
  0.0669
 -1.3185
  0.8753
 -0.7780
  0.1769
  1.0233
  2.5029
 -0.7246
  1.1648
  2.6675
  0.4026
  0.9266
  1.8093
  0.0014
  0.1190
  0.3722
  2.9085
  2.6524
  2.0858
  0.9074
  0.0700
  1.4690
  0.7560
  2.5552
  2.0674
  0.8338
  2.3582
  1.9860
  0.1429
  2.2071
  0.8133
  0.6019
  0.8096
  3.4619
  1.4758
 -0.0483
 -0.2810
  1.9345
  0.0434
 -0.1922
  2.2562
  0.6638
  0.9665
  1.8556
  1.0036
  1.7423
 -0.4444
 -0.8459
  0.4620
 -0.1715
  0.6764
 -1.5345
 -0.3818
  0.4261
  1.1936
  0.6490
 -0.2063
  1.0656
 -0.4593
  0.3776
 -0.0002
  2.7203
  1.0095
  1.3506
  2.0026
  0.4156
  1.4430
 -0.4150
  1.9555
  0.0726
 -0.3243
  0.6048
  1.3693
  1.2346
  0.0627
 -1.0375
 -0.1501
 -0.4348
  0.4875
  1.2001
  0.4261
  1.7066
  1.2760
 -0.0407
 -0.1241
  0.2665
  0.2964
 -0.8825
  0.7213
  0.9371
  0.0370
  0.0658
  0.2492
 -0.1146
 -1.7164
  0.0840
 -0.6528
 -0.6584
  0.6804
  0.8775
  2.9538
 -0.9064
 -0.9397
 -0.6271
  0.4785
  0.5285
  0.9708
 -0.1876
  1.7643
  2.2496
 -0.5079
  1.8836
  0.8957
  1.5167
  1.7605
 -0.2344
  2.5421
 -0.4228
  0.7211
  1.6948
  0.1445
 -1.2011
  0.3298
  1.0876
 -0.5933
 -1.1138
 -0.7014
  1.7402
  0.5912
  0.4109
  0.3522
  1.8081
 -0.9342
 -0.6959
  1.4430
  1.8220
 -0.6833
  1.9929
  0.0383
  0.7828
  2.2689
  1.6437
  3.5126
  2.4119
  0.2112
 -0.0266
  0.7680
  0.4982
 -2.6158
 -0.5814
 -0.4984
 -0.3141
 -1.2593
 -1.2291
  1.2801
  0.5520
  0.7976
  1.3704
  0.4787
 -0.3668
 -0.1948
  1.0122
  1.4441
  0.4957
  1.5920
 -0.3108
  0.1194
  1.3083
 -1.0519
  0.6383
  1.2136
  1.1385
  1.6649
  0.6027
  1.0080
 -0.5883
  0.5122
 -1.1255
 -0.4182
 -0.4221
  2.2142
  1.9086
  1.1057
  0.8635
  0.4206
  1.5724
 -0.5363
  0.5679
  0.9200
 -0.5785
 -0.4679
  1.5679
 -1.3381
 -0.9734
  1.9041
 -1.7327
  0.0288
  0.7421
 -1.9284
  1.2027
 -1.7364
 -2.3182
  0.3145
 -2.5309
  0.1126
  0.2081
  0.0700
  0.4408
  0.3099
 -0.9952
 -1.1449
  1.7547
  0.6173
 -1.3291
 -0.8441
  0.6958
 -1.2794
  1.3024
  1.2929
  0.5182
  1.3802
  0.5259
  0.7617
  0.2457
 -1.2972
  0.4082
  3.1421
  2.3131
  0.9631
  1.7340
 -0.6498
  0.2277
  0.2243
 -0.2437
  0.8870
  0.8880
  2.3831
  0.6903
  1.3201
  1.3201
  1.0680
  0.9592
  2.3162
  1.7104
  1.5772
  2.3729
 -0.2852
  1.6135
  0.3485
  1.9241
 -0.0653
  1.0096
  0.4342
 -0.9608
  0.0231
  0.0852
 -0.1963
  1.2215
  0.4070
  1.5712
  0.0187
  1.2892
 -1.1446
 -0.7839
  0.1427
  0.6063
  1.8168
 -0.0998
  0.5051
  1.2812
 -0.1967
  1.9741
  1.5189
  0.3622
 -1.1000
  0.7341
 -0.6648
  0.7692
 -0.2280
 -0.5938
  1.5729
  0.8851
  2.8025
  0.3572
  0.3573
  1.1435
  1.6392
  1.8732
  1.1488
 -2.3921
 -0.0590
  0.6575
 -1.1814
 -0.8992
  1.0091
  1.0477
  0.1585
  1.0646
 -0.7375
 -0.5487
  0.3127
 -0.5410
  1.9615
  2.1235
  0.0733
  0.7208
  0.3454
  0.2319
 -0.8692
 -0.4478
 -0.1794
  0.9802
 -0.5258
  0.1228
  0.8461
  0.2495
  0.5634
  1.5114
  0.2150
  0.9233
  1.0723
  1.1624
  1.1008
  1.8826
  0.4425
  0.5620
  2.9467
  2.6389
  1.3155
  1.4393
  1.5799
 -0.0019
  0.2298
  0.6321
 -0.8307
 -1.4008
  0.7949
  0.9641
 -1.0295
 -0.5970
  3.0968
 -0.0168
  1.3568
  1.4122
 -0.3359
  0.8666
  0.6785
 -0.7330
 -1.3601
 -0.0147
 -0.3118
 -1.3000
  1.4497
  0.9353
 -0.1869
 -0.9714
  0.9998
  1.6483
  0.2433
 -1.1139
 -2.5569
 -0.5913
  1.7977
  1.5754
  1.0996
  0.8350
  1.1932
  1.2901
  0.3861
  1.0734
 -0.4655
 -0.0264
  0.7177
 -1.3060
  0.1460
  1.7198
  0.1870
 -2.7735
 -1.2985
  1.0230
  0.4165
  2.1750
  0.1954
  2.4973
 -0.3486
  0.4482
  1.6480
  0.8160
  0.0919
 -0.4254
  1.0955
 -0.6492
  0.2026
 -1.4999
  1.4384
  0.4681
 -0.9257
 -1.7147
 -2.5246
  0.0926
 -0.2452
 -0.4378
 -0.9803
 -1.1820
 -0.1028
 -0.6195
 -0.3900
 -0.8174
 -1.2883
 -0.3594
  0.0748
 -0.3108
 -1.7206
  0.1868
 -0.3439
  1.4403
  1.9461
  1.9956
  2.1796
  0.9509
  0.6761
  1.8649
  1.5563
 -0.5308
  1.0781
 -0.1877
  0.3570
 -0.1201
 -0.0941
  0.5303
 -0.7774
  1.9826
  2.1053
  1.0106
 -1.4758
 -0.5306
  0.3133
 -0.6305
 -1.7655
 -0.5959
 -1.3359
  1.2424
  1.6090
  3.2447
  1.3255
 -1.0877
  2.0900
  0.8215
 -0.0381
  1.6115
 -0.1425
  1.1148
  0.5374
 -0.2688
  2.0551
 -1.0620
 -0.3746
 -1.6411
 -1.7673
  2.8721
  1.1128
 -0.0151
  1.2690
  2.0674
  1.1388
  0.8877
  1.8631
  0.8026
  1.7818
  0.4553
 -0.5434
 -0.1856
  0.4484
  1.4222
 -1.1165
 -0.8478
 -1.3464
  0.0426
  1.2164
  1.3482
 -0.0287
  1.7752
  1.0850
  1.8104
  1.2963
  0.8468
  2.1019
  0.6026
  0.4922
  1.6687
 -0.9623
  0.4821
  1.3928
  0.6859
 -0.6155
 -1.4870
  0.8264
  1.1843
  2.2724
 -0.1313
 -0.5798
 -0.5305
  1.5231
  0.4801
 -0.0530
  0.4952
  1.0021
  0.4519
  0.9401
  1.1346
  1.9443
  0.4354
 -1.3614
 -0.2636
  1.7986
  0.5548
  1.4104
  1.7877
 -0.4356
  1.5723
  1.1481
  0.0287
  0.9538
  2.1557
  0.6300
  1.1442
  1.3884
  2.6236
  0.2918
  0.4192
  0.9605
  1.0108
  0.0473
 -0.3595
  1.2057
  0.9512
  0.2783
  1.3469
 -0.7479
  2.6885
  1.7795
  1.5511
  1.7046
  0.2699
  1.4941
  0.1209
  2.0815
  2.2881
  0.2995
  0.4977
  0.6139
  0.3889
 -1.5647
  0.2611
  1.1646
 -1.2061
  0.8750
  1.8458
 -1.2247
 -1.8595
  0.8801
 -0.4823
  0.9153
 -0.6160
 -0.6431
 -0.0773
  0.1348
  1.7152
  1.6877
  1.0753
  3.3871
  0.6378
  1.7731
  1.8698
  0.3064
  0.2350
  1.1651
  0.0872
  1.1022
  2.1421
  0.9280
 -0.2845
 -0.1643
  1.5056
  1.4655
  1.5593
  1.0146
  1.0397
 -0.8607
  0.9705
  1.6921
  0.0492
  3.5426
  0.6071
  2.0895
  2.0954
 -0.5580
  0.6715
  0.6434
  0.8130
  0.2129
  0.1331
 -0.6483
  0.6460
  0.7078
 -0.9100
  0.8729
  1.3329
  1.1246
 -1.6368
  0.2240
 -0.8539
  1.0340
  0.0059
  0.0645
  0.0000
  0.5199
 -0.0590
  0.4647
  0.0173
  4.1378
  0.1355
  1.3654
  2.4419
 -0.0460
  1.1231
  2.3437
  1.2047
  2.6619
 -1.5441
  0.4802
 -0.4931
  1.7857
 -0.2184
 -1.2108
  1.8622
  1.0772
 -0.9211
 -0.9186
 -1.1765
 -1.6540
  1.3813
  0.2931
  1.7008
  0.9493
  1.2573
  0.5898
  1.4629
  1.0551
  0.5578
  2.2541
  0.9273
  0.6036
  1.6944
  0.9055
  2.4639
  0.0369
 -1.0895
  1.0692
  0.9706
  0.6696
  0.6041
  0.0126
 -0.1171
  2.0664
  0.2335
 -0.8796
 -1.0240
 -0.6775
  0.8414
  1.2477
 -0.2663
 -1.2557
  0.2164
 -1.3471
 -0.0973
 -0.2849
 -0.2260
 -0.4460
 -0.4107
  1.8497
  1.1919
  0.8568
 -1.0824
  1.1103
  1.5596
  1.6055
 -0.0707
  1.3027
  0.7248
  1.4655
  0.6516
  3.1232
  1.5948
  1.2866
  0.3299
 -0.5539
 -0.3372
 -0.9842
 -0.9596
 -0.8048
 -0.4622
 -1.5663
 -0.4434
 -1.4678
  0.2319
  1.6665
 -1.8161
  1.9452
  0.2318
 -1.1870
  0.0593
 -0.5274
 -1.4261
  0.7568
  1.6505
  0.1803
  0.1272
  0.8147
  0.7327
  1.8867
  0.6276
 -0.1998
 -0.0282
  0.2524
  0.6717
  1.1783
  0.1911
  0.2504
  2.0430
  3.1902
 -0.0776
 -0.5735
  0.5612
  1.7963
  0.3605
 -0.7191
  0.5289
  1.0729
  1.2431
  0.4928
  0.5446
 -1.0114
  0.9117
  0.1206
  2.7967
  0.2064
  0.7378
  1.2409
  0.4278
  0.5523
  0.4669
 -0.5119
  1.0650
  1.1846
  2.5311
 -1.5872
  1.2682
  0.0577
 -0.0576
 -0.5702
  0.8217
  0.1534
 -0.1742
  0.8708
  1.1023
 -0.4637
 -0.1919
  0.7822
  2.0999
  2.4374
  0.9372
 -0.8196
 -0.4238
 -0.5186
  0.9497
 -0.0764
  2.1818
  0.0152
 -0.4447
  1.2797
  2.1202
  1.2519
  0.2575
  1.5332
  1.7261
  0.9660
 -0.1214
  0.9181
  0.1246
 -0.1602
 -0.7435
  1.7189
 -0.1178
  2.0524
 -0.2301
 -1.0207
  0.8019
  1.1826
  0.1397
 -0.3503
 -0.7584
 -0.7281
 -0.2079
 -0.0163
  0.0857
  0.1085
 -0.0733
  0.2882
  0.4932
  1.5389
 -0.5286
 -0.6899
  2.1208
  1.8397
  1.8054
  0.3421
  0.8193
  1.7020
 -0.5194
  0.3327
  0.8111
  2.0562
  0.8626
  1.4978
  0.6971
  2.6167
  2.5756
  2.6988
  1.3771
  0.3514
  0.6649
  0.3809
  0.1867
 -0.0770
 -0.8209
 -0.3938
 -0.5371
 -1.1823
 -2.3123
  0.1849
 -0.2595
  2.0956
  0.4296
  2.6756
  1.5995
 -0.9034
  0.1934
  0.5988
  0.6414
  0.9034
  2.3871
  1.8372
  1.0032
  1.6225
  2.5164
  0.1486
 -0.6727
 -0.4448
  0.8392
 -0.4444
  0.7762
  0.7419
  0.3942
  1.1039
 -2.0690
  0.2659
  0.9499
 -0.2231
  0.1007
  0.6660
 -1.1373
  0.3872
  0.2729
  2.0304
  0.7941
  1.5400
 -0.3539
  0.3569
 -0.6725
  1.9666
  2.9626
 -0.3704
  0.5087
  1.3855
  2.2198
  2.4327
  1.1795
 -0.6483
  1.7669
  0.2360
  0.7926
  3.4981
  0.5938
  1.9737
  0.1893
 -0.9985
  0.3177
 -0.2497
  1.5894
  0.3598
 -0.0422
 -0.0691
  1.5256
  0.6704
 -2.1057
 -0.4446
 -0.5403
 -1.2415
 -0.2026
  0.5803
  0.2313
  0.0523
  2.7597
 -0.3004
 -0.0333
  0.2566
  1.2117
  1.5982
  1.7196
  1.2698
  1.5682
  1.7852
  0.5070
  0.5743
  0.4259
  0.9208
  1.0892
 -0.8820
  1.0056
 -0.7497
  1.5424
  0.4067
  0.6596
  0.0712
 -1.7007
 -0.3741
 -1.3537
  0.4741
 -0.0667
  1.1731
  1.6158
  0.7193
  2.5951
  1.4670
  0.1722
 -0.1117
  1.0107
  0.2642
  0.0127
  1.1618
  1.7342
 -0.4463
  2.6419
  0.1945
  0.8419
 -1.6056
  0.0008
  0.3622
  0.6107
  0.3643
  1.8802
  0.6026
 -1.0512
  0.4012
 -1.7121
  0.2198
  0.6024
  1.7914
  3.7950
  0.9236
  1.7788
 -0.1673
 -1.6304
 -0.3523
  0.0050
 -0.3969
  1.5949
  0.3043
  1.7296
  1.4067
  1.7941
  2.4476
  1.9092
  0.2692
  2.2265
  0.7383
  0.9911
  0.8484
  1.7708
  0.2177
  1.9507
 -0.3868
  0.9485
  0.3335
 -0.1776
 -0.0875
 -0.1719
 -1.0789
 -0.0165
  0.7013
 -0.9552
  0.3222
 -0.4210
  1.2034
  1.3339
  0.0028
 -0.6229
 -0.1997
  0.3628
  0.2280
 -0.5728
  0.4162
  0.4965
  2.7331
 -0.1605
  1.2052
  1.4011
  0.7328
  0.9220
  1.6832
  0.8656
  1.8387
  1.2035
  0.8110
  0.4278
  1.6967
  1.6244
 -1.2180
  2.7646
  0.0860
 -0.3358
  3.0721
 -0.8138
  0.5213
 -0.0307
 -1.2004
  1.2536
  0.6195
  3.1085
 -0.8137
  0.9209
  2.0581
  0.3110
 -0.0149
 -0.3966
 -0.5084
 -0.1236
 -0.0713
 -0.2925
 -1.1347
  0.6223
  1.6695
 -0.0342
  1.8418
  0.8677
  2.3857
  0.4346
  1.5419
 -0.4648
  1.2683
  2.1043
 -0.3357
  0.6484
 -2.5766
 -0.8915
  1.8754
  0.4659
 -0.1023
  0.1483
 -0.0119
 -0.8020
  2.0537
 -1.4150
 -0.6941
 -0.1020
  0.2264
 -0.3382
  0.6309
  1.5824
  0.9447
 -0.1324
  0.7557
  0.4823
  1.4618
 -0.6722
  1.4116
 -0.1434
  1.1562
  0.7400
  1.4934
  1.3621
 -0.3221
  2.1511
  1.3660
 -0.2240
  0.9053
  0.0570
  0.4776
  1.0100
  0.5017
 -1.0324
 -0.4912
  0.8546
 -0.6501
 -0.1693
 -1.0255
  0.6689
 -1.1850
 -1.2595
  0.2101
 -0.5356
  0.6528
 -0.6043
  1.5506
  2.6271
 -0.3287
  0.4624
  1.9280
  1.7394
 -0.8694
  0.2120
  0.6048
  0.3478
 -0.4117
  1.3762
  0.9907
 -0.8272
  0.3073
  0.9041
  0.0009
 -0.4664
 -0.2371
  0.1303
 -1.0031
 -0.1272
  0.9600
 -1.2320
  1.3953
 -0.1469
 -0.2911
  2.4423
  0.4116
  0.9900
  3.3576
  0.5975
 -0.6889
  1.6307
  1.3967
  1.7393
  0.1622
  0.7890
  1.0319
  0.9039
  1.9617
 -0.2844
 -0.3611
  0.5175
 -0.4140
 -1.7529
  0.6692
 -1.9781
  0.7128
  1.1545
 -0.3046
  0.5491
 -0.3386
 -0.5100
  0.1444
 -0.5793
  1.3093
 -1.0421
 -0.0605
  0.8795
  0.3433
 -0.9889
 -1.3185
  0.6620
  0.0712
  1.5584
  1.9963
  0.4404
  1.3362
  3.6901
 -0.3199
 -0.1192
  0.4426
  2.0164
  1.1822
  1.9496
  0.4120
  1.4882
  0.4757
  1.1174
  0.6938
  0.3518
  0.3796
  1.3098
 -0.8243
 -1.9243
  2.2994
 -0.0104
 -0.2878
 -0.7197
  0.4761
 -0.9909
  0.3758
 -2.0988
 -1.6809
  0.9225
  0.0144
 -0.6099
  0.3386
  1.6657
  1.4473
  2.7297
  2.1442
  0.7553
  0.6745
  2.6484
  0.4621
  1.0683
  0.1898
  0.4845
 -2.1305
  0.6315
 -0.1220
  0.2578
 -1.7022
 -2.0767
 -0.3568
 -1.6947
 -1.1446
  0.0942
  0.0470
  2.2358
  1.4750
  0.8525
  1.3782
  1.7242
  1.6467
  1.0643
  0.5185
  1.9038
  0.8940
  0.0351
  2.3246
  1.5556
  0.8035
  1.6692
  0.4396
  4.4769
  1.7888
 -0.1205
  1.6984
  2.0698
  0.8837
  0.7553
 -1.9890
  0.6949
  0.2295
 -0.0593
  0.9122
 -0.3009
  0.0198
  1.4543
  0.5160
 -0.9350
 -1.0331
 -1.0056
 -0.1024
  0.6146
 -0.8663
 -1.1514
 -2.0480
 -1.3691
 -0.9955
  0.4927
  1.4107
 -0.9475
  2.4940
  2.2275
 -0.6705
  0.5120
  1.4465
 -1.8287
 -2.4410
 -0.2625
  0.7461
  1.5526
  0.2513
  1.5087
  0.4185
  0.1956
 -0.1789
  0.8815
  0.2889
 -1.1592
 -1.1133
  0.2491
 -0.3528
 -0.5482
  0.0873
 -1.0579
  0.6555
 -0.4569
  1.1158
 -0.9699
 -0.3535
  0.1991
 -0.7367
 -0.1986
  2.5957
 -0.0552
  0.9662
  1.6371
  0.0452
  0.5871
  0.8718
  1.3577
  1.8044
  0.0331
 -0.8045
  2.6733
  0.5250
 -1.3147
 -0.6681
  0.5136
 -0.0424
  2.9509
 -0.0089
  1.3850
  1.1341
  1.5861
 -1.0464
 -0.2161
  0.8173
  0.0513
  0.1353
  1.6475
 -0.6503
  1.6460
  1.7810
  0.3489
  0.6847
  1.5219
  1.1339
 -0.7325
  0.5175
  0.5016
  2.4935
  0.1017
  1.1055
  0.4356
 -0.3227
  1.0059
 -1.1049
  1.4725
  1.9796
  1.6012
  2.2939
  0.3773
  0.7816
  0.3876
  0.0905
  0.3030
  0.6097
 -0.0799
 -1.0714
  1.3108
  1.2692
 -1.4491
 -1.4146
  0.5388
  0.8587
 -0.1556
 -0.0644
  0.2185
 -0.2026
 -1.0868
 -1.0541
  0.6750
  0.0913
  2.8960
  0.6120
 -0.4947
 -0.8764
 -2.3375
  0.0821
  0.8625
  1.2608
  0.7078
  1.4241
  1.1452
  1.9344
  0.6502
  0.9904
 -0.1363
 -1.2166
 -0.2216
 -0.6417
  0.7381
 -1.1045
  0.3344
  1.5214
  0.8266
 -0.6975
 -0.7025
  1.4408
  0.2082
  1.2519
  0.2243
  1.0953
 -1.4436
  1.5251
  1.5502
  0.1758
  0.7644
  1.0260
  0.1071
  1.5152
  0.3888
  0.9749
  3.7631
  1.3041
  0.1994
  1.5945
 -1.9200
 -0.1164
  0.5832
  0.1273
  1.1672
  0.7392
  0.6612
  1.0506
  0.7579
  0.6627
  0.0111
  1.0952
  1.2410
  1.9064
  1.2124
  1.3264
  0.4543
  0.2491
  0.5180
  0.9939
 -0.9217
  0.6632
  2.1703
  1.8909
  1.4457
  1.2656
  0.6335
  0.4694
  2.7396
  1.3135
 -0.6177
  1.3653
  1.2780
  0.2040
  0.8840
  0.7988
  0.5120
  0.5859
 -0.3259
 -0.0108
  0.2815
  1.1156
  1.3802
 -1.0876
  1.6615
  1.1539
 -0.2291
  0.1744
 -1.0758
  0.1539
  0.3207
  1.4800
  2.2796
 -1.7111
 -1.0085
  1.8129
 -1.4999
 -0.4094
  0.2356
 -1.2424
 -0.5863
 -0.3939
  2.2806
  1.4341
  0.3404
  0.0422
  0.3326
  0.8412
  1.7547
  2.0467
  3.2925
  0.2535
  2.1943
  0.0224
 -1.0217
 -1.5855
 -1.3602
  0.0285
  0.1496
  1.3903
  1.1353
 -0.0269
  1.0484
  0.4365
 -0.3158
 -0.1675
 -1.8982
  3.0886
  0.8535
  1.0929
 -0.1776
 -0.4311
 -1.2962
  1.2198
  0.1717
  0.3404
  1.8743
  0.2641
 -0.6729
  3.1209
 -1.3167
  0.2565
  0.8039
 -0.1869
 -0.9401
  1.1344
  0.8696
  1.1057
  2.1742
  0.5348
  0.5264
  1.9504
 -1.7199
  1.9847
  0.5998
  2.6173
 -0.0467
  2.0665
  1.6796
  1.5095
  1.0173
  1.1920
  3.2402
  1.0237
  2.2328
  1.6592
  0.6980
  1.5845
  0.5783
  2.6010
  0.0029
  0.7022
 -0.2937
  1.5400
  0.2525
  1.8153
 -0.8929
 -0.4243
 -0.1979
 -0.3276
 -2.4956
  1.1456
  1.3251
  0.0393
  0.3283
  0.5964
 -1.5539
 -1.5245
  1.1651
 -0.3873
  0.4644
 -0.2470
 -1.8322
 -0.6315
  0.9814
  0.3860
 -0.1852
 -0.0143
 -2.3180
  1.2335
  1.4970
 -0.9655
 -0.2138
  0.4337
  1.0913
 -0.0093
  0.1433
 -0.1032
  0.8996
 -0.1951
  0.4455
  0.1065
  0.7426
  1.4072
 -0.2894
 -1.0441
 -1.6604
  0.4926
  0.1576
  1.3407
  1.2800
 -0.9168
  1.7486
 -1.7488
 -0.0182
  1.3470
  1.0349
  0.1255
  1.6362
  1.4296
 -0.1301
  1.0420
  0.7869
  0.8465
  0.8988
  1.6134
  1.0627
  1.1829
  0.5421
  0.1862
  1.5098
  0.7507
  2.1039
 -0.2606
  0.3153
  2.2701
  0.4822
  0.5783
  1.3357
  0.0330
  2.2352
 -1.5364
  1.1346
 -0.7312
  0.5404
  1.6110
  0.6950
 -0.0182
  2.0260
  1.0283
  1.8420
  0.2620
  2.4129
  1.3451
 -0.0781
  1.2677
  1.8489
  2.3972
  2.3908
  1.5708
  3.0520
  2.5403
  0.7712
  0.0461
  0.8669
  0.1310
 -0.4093
  0.5958
 -1.4363
 -0.4410
 -0.6417
 -0.4436
 -2.0495
  0.5042
  0.0570
  0.6244
 -1.4279
 -0.8186
 -1.2651
  0.2421
  0.1598
  0.3525
  0.6485
  1.3274
  2.0457
 -0.5357
  0.0296
 -0.7988
 -0.5027
  0.3328
 -0.7630
  1.5965
  0.9160
 -0.0870
 -1.0944
  1.7733
  0.7829
  0.7852
  0.6609
  1.3158
  0.1778
 -0.0915
 -0.3964
  1.5792
 -1.2953
  0.6502
 -0.2258
  1.3728
  0.5990
 -0.5334
  0.1386
  0.4469
 -0.3256
  1.0199
  0.7073
  2.0667
  1.7699
  1.2874
  0.2657
 -0.4100
  0.1803
  0.5804
  1.4055
  0.0111
  0.3264
  0.1424
  3.0684
  0.1530
  0.4575
 -0.4570
 -0.8071
 -1.1085
  0.7804
 -0.3769
  1.6507
 -0.6642
 -1.1415
 -0.9749
 -1.3107
  1.9003
  0.9407
  1.5142
  1.3973
  0.4339
  0.9401
  1.0267
  1.4933
  0.8625
  1.1322
  1.9362
  1.9916
  0.6271
  2.0667
  2.1391
  0.5136
 -0.1143
 -0.4907
 -0.0717
  0.3463
  0.2727
  0.6394
  0.5504
  0.2556
  0.6207
  0.4065
  0.6836
  2.4603
  1.9054
  0.9640
 -0.4100
  0.8968
  0.3383
  0.1916
 -2.0774
 -0.5320
  1.2723
  0.0578
 -0.7592
  2.2646
  1.2979
  0.9818
  1.7583
  0.4592
 -0.1837
  0.1604
  1.8136
  0.0782
 -0.2059
 -2.1169
 -0.0813
 -0.4820
 -0.2970
  0.8097
  0.8234
  1.4291
  0.2690
  2.1060
  0.5427
 -0.1891
  0.8239
 -1.5941
  3.2874
  0.4677
 -0.1070
  0.0177
  2.0204
  0.5323
  0.1754
 -0.1707
  3.3876
  1.4342
 -0.5878
  1.3895
  1.9137
 -1.2977
 -0.7009
  0.1937
 -0.0035
  0.3224
 -2.8941
 -0.3020
  1.3193
  0.7741
  0.5387
  1.1577
  2.8716
  2.6101
  1.2104
  2.8556
 -1.0172
  0.0243
  1.6348
  1.0482
  0.4285
 -0.1896
 -1.4156
  2.0362
 -0.2446
 -0.2349
  0.3232
  1.4948
  1.5257
 -0.3014
  0.3072
  0.2349
  1.3178
  0.0631
 -0.3064
  1.3051
  0.0416
  0.3128
 -0.0046
 -0.2702
  1.7477
  1.9460
  2.7832
  1.8657
  2.3522
  1.6165
 -0.1730
 -0.5655
  1.0532
 -0.8025
  2.1758
  1.4156
  0.1163
  1.1524
  0.9137
  0.8551
 -0.0100
  1.4772
  0.4579
  0.4876
  0.5338
  0.3912
 -0.1739
  0.3485
 -0.4219
 -0.2344
 -0.5852
 -2.5845
 -0.7903
  1.1720
  0.1977
  0.3628
  1.2447
 -1.5008
  0.2655
 -0.9282
 -0.3222
  0.5915
  1.4185
  0.9108
  0.1232
  1.3689
 -1.5137
  1.9290
  0.4253
  0.5072
  0.8117
  1.6375
  1.0934
 -0.2443
  1.7873
  1.8030
  3.0388
  1.5748
  1.5768
  1.5323
  1.8591
 -0.4260
  1.5920
 -0.3699
 -0.6282
 -0.7825
 -1.3464
 -0.2376
  1.8355
  0.1069
 -0.7044
 -0.0542
 -0.7121
  1.5210
  1.7094
  0.9637
  0.2786
  1.1053
  0.1538
 -0.9021
 -0.7895
  0.8948
  0.4275
  1.5102
  1.2940
 -0.4354
  0.2464
 -0.0647
  0.0089
  0.4430
 -0.6453
 -1.0938
  1.4767
 -1.0690
  2.4844
  1.1645
  0.6536
  0.0727
  1.4173
  1.3244
  1.6645
  2.4112
  1.7466
  1.4463
  2.6022
  1.2655
  1.2139
 -0.6171
  0.6223
  3.6917
  2.8420
 -0.1472
  1.0051
  0.4670
  2.9791
  1.5567
  1.6677
  0.5057
 -0.5116
 -2.3569
 -0.1297
  1.8901
 -0.2794
 -1.1341
 -0.5486
 -0.3043
 -0.5839
 -0.3817
  0.3127
 -0.7209
 -0.2999
 -0.7527
 -1.2488
  0.3438
 -0.9340
  1.0086
  0.7321
 -0.1372
  0.9753
 -1.7896
  0.3201
  1.4807
  1.8327
  1.0522
  1.0265
  0.5847
 -0.2252
 -0.3896
  0.8775
 -0.1468
  1.3870
 -0.2504
  0.4431
  1.4483
  0.4093
  1.4904
  0.4576
  1.2077
  1.6620
  0.8077
  1.8053
  0.7161
  2.4929
  0.5509
  0.8348
  2.2478
  1.8418
 -1.9399
 -0.2695
  0.7887
 -0.7948
 -0.2231
 -1.0391
 -1.9755
  1.8624
  0.8098
  1.2964
  1.3673
  0.8999
  0.5951
 -1.1058
  0.0333
  0.9262
  0.5710
 -0.4284
 -0.7852
 -1.3974
 -1.5191
  0.7884
  0.5349
  2.6747
  1.6850
  2.4541
  0.5391
  0.5678
 -0.3843
  1.6227
  0.0262
  0.1685
  1.0297
  1.2574
  0.5560
  1.1544
  1.0829
 -0.2152
  0.0156
  0.5524
  1.5422
  1.4950
 -0.1620
  0.2186
  1.3553
  1.4332
 -0.2974
  0.4511
  1.1137
 -1.8690
 -0.2772
  0.7924
 -0.4872
  0.0747
 -0.9092
 -0.9512
  0.0608
  1.6924
 -0.4051
  0.3471
 -0.2804
  0.4888
 -0.8382
  0.2907
  1.2709
 -1.2714
 -0.8459
  0.0931
  0.5774
 -0.6496
 -0.0372
 -0.0105
  0.1881
 -0.7229
  0.0792
 -0.1622
  1.0479
 -0.0960
 -0.1600
 -0.4500
 -1.0594
  0.6595
  0.6013
  1.0189
  1.1394
  0.5098
  1.3369
 -1.0730
 -0.2536
  3.0535
  0.7173
  2.0727
  1.4629
  0.3629
  1.9538
  1.5533
  1.7432
  0.3150
  3.0895
  0.1441
 -0.5046
  0.2766
 -0.0098
  0.7529
 -0.1762
  2.4103
  0.1211
  1.2132
  0.5893
  1.6663
  2.1484
 -0.3480
  0.9909
  0.9513
  0.7870
 -1.7721
  0.7912
 -0.2925
  1.5361
 -0.0641
  1.4937
 -0.2837
 -0.5852
  1.7148
  1.4834
 -0.5957
 -0.5563
  1.4202
  1.1139
  1.0016
  0.3591
  0.8338
 -0.0200
  0.1794
  1.0058
 -0.4821
  0.7581
  0.6138
  1.8949
  0.3902
  0.5157
  2.7695
  1.2023
  0.4279
  2.5293
  1.8040
  0.7622
  1.3580
  2.0576
 -1.4054
  0.9415
  1.5915
 -0.1176
 -0.4397
 -0.4349
  1.2321
  0.4294
 -0.9560
  0.3796
 -0.7256
  0.7133
  1.1080
  0.1468
  0.0411
  2.2929
  0.6316
  0.5853
 -0.9802
  0.9908
  0.2422
 -2.0519
 -0.8010
  0.1548
 -0.1160
 -0.1558
  0.6486
  0.7693
  0.4834
  0.0208
 -0.0795
 -1.0668
 -1.3892
  0.0797
  0.8897
  1.3288
  0.5291
 -0.5492
  1.2521
  1.6147
  0.1973
  0.8293
 -1.6394
  0.2142
  1.1838
  1.3407
  0.6382
 -0.0335
  1.5385
  0.6899
  0.9429
  0.3169
 -0.4660
  0.8745
 -0.1279
  1.3947
 -0.5098
 -0.1902
  1.1792
 -0.4725
  1.9196
  1.3379
  0.0673
  0.2855
  1.0812
 -0.1560
 -1.0904
 -0.0693
  0.0418
  1.1500
  0.1645
  0.4687
  0.8862
  1.7932
 -1.6829
  2.3866
  0.6344
  0.2574
  2.6439
  0.7563
  0.1801
  0.5222
  0.4995
  0.2530
  2.3603
 -0.1601
  1.8622
 -0.3433
  1.4434
  0.6389
  1.1026
  0.8460
  0.2153
  0.2270
  0.8038
  1.5537
  1.0169
  0.4464
 -1.5309
  1.3133
  1.8218
 -1.9623
  0.1950
 -0.5387
 -0.0226
  0.2494
  1.5253
 -1.3937
  1.2874
  2.7512
 -1.2288
 -0.0715
  0.6471
  0.1562
 -0.0432
  0.7117
 -0.3367
 -1.0532
  0.8384
 -0.1312
  0.6571
  0.5534
  0.2304
 -1.5987
 -1.9738
 -1.1655
  0.0922
  0.2188
 -0.1341
  1.1572
 -0.6883
 -0.2567
  3.1340
  1.5854
 -1.8538
  2.9004
  3.1624
  2.6021
  1.1631
  0.8688
  1.1366
  1.4786
  1.4542
  0.2267
  0.8324
  0.3476
  1.8441
  1.5986
  1.4828
  2.1348
  0.4264
 -0.2599
  0.0124
 -0.4973
 -0.0989
 -1.0021
  0.3867
 -1.2221
  0.9612
  0.6209
  0.9313
 -1.1496
  0.0493
 -0.9382
  0.6532
  0.4046
  1.7479
 -0.7891
 -0.0734
  0.4484
  0.7912
  0.8431
 -1.4191
  0.3306
 -0.8075
  1.1079
 -0.4095
  0.7746
  0.1153
  0.6681
  1.2301
 -0.1471
  2.3321
  2.0279
  1.7329
  1.5797
  1.3264
  1.0826
  1.1672
  2.1596
  0.2233
  1.0472
  1.4141
  0.3688
  3.2731
 -1.2973
  1.7902
  0.7519
 -0.6775
 -0.1660
  0.5604
  0.2420
 -0.2340
 -0.0021
 -0.2163
  2.1085
 -1.0034
 -0.8599
  0.4535
 -0.5121
  0.3242
 -0.2484
  1.7565
 -0.1810
  0.7489
  2.1222
  0.0898
  0.9413
 -0.4087
  0.1029
  0.0839
  1.0873
 -0.4155
 -0.1335
 -1.0481
  1.6287
  1.0319
  1.3821
  0.2540
 -1.0321
 -0.0868
  0.9426
  0.0244
  2.0310
 -1.4661
  1.4249
  2.4610
  1.2073
  0.8903
  2.5020
 -0.0626
  1.2325
  0.3687
  1.2880
 -0.0536
  1.2178
 -0.1083
  0.7185
  1.9836
  1.9895
  0.5832
 -0.2703
 -0.1394
 -1.3790
 -1.2633
 -0.0365
  0.0927
  0.5114
 -1.0059
 -0.1422
 -1.9468
 -0.7289
 -0.5217
  1.7765
  0.7512
 -1.8069
  0.3108
 -0.4054
 -0.7866
 -0.1028
 -0.1570
 -0.5415
  1.0934
 -0.4104
 -0.5743
  1.0649
  1.2204
  1.1739
 -0.8054
  1.1340
  1.7302
  2.4039
  0.8643
  1.7043
  1.0529
  2.3288
  3.1260
  2.0885
 -0.1631
  2.7071
  0.4836
  0.1616
  2.0697
  1.0634
  1.5721
  1.5385
  1.4192
  3.0171
 -0.8831
  0.3061
  0.4398
 -1.8922
  0.5826
 -0.8189
  1.4337
 -0.4133
 -1.2969
 -0.4322
 -1.2330
 -0.6419
  0.5336
 -0.2887
  2.4023
 -0.5586
  0.7574
 -0.7170
  0.4473
  1.0097
 -0.2260
  0.4884
  0.1430
 -0.4771
 -1.8925
 -0.4507
  0.2213
 -0.5899
  0.2818
  0.3618
  0.8276
  0.8201
 -0.8813
 -0.5491
 -1.2701
  1.7920
  2.1862
  0.7757
  0.9921
  0.9338
  0.0809
  1.6657
  1.8268
  0.6505
 -0.8575
 -0.0190
  1.3977
  1.4767
  2.6754
 -0.2300
 -0.1236
  2.0332
 -0.6538
 -0.8092
  0.3628
 -0.6275
 -0.0017
 -0.3331
 -0.0233
 -1.2726
 -1.3070
 -2.8459
  2.0920
 -1.1198
 -0.4835
 -0.5981
 -1.0269
 -1.6891
  0.5556
  0.6732
  0.2691
  0.5301
  0.2409
 -0.4941
 -0.0441
  1.5063
 -1.1230
  1.3220
  2.4981
  0.9180
  1.3256
  0.7154
  1.1433
  1.5460
  0.9967
  0.1360
  0.7045
 -0.0095
  2.5648
  0.4124
  0.7204
  2.6875
  0.2029
  0.6200
  1.0284
  1.3347
  2.7369
  0.6882
  0.2671
 -0.4818
 -0.4306
 -1.5418
 -1.0259
 -0.6642
  0.8190
 -0.4672
 -1.3631
  1.1023
  1.0422
  0.6897
 -0.6886
  0.7235
 -1.5059
 -1.0558
 -2.0258
  0.8143
  1.7396
 -1.1647
 -0.0069
  0.9067
  2.0647
 -2.2128
 -0.2625
 -0.2348
  0.0839
 -2.0861
  1.6168
  1.8230
 -0.0403
 -1.6485
 -0.1074
  0.4361
  0.1800
 -0.2912
 -0.3069
  1.0963
  0.9860
 -0.6466
  0.5590
  0.0676
  0.4222
 -0.6195
  2.6005
  1.2032
 -0.1152
  0.7673
  0.4395
  0.0525
  3.3842
 -0.1220
  2.4949
  0.2756
  1.6650
  0.0304
 -0.5331
 -0.7718
  0.7244
 -0.4942
 -0.8868
  1.4085
  1.6292
 -0.0684
 -0.6152
 -0.1432
 -1.0092
 -0.3900
 -1.6158
 -0.5118
 -1.7379
  0.8450
  0.3790
 -0.5868
  1.2743
  0.5662
  1.5342
  1.8803
  0.7103
 -0.1744
  1.3298
 -0.0455
  2.0123
  2.1661
  0.7630
  0.1299
  1.9721
 -0.8603
  0.8245
  0.7655
  2.4525
  1.7631
  1.1402
  2.0192
  2.8706
  0.5191
  0.3079
  0.7055
  2.4001
  0.4194
  1.9448
 -1.7994
  0.5777
  0.7101
 -1.9289
 -0.8867
  0.0824
  0.2545
  0.4311
  0.4522
  1.9328
 -1.1528
 -0.7003
  1.4608
 -0.9918
  0.9723
  0.7158
  0.7401
  2.3548
  0.2212
 -1.8905
 -0.5717
 -1.5334
  0.7415
  1.4221
 -1.2245
  1.2510
  0.8446
 -1.4908
  1.7167
  0.0549
  0.8011
  0.9911
 -0.9935
  0.8396
 -1.0699
 -0.5040
  0.5107
  0.1145
  0.4912
  0.8992
  1.3960
  0.6420
  1.3104
  2.3389
  1.1635
 -0.4555
  1.3552
  0.3371
 -1.5688
  0.4805
  0.1190
  2.0136
  0.8790
  1.5689
  1.8195
  0.4575
 -0.0119
  1.4426
 -1.5308
  0.2984
  0.2502
  0.1518
  1.5217
 -0.3208
 -0.1088
 -0.3128
 -0.4383
 -0.3105
 -0.0533
  1.0537
 -0.8797
  1.7778
 -0.8837
 -0.2216
  2.2416
  2.1598
  2.5002
  1.9269
 -0.6057
  1.6280
  1.4513
  2.0292
  0.5035
  0.2712
  1.6066
  0.8998
  1.2470
  0.2247
  0.4452
  0.5700
  1.1001
  1.1008
  0.9057
  0.1839
  1.4572
  0.7518
  0.6090
  1.2043
 -1.1519
  0.3261
  2.2318
 -0.1930
  1.0985
  0.8454
 -2.0136
 -0.9738
  0.6497
  0.0973
  1.6262
 -1.9708
  0.5527
 -1.0739
  0.5348
 -1.0997
 -0.5800
 -0.2207
  0.2476
  0.3391
 -1.0181
 -0.1945
 -0.5975
  1.6146
 -0.0360
 -0.3068
  0.2107
 -0.7894
 -1.8640
 -0.4169
  0.7008
  0.3655
  0.5651
  0.5329
 -0.7718
  0.4861
 -0.5814
 -0.0068
 -0.0464
  0.0665
  0.3725
  2.8042
  1.9437
  0.4936
  0.5704
 -0.7525
  1.9315
  1.6334
  1.5611
  1.8531
  1.5383
  2.1041
 -0.0387
  0.5999
  2.1075
  1.6470
  0.2429
  1.5560
  1.1470
  2.6831
  0.2161
  3.7020
 -0.4234
 -0.9785
 -0.0641
  0.0003
  0.5258
 -0.0129
  0.8594
 -0.7987
  1.1971
  1.2377
  0.2303
 -0.2817
  0.9144
  0.8039
  1.8654
  0.8317
  3.5761
  1.2835
  2.7451
 -0.3493
  3.6506
  0.9725
  0.8066
  0.2886
  1.3716
  1.0500
  1.7711
  1.2807
 -0.1342
  0.6610
 -0.7561
  0.4964
  0.5114
  0.7618
 -1.6087
  2.0261
  1.9428
  1.7669
  2.2918
  0.0332
  0.5077
 -0.9214
  0.4744
  0.7612
  1.2358
  0.2184
 -0.5395
 -1.6624
  0.6785
 -2.6201
 -0.4495
  0.8749
  0.6771
  1.2379
 -1.9280
  1.8429
 -0.6623
  0.6633
  0.0195
 -0.3562
 -1.0078
  0.3965
 -0.0235
 -0.6519
 -1.0644
  1.8561
  1.8827
  1.4523
  1.7001
 -0.5889
 -0.7120
  2.6150
 -0.9384
  0.4249
 -0.8328
  0.6025
 -0.9109
  0.2761
  1.0010
  3.0467
  0.9389
  1.1662
 -0.0887
  2.7199
 -0.3188
  2.7667
 -0.2248
  0.6690
  0.4206
  0.1420
  1.6973
  2.5447
  0.8298
  2.4594
  1.0135
  1.9832
  2.0852
  3.8609
  0.5971
  1.7439
  0.5565
 -0.4710
 -0.1555
  0.9522
  0.4270
  1.7568
  0.3076
  0.2593
 -2.3374
 -0.2197
 -2.3488
  1.1738
 -1.3928
  1.9129
  0.8371
  0.5482
  0.4787
  0.7969
  1.1788
  1.4314
  1.0995
  0.5772
  1.0942
  1.5971
 -0.1045
  0.7713
 -0.2098
 -0.1878
 -0.0188
 -0.5522
  1.6412
  0.9590
  1.5387
  0.5477
 -0.1625
  1.2712
  1.5655
 -1.1169
 -0.2504
  0.2462
  1.7831
 -1.7358
 -0.7036
  1.0210
  1.1533
  0.0875
 -1.5465
 -0.9499
 -0.7919
  0.3797
  0.9449
  1.8356
 -1.2647
 -1.1874
  0.4520
 -0.0938
 -2.2516
 -0.3261
 -0.2233
 -0.6450
  2.0489
  0.2858
 -0.4093
 -0.5503
 -1.4385
 -1.8579
  0.5549
 -0.8884
 -0.1917
  1.0477
  0.2593
 -1.0415
  0.1815
  0.2117
 -0.0686
  0.8808
  0.0470
  2.5044
 -0.6097
  1.2777
  2.2870
  0.3199
  1.7351
  2.0242
  2.1028
  1.3710
  1.4662
 -0.1187
  1.9727
  1.6029
  1.1512
  0.3591
  0.9109
  1.9073
  2.8008
  2.6560
  1.5460
  1.0783
  0.9066
 -1.1586
 -0.0776
  0.0354
 -0.1256
  0.0622
  1.9970
 -0.9491
  1.1662
  0.6183
  0.7785
  2.4502
  1.6334
 -0.8394
  1.4778
  0.2673
  0.3653
  1.0701
  0.2513
  1.1031
  1.5847
  2.1929
  1.4245
 -0.3859
  1.9877
  0.7626
 -1.7752
  1.9523
  0.8718
 -0.5655
  0.5633
  0.6976
  1.6442
  1.5483
  0.6827
 -0.0432
  0.4062
  3.5237
  1.7477
  2.3267
  0.9540
  1.6615
  2.7925
  0.4667
  1.8013
  0.6912
 -0.0821
 -1.2321
 -1.9666
  0.3609
  1.6417
 -1.3599
 -0.6310
  0.2801
 -0.0238
  0.8708
 -0.6098
  0.1117
  0.8708
  0.6512
 -0.0662
  0.1672
  0.2900
  1.4404
  0.7190
 -1.3648
  0.4860
 -0.0565
 -1.2363
  1.4457
  3.0362
  0.4919
  1.3521
  1.2877
  0.4750
  2.7448
  0.5484
  2.3284
  0.8682
  1.9640
 -0.4344
 -0.7072
  0.4137
  0.6292
 -0.5906
  0.0935
  0.6537
  1.6287
  0.4727
 -0.4705
  2.4667
  0.0063
  1.0666
  0.5820
  1.4579
  0.7289
  1.2961
  2.5506
  0.2875
  0.8416
 -0.2085
 -0.3592
  0.8128
 -1.7722
  1.1258
 -0.2726
  0.5459
 -1.0665
 -0.7349
 -1.9992
  2.4051
  0.4652
  0.1482
  1.9183
 -0.5237
  2.0534
  1.4931
  0.9052
  1.8488
  1.9809
  0.1732
 -0.0739
 -0.3123
  0.5244
  0.7965
  1.0999
  0.7354
  2.3114
  0.9648
 -0.4500
  0.7997
  0.2606
  2.7622
  0.8078
  0.0737
  2.3617
  1.2009
  2.2449
  1.1454
  0.2243
  0.8199
 -1.0971
  0.1955
  1.4960
 -0.1562
 -0.5091
 -0.8335
  0.2311
 -0.3100
 -0.3443
  2.2643
 -0.5603
 -1.5265
  0.6500
 -0.8990
 -1.0666
  0.1181
  0.5713
 -0.2804
 -0.6967
  0.4550
 -0.7567
  1.4244
  0.6333
 -0.0599
 -1.2478
  1.7776
 -1.4554
  1.1644
  0.1458
  1.9338
  3.2939
  0.9596
  0.0067
  0.4529
  1.6103
  2.1864
  1.1740
  1.0471
  0.9263
  2.8761
  1.9958
  0.5966
  0.5881
 -0.3314
 -0.2212
 -0.0230
  1.8355
  0.7936
  1.3632
  2.7717
  0.5721
  0.3583
 -0.1936
  1.6767
  1.5180
  0.7816
  0.7077
  0.1920
  1.9352
  1.1824
  0.7044
  1.0567
  0.6816
 -0.5668
 -2.8734
  1.8152
  1.5217
  0.0601
 -0.0982
  1.0049
 -0.5325
 -1.3526
  0.2557
  2.6593
  2.3461
  1.6745
  1.3303
  0.3707
  0.0605
  1.8826
  1.7895
  0.9996
  1.4209
  0.8260
  2.0267
  0.1407
 -0.8803
  0.3525
  0.2514
  1.5518
  3.0510
  1.3362
  2.2214
  0.4109
  1.5241
  2.1755
 -0.4756
  0.1818
 -0.4655
  2.4242
 -0.8267
 -0.9794
  2.1438
  1.3629
  0.6857
  0.9541
 -0.4125
  0.0619
  0.6034
 -1.1607
  2.0992
  0.5471
  0.4778
  0.5390
  0.1683
  1.4864
  0.5173
 -0.0277
  0.3716
  1.9539
 -1.7593
  0.4325
  0.0767
 -1.6166
 -1.1375
  1.1167
  0.5647
  1.9726
  0.4031
  0.4293
  0.9291
  1.7818
  1.2414
  2.6056
  1.9066
  0.2968
  0.5038
  1.7371
  1.2020
  1.0076
  0.2036
  0.2057
  0.9303
  1.0438
  2.2711
  0.4873
  0.3147
  1.2508
 -0.7386
  0.8381
  1.3598
 -0.0621
  1.8053
  0.6831
  1.4555
 -0.5853
  0.3532
 -1.7368
  2.2790
  0.6614
 -0.8643
  0.3476
  0.0192
  1.3550
 -0.0221
 -0.1939
 -1.3202
  0.2362
  2.5507
 -2.2400
  0.0163
 -0.1940
  2.6369
  1.1601
 -0.0917
  0.3585
  2.6288
  0.6612
  3.5829
  0.6776
  1.2317
 -0.8212
  0.9182
  0.7538
  0.8468
  2.1407
  2.0916
  0.3506
  0.6370
  1.8844
  0.2989
  1.9831
  1.7690
 -0.0142
  1.7071
 -0.0979
  0.5921
 -1.0877
 -0.2591
 -0.1192
  0.6119
 -0.1604
  0.8847
  2.3112
  2.3672
  2.0429
  0.2373
  1.7020
  0.5146
  0.1226
  0.7145
 -0.3320
 -0.7394
 -0.1155
  0.8220
 -2.8656
 -0.3125
 -1.0315
  0.5440
 -0.3549
 -1.1652
  0.8009
  0.1497
  0.6537
  0.9155
  3.2016
 -0.4384
  0.3339
 -0.1100
  1.1433
  0.7088
 -0.6598
  0.7762
  1.1458
  0.0642
  2.2126
  1.1639
  0.2962
  2.0178
  1.4824
  0.0249
  2.6248
  2.5233
  1.5852
  2.3267
  0.8191
  0.9973
  0.5750
 -0.2456
 -0.5668
  0.5065
  1.4876
  0.8591
  1.1850
  1.8913
  0.9896
  1.4015
 -0.0126
  0.2585
 -0.7415
 -0.7792
 -0.9080
 -1.2071
  0.0230
 -2.4561
 -0.1847
  1.0029
  2.2054
  1.7842
  2.3290
 -0.2905
  0.4917
 -0.3519
  2.7416
 -1.0177
  1.2932
  1.9323
  1.3755
  0.7414
  1.9546
 -0.1017
  0.7094
  1.1548
  0.4514
  0.7531
 -0.0746
  0.1919
 -0.5055
  1.4016
  1.2783
  1.3012
  0.8263
  2.3029
  0.4269
 -1.3226
  2.9883
 -0.3776
  1.3827
  0.6385
  1.5129
  1.0569
  1.4394
  1.2700
  1.4664
  2.1039
  0.8396
  0.7724
  0.4940
  0.4168
  3.4394
  1.6930
 -1.8066
 -0.7448
 -1.1568
  1.1043
 -0.3138
 -1.4610
 -1.3952
 -1.3937
 -0.0177
  1.0354
  0.8422
  2.1298
 -0.0157
 -0.0520
  0.9976
  2.1178
  1.2298
 -0.1663
  1.5986
  2.0715
  0.6655
 -0.0129
  0.4469
 -1.8062
  1.7767
  2.0704
  3.1114
  0.1654
  3.1167
  1.1675
  1.6402
  2.6596
  1.5045
  1.4171
  1.9969
  0.9152
 -0.4595
 -1.7809
 -1.0556
  0.4706
  1.2743
  1.7461
  1.8876
  0.7929
  0.1276
  0.8848
  3.0655
  2.3857
  2.4022
  1.0778
  0.4821
 -0.5569
  0.8440
 -0.2530
  1.0220
  0.8689
  1.7864
 -0.4671
 -0.4158
  0.9217
  0.8659
  0.1101
  1.5543
  1.7426
  0.6851
  1.3159
  1.2964
  0.8115
  1.2477
  1.9245
  0.5888
  0.0473
  0.4629
  1.7692
  0.0128
  0.3681
 -0.0145
 -0.3195
 -0.7533
  0.3479
  2.2645
 -0.9712
  1.4569
  2.0601
  1.0863
 -0.9583
  1.3606
  0.9370
 -1.0785
  1.1226
  2.4091
  1.4620
 -0.3503
  0.8848
  2.4002
  1.6424
  0.2620
  1.8617
  2.5300
  3.0005
 -1.5915
  2.0563
 -0.7669
  0.7008
  1.4222
  1.7879
  0.8176
 -0.2729
  0.4969
  1.3067
  0.8078
  0.4442
  1.3717
  1.6234
  0.7852
  1.3541
  0.8162
  1.6322
 -0.0875
  1.8238
  0.1807
  0.8250
  2.3860
  0.4648
  2.0181
  2.2088
  0.3221
  0.6401
  0.1299
 -0.1998
 -1.1023
  1.8496
 -0.1814
  2.1868
 -1.0986
 -0.7129
 -0.2885
  0.1864
  0.4363
  0.0578
  0.8682
  1.0330
  0.3281
  1.6545
  0.8548
 -0.3302
  2.5085
  0.9887
 -0.7300
  0.5193
  0.5745
 -0.0715
  0.9361
  1.6869
 -0.6209
  0.0569
  1.4783
  1.8422
  2.8073
  1.6954
  1.5393
  1.8066
  2.2483
  2.4179
  0.1006
  1.8085
  0.3697
  2.5552
  0.2439
 -1.7455
  0.5191
  2.1053
  1.7140
  0.5037
  0.9079
 -0.2739
 -0.4462
  0.5953
 -0.0875
 -0.2572
  1.0707
  2.2568
  1.6424
  0.2564
  2.2654
  2.5224
 -0.9397
  1.0906
  1.2586
  0.8715
  2.6548
  1.7738
 -0.2823
  1.4004
  0.3161
  3.7019
  2.4087
  3.5276
  1.1414
 -2.0927
 -0.5732
  1.4016
  2.4132
  0.0382
 -0.6262
  2.8688
  0.8210
  2.5827
  2.1605
  0.3908
  1.9893
 -2.0397
  2.0427
  0.6536
  1.5885
  0.6643
  1.4345
  2.0935
 -0.3838
  0.3087
 -0.5504
  0.5335
 -0.1382
  0.4680
  0.3849
 -2.5514
 -1.5377
 -0.3971
 -0.1601
  1.3828
  1.7277
  1.4583
  1.5161
  1.1622
  2.0661
  2.1292
  0.0044
 -1.4473
 -0.7931
 -1.6190
  0.4809
 -1.6652
  1.0081
 -1.3866
 -0.5566
  0.7950
 -0.0094
 -1.4374
  0.4977
 -1.6921
 -0.0309
 -0.4680
  1.7980
  1.8789
  1.2792
 -0.1198
 -0.0830
  1.6186
  0.7868
 -0.0088
  0.4496
  1.5473
  2.8741
  0.0586
  0.3377
  1.2518
  1.6952
  1.2943
  2.2865
  0.5218
  0.0577
  1.1255
  0.1215
 -0.7247
  1.7781
  0.7834
 -1.1465
 -0.9035
 -0.6152
  3.5951
  2.0457
  0.1277
  0.1978
 -0.2624
  0.4073
  2.1105
  0.8708
  1.0649
  1.8453
  2.3927
  2.4552
  1.7728
  2.3556
  1.1943
  0.9894
 -0.2268
  1.2185
 -0.6336
  2.2817
  1.2364
  1.3369
 -1.0040
  1.3971
  0.5499
  1.4408
  2.0269
 -0.3123
  1.7283
  1.9062
 -0.1382
  1.3675
 -1.2599
  0.8864
  1.0824
 -0.1651
  0.2770
  0.8174
  1.5959
  2.2754
 -0.7358
  1.9965
  3.2023
  0.5272
  0.4233
 -0.5924
 -1.9205
 -1.4772
 -1.2070
  0.0516
  0.8794
  2.8627
 -1.2011
  0.1484
 -0.2973
  0.6361
  0.4279
 -0.9844
 -0.8717
  0.5139
 -1.1109
 -0.0628
  1.1003
  1.3237
  0.2620
  0.7842
  3.2762
  0.6984
  0.1812
  0.1593
  0.4195
  1.8723
  0.6160
  1.1012
  0.9497
  3.1727
  1.0809
  0.3595
  0.1616
  0.5335
 -0.7293
  0.8295
 -0.3565
  2.5550
  1.4017
  2.6955
  1.3820
  1.8544
  1.6116
  1.0413
  2.4931
  2.3058
 -1.0172
 -0.1015
 -0.8356
  1.2927
 -1.1277
  0.6627
  0.0953
  1.2012
 -0.9634
  1.6146
  0.9587
 -0.5620
 -0.1586
  2.0034
  0.7296
  1.4347
  1.3171
  1.3805
  1.9071
  1.5903
  0.0513
  0.7894
  0.6318
  1.6109
  1.1493
  1.1562
 -0.4095
  2.6911
  0.4682
  1.4386
  1.0540
  0.7904
  1.2928
 -0.4302
  0.3031
 -0.9125
  1.6658
  1.1929
  1.9077
  1.4130
  1.0021
  1.9525
  2.7609
  0.3304
 -0.9063
  1.8720
  0.8363
  1.6236
  0.7375
  2.0815
  0.9570
 -0.3423
  1.5206
  2.2356
 -0.9859
 -0.5213
  0.2971
  0.8297
 -1.2778
 -0.6695
  1.0294
  0.2878
  1.3149
 -0.6401
 -0.7984
 -0.1869
  0.6261
  0.5621
 -1.1915
  0.4557
 -0.4337
  0.7420
  0.7364
 -0.8170
  0.6208
 -0.2054
  1.4715
 -0.3487
  0.8129
  0.0725
  2.2333
  0.4630
  1.9997
  0.7689
 -0.2283
  1.6911
  0.4627
  1.3006
  0.6752
  1.2714
  0.4522
  1.5202
  0.3381
  0.7460
  0.6634
  1.2539
  0.5474
  1.3246
  0.5679
  0.6545
  0.9754
  1.8302
  2.8653
  2.2761
  0.2459
  0.2422
 -0.0547
 -0.1661
  0.1796
  0.9097
  1.8410
  1.7654
  1.7344
  0.0340
  0.8955
  2.2704
  1.0818
  1.8021
  0.8768
  0.4412
  2.1166
  1.2156
  1.6085
 -0.1821
  0.0681
  1.5657
  2.2373
  1.4264
 -0.7371
  0.0192
  1.7881
 -0.0272
  2.3230
  2.0493
  2.3929
  0.2850
  1.9886
  1.0253
  3.5733
  2.4045
  1.3202
  0.9075
  3.2705
 -0.0236
  1.1770
  0.7559
  1.9281
  2.4545
  0.4942
  2.0598
  0.3893
  0.5199
 -0.4951
 -1.1558
  0.2130
 -0.7820
  0.2241
 -0.1810
 -0.8585
 -1.6114
  0.6884
 -0.4772
 -2.5293
  0.0550
 -0.0010
 -0.2301
 -0.3899
 -1.9274
  1.3185
  0.9857
  0.1351
  0.4756
 -0.5660
 -0.0247
 -0.1681
  2.8950
  0.8313
 -0.1904
  1.3659
  2.0753
  2.4624
  1.8682
  1.0907
  0.4394
  0.6438
  1.1707
  1.4323
  1.1544
  0.9701
  0.3714
  1.0603
 -0.2656
  0.9214
  0.5593
  0.8895
  0.5018
  0.8293
  0.4917
  0.2691
  1.0417
  0.0207
  2.4349
  2.4158
  0.7528
  2.2290
 -0.2650
  0.6712
 -0.4558
 -0.0514
 -0.1629
 -0.9351
  1.3515
  1.4886
  2.2013
 -0.7325
  0.4300
  0.0924
  2.3328
  2.6163
  0.7464
  2.1073
  1.9554
  1.7701
  0.8924
  1.3298
  2.2875
  0.8607
  0.6487
  1.1655
 -0.8991
  1.4793
 -0.1733
  0.9266
  1.2021
  1.4748
  2.2099
 -0.3016
  2.1656
  0.9316
  0.9834
  1.3019
  2.7277
  3.5947
  1.1527
  1.4152
  0.2894
  2.1302
 -0.3782
  1.0189
 -0.1522
 -1.3880
  1.0820
 -2.5384
 -0.8672
  0.8633
 -0.6849
 -0.6449
 -0.7962
  0.8690
 -0.7484
  0.5138
  0.8376
  1.2029
 -2.8677
  1.3620
 -0.1583
  0.2792
  0.2404
 -0.9063
  0.9363
 -0.7369
 -1.5827
  0.5527
  0.8685
  1.0648
 -0.1257
  0.2490
  0.2817
  1.2099
  1.1844
  2.2103
  0.6807
  1.8861
  1.0110
  0.6726
  1.6519
  1.2776
  0.3156
 -0.9928
  1.0844
  0.9551
  2.0562
  1.4274
  1.3785
  2.6227
  1.1352
  1.0050
  0.9654
  0.3774
  1.1710
  1.7877
 -0.1212
 -0.5277
  0.2149
  1.7243
 -0.1763
  0.2273
 -0.3153
 -1.0668
  0.8789
  0.0354
 -0.7653
  3.5419
  1.4534
 -0.1571
  1.2217
  2.2238
 -0.1360
  1.7644
  1.3325
 -0.4009
  0.9956
  1.2473
  0.8494
  0.9185
  0.5906
  1.5161
  0.2358
 -0.6289
 -0.7384
  0.1658
  2.0541
  2.5633
 -0.5871
  1.2303
  2.4515
  2.2904
  0.9250
 -0.5528
  1.3179
  1.0425
  2.3094
  0.1324
  0.9304
  0.7051
  1.9657
 -0.2327
  0.9453
 -0.8730
  1.2341
  0.7887
  1.7608
  1.6215
 -0.2927
 -0.0868
  1.9111
  0.3967
  0.4353
 -0.2920
 -0.0651
  0.1005
 -0.6044
 -0.8784
 -0.7343
  0.3699
  0.4422
 -0.2137
 -1.2593
  0.4307
 -1.3048
  0.9015
 -1.6494
 -0.6965
 -0.0722
  1.7557
  0.3961
 -0.7478
  1.2945
  1.6812
  1.3610
  0.9619
  0.4514
  1.8709
  0.9140
  1.6183
  1.5050
 -0.5750
  1.3534
  0.7106
 -1.9732
  2.1827
  1.4108
  2.1082
  1.0609
  1.8037
  1.5609
  0.3914
  0.4252
  0.9675
  2.7937
  2.0040
  2.9681
  0.2671
  1.8210
  1.5407
 -0.1238
 -0.3322
 -0.5009
 -1.5863
 -0.1066
  0.5240
  0.0128
  1.0312
 -0.0881
 -0.5103
 -0.2339
  1.9325
  1.3871
  2.1985
  0.9863
  1.0421
  1.3066
  1.0570
  2.4449
  1.1961
  1.2329
  1.1702
  1.7812
  1.0046
  1.6029
  1.6685
  0.3201
  1.6332
  1.0229
  2.0810
  0.7492
  1.1153
  0.9258
  0.7754
  1.1067
  0.2447
  1.1502
  0.8803
 -0.2724
 -0.0114
  2.2572
 -1.0731
 -0.1201
  1.8158
  1.8899
 -1.0971
 -1.5086
  0.0558
  1.6533
 -0.1962
 -0.1348
  0.7536
  1.2970
  1.2368
  0.2444
 -1.7988
  1.3968
 -0.6518
  0.2085
  0.8217
 -0.6529
 -1.4244
 -0.5268
  0.4976
 -0.9066
 -0.0493
 -0.8449
  1.3101
 -1.1616
 -1.9163
 -1.5777
 -0.3451
 -0.4210
  0.6731
  1.6048
  1.6151
  0.3027
  0.9852
  0.1670
  1.9296
  2.8052
  0.6984
  0.3045
  2.9606
  1.1189
  1.2361
  1.1036
  0.8996
 -0.2447
  0.5598
  0.0508
  0.7789
  1.0315
  1.9231
  1.2261
  1.7009
  0.5565
  2.4648
  2.2413
 -0.2948
 -0.1337
 -0.1466
 -0.7480
 -0.8790
 -0.1624
 -0.4288
  1.2677
  0.2821
  1.4182
  1.5142
  1.1682
  0.9539
  3.7675
  1.3929
  1.5750
 -2.1483
 -0.4490
  0.0554
  0.6307
 -0.5958
  2.0325
  1.5641
 -0.1036
  1.7099
  1.8199
  0.5303
  1.8601
  0.4553
  2.5275
  0.3289
  1.7261
  1.4888
 -0.8172
 -0.1705
  1.6354
 -0.0622
  1.9866
  1.1522
  1.2288
  1.3147
  1.3378
 -0.0084
  0.0624
  1.4436
  0.9165
 -0.3476
 -1.2034
 -0.5851
  0.2402
 -1.4041
 -0.2351
  1.2286
 -0.3182
 -0.8288
  0.9644
 -1.3051
  0.5747
 -1.1071
 -0.4599
 -0.1104
 -0.0320
 -0.1430
  0.7527
 -0.0126
  0.1632
  0.4571
 -0.6446
  1.1446
  0.7786
  1.1864
  0.0749
  0.2789
  1.1011
  0.8180
  0.6706
  0.9353
  0.8045
  1.0628
  0.9694
  0.8275
  2.1853
  1.2568
  0.0967
  1.3937
  1.3563
  1.9184
 -0.3573
 -0.2029
  0.8435
  1.7639
  1.7224
  0.2113
  0.4601
  2.9815
  1.0939
  2.1098
  0.4816
  0.8236
  2.0968
  0.0832
  0.7878
  0.0052
 -0.2382
 -0.6540
  0.6729
 -0.6162
  1.2724
 -1.8365
 -0.8652
 -0.1780
  1.3890
 -0.4402
  1.2693
  2.7102
  2.2037
  1.6514
 -0.0174
 -0.3569
  0.6363
  2.3369
  0.3380
 -0.7168
  0.5091
  1.1793
  0.2325
  0.7437
  0.5871
  3.0322
  1.3850
  1.4549
 -0.0516
  1.6257
  1.4175
 -0.0178
 -0.0532
  0.2500
  1.2686
  0.8772
  1.7521
  2.2113
  0.0233
  1.0534
  0.9221
  1.7986
  0.6301
  2.3112
  2.2666
  0.9844
 -1.8344
  0.2634
  0.7792
  0.3276
 -0.7662
  1.8076
 -0.0592
 -0.6122
  1.6886
  0.2688
  0.9392
  1.5567
 -0.2133
  0.9774
  1.2162
 -0.0234
 -1.0974
  0.5614
  0.6830
  0.1502
 -0.0170
  2.0480
 -1.5526
  1.0818
 -0.6424
  0.6186
  1.4146
  0.3376
  2.2755
  0.8659
  1.4737
  1.1523
  0.1744
  0.7196
  1.3679
  2.3143
  0.0322
 -1.2119
  0.5105
  0.2998
  1.5807
  2.3680
  0.7474
  1.3036
  1.2541
  0.7639
  2.5900
  1.4443
  0.0839
  1.1287
  0.4348
  0.0187
  0.2361
  0.3870
 -2.2344
 -0.5336
 -2.1278
  1.0904
  0.0883
  0.1665
  0.1070
  1.9043
 -1.0906
  2.6477
  0.5150
  1.2010
  0.6288
  1.3770
 -0.0298
 -0.2496
 -2.0294
  2.2656
  0.3766
  0.0785
 -0.3721
  1.1381
 -0.3518
  1.9479
  1.0856
  1.4941
  1.5675
  2.6301
  0.6128
  1.3585
  0.8999
  1.1392
  1.9187
  0.5566
  2.8456
  1.1555
  2.6582
  1.8947
  3.0193
 -0.6510
  2.8785
  1.8772
  0.4844
  0.8612
  1.8900
  0.7892
  1.1272
  0.8738
 -1.7467
 -1.6843
  1.1454
  0.3550
 -0.1279
  1.4295
 -1.5281
  0.6976
 -0.4138
 -1.6031
  0.0238
 -0.9090
 -0.8296
  0.0698
 -0.7664
  0.1257
  0.0856
 -0.7663
 -0.1730
  0.0501
 -3.8309
 -0.9095
  1.5877
 -2.0618
  1.1629
 -0.5346
  1.4843
 -0.4801
  0.3676
  1.2108
  0.1150
  2.2080
  0.1186
 -0.9722
  2.8580
  3.1546
  0.2030
 -0.2574
  0.7604
  1.4307
  1.9173
  1.9908
 -0.4139
  2.0059
  1.5550
  0.9829
  1.7817
  1.2896
  0.9617
 -0.4163
  2.0566
 -1.1827
  0.1315
  0.6947
  2.3832
  0.9718
  2.9559
  2.1636
  1.4657
 -0.4206
  1.0273
  1.6162
  0.1423
  0.8409
  1.2249
  0.1738
 -0.7658
  1.4246
 -0.5914
  0.9685
 -0.2473
  0.6403
  0.1202
  2.8305
  1.8404
  0.7195
  0.8450
  0.8665
  0.3925
  1.2819
  1.2943
 -0.0692
  0.7312
 -0.5107
  2.2792
  1.3112
  1.0806
  0.8045
 -1.4517
 -0.2039
 -0.9176
  0.0259
  1.9817
  0.8717
  1.9949
  1.2515
  1.4779
 -0.7946
  0.7333
  1.1782
 -0.5233
 -0.8781
 -0.4868
 -0.0270
  1.1296
  1.9269
 -0.9078
 -1.3374
  1.3696
 -0.6361
  0.9839
 -0.5358
  0.2322
  1.2186
  0.1677
  0.6237
 -1.3246
 -0.3518
  0.2687
  1.9992
 -0.9640
 -2.0997
  0.4358
  0.3715
 -0.3563
  1.2112
  1.7338
 -0.3579
  0.7963
  0.8509
  0.6421
 -0.7729
  1.4243
  0.0560
  0.7998
 -0.4061
  1.3181
  1.7719
 -1.1845
  2.0268
  1.3099
  0.5721
  0.0975
  1.6786
  0.9948
  0.7494
  1.3049
 -1.3067
  1.1929
  0.6630
 -0.4178
  0.1025
  0.7211
  1.5472
  1.9725
 -0.0345
 -0.6443
  0.2606
  2.9546
  0.7146
  0.6745
  2.6818
  1.5451
  0.9759
  1.2438
  0.1962
  0.9226
  1.2769
 -0.9054
  2.0040
  0.8901
  1.7423
  2.9928
  0.9216
  1.4177
 -0.2905
 -1.6977
 -0.2272
  1.0012
  0.2904
  0.2847
  1.2466
  0.5009
  0.5594
 -0.1132
  0.7887
  0.2414
  0.7482
 -0.4766
  1.7148
  0.2320
 -1.1026
  1.4615
  1.3828
 -0.2220
  2.2326
  0.1381
  0.4174
 -0.1928
 -1.5015
  2.1823
 -1.0653
  1.2535
 -0.7977
  1.3084
 -1.9015
 -0.0804
 -2.2161
  1.0539
  0.1006
  0.0073
  0.2243
 -1.5358
 -1.1566
  0.4239
  0.1690
 -2.0356
 -1.3256
 -1.2699
 -0.7795
 -0.5868
  0.8906
 -1.7456
  0.0718
  0.4217
 -0.9233
 -0.0315
  2.0930
  1.3776
  1.2565
  1.4343
  1.7304
  1.2928
  0.1174
  2.0467
  2.0676
 -0.6263
  1.4146
  1.4532
 -0.2551
 -0.2052
 -0.8354
 -0.1847
 -2.1156
  0.3712
 -0.5154
 -1.5088
 -0.9444
  0.4847
  0.2197
 -0.2735
  1.4136
  1.3569
 -0.3926
  2.3640
  2.0387
  0.5395
  0.6080
  2.3219
  2.1205
 -0.2748
  1.3583
  0.8915
 -0.0852
  0.1396
 -0.1253
  0.3358
  0.6581
 -1.5193
  0.2624
  0.8750
  0.5877
  1.8113
  1.7156
 -0.1321
  0.9006
  1.8299
  1.0458
  0.1625
 -0.6689
 -0.5542
  0.4071
  0.5860
  1.2610
  1.5487
  1.9484
  1.2479
  0.1724
  1.1571
  2.3641
  3.2412
  0.4163
  1.1510
  1.3440
  2.0566
  1.7983
  1.3582
  2.7131
  0.5872
 -1.8136
 -0.5520
 -0.6063
 -0.7519
  2.2885
 -2.7624
 -1.1484
 -0.5734
 -1.2779
  0.4022
 -0.3527
  2.9046
  0.7530
 -0.3082
  1.6597
 -1.0119
 -1.0430
 -0.0119
  0.5274
  1.8101
 -0.1875
  0.1456
  0.4268
  1.1257
  0.4969
 -0.6699
 -0.3003
  2.6721
 -0.7229
  1.8091
  1.7066
  0.3501
  0.5784
  2.1265
 -0.0519
  1.0547
  3.0125
  0.9430
 -0.0821
  1.7452
 -0.9120
 -1.4981
  0.0498
  0.0212
  0.9317
  2.5221
 -0.5227
 -1.2694
 -0.5955
 -0.9830
  1.2788
  0.0168
  1.4585
  0.5513
  0.2314
  1.4177
  1.2202
 -0.8995
  1.2640
  1.7321
  1.8585
  1.3984
 -0.8614
 -0.7038
 -0.2470
  0.1571
 -0.0277
 -0.2841
 -0.3932
  0.0244
 -0.0089
  0.1567
  1.1270
  2.0694
  0.1430
  3.2303
 -0.6483
 -0.4757
  0.8070
  3.1268
  2.9646
  1.5943
  1.0870
 -0.1138
  1.6221
  1.7954
  0.7984
  0.2214
  2.1570
  3.5287
 -0.9524
  2.4734
  1.2707
  1.2699
  1.3358
 -0.6994
  2.7280
  0.6274
  2.7493
 -0.7191
 -0.4697
  0.4051
  1.3262
 -0.6113
 -0.7692
 -0.3166
 -0.5609
  1.1756
 -1.6216
 -0.8135
 -1.1552
 -0.6612
 -0.1363
 -0.2526
  1.1955
 -1.4716
 -0.4652
  1.7492
 -0.0201
 -0.0116
 -1.4242
 -0.4995
 -0.3115
 -1.2125
  0.5780
 -0.0437
 -1.0758
  0.1629
 -1.5336
  0.2630
  1.5366
  2.4515
  1.7889
  1.1686
  0.2914
  2.0238
  0.1669
 -0.1130
  0.6544
  1.5652
 -0.8881
 -1.7818
 -0.4497
 -1.1852
  0.5483
  1.1342
 -0.3453
 -0.0073
  0.4535
 -0.8098
  0.9184
  2.3778
  2.4125
  0.4540
  1.2443
  0.1279
  1.8489
  2.0308
  0.3099
  0.4955
 -0.6419
  0.9874
  1.5724
  1.9875
 -0.2719
 -0.5782
 -0.5712
 -0.1881
  0.1067
  0.6983
  2.4521
  1.0355
 -0.7590
 -0.0025
  1.8534
  3.1972
  1.4240
  0.5526
  1.3436
  2.0069
  0.3319
  0.5100
  1.8645
  0.7004
  1.3077
  0.1095
  0.2798
  1.4472
  0.7224
  1.8775
  0.8817
  0.5132
  0.6116
  1.6772
  1.1585
  3.2538
  2.3467
  1.2920
 -0.8726
  0.4438
  0.2226
 -0.7910
  0.6223
  1.1469
  0.2907
 -1.4775
  2.2753
 -0.1665
 -0.6647
 -0.4018
 -0.3716
 -0.9036
  2.0191
 -0.5720
  0.2666
  0.0533
 -0.5168
 -0.9268
  1.6993
 -0.5090
 -1.6446
 -0.2711
 -0.4766
 -0.2625
  0.1989
 -1.8224
  0.1569
  0.4491
  0.7869
  2.0635
  0.2618
  1.2936
  0.5267
  1.6318
  2.0094
  1.3112
  1.4698
  1.7239
  1.5535
  0.4193
 -1.2681
 -0.3643
  0.4154
  0.6628
  1.2838
 -1.1801
  0.1762
  1.7080
 -0.9267
 -1.9458
  0.3239
  2.9242
 -1.0681
  0.7165
  1.4305
  2.0130
  2.9180
 -1.1842
 -0.5184
 -0.4488
 -0.5468
 -1.1370
 -0.8240
 -1.4765
 -0.0657
 -0.0152
 -0.5404
  0.0706
  0.7630
  2.0074
  1.6130
 -0.5120
  2.0299
  0.7063
  0.4875
 -0.6537
  2.7704
  1.0621
  1.2433
  1.2571
  1.1764
  1.1771
  0.5755
  1.9833
 -0.8436
  1.7422
  1.4819
  0.0651
  0.1537
  0.9760
  0.8924
  2.1472
 -0.2993
  1.8130
  1.1245
  1.2744
  1.4575
  2.3975
  0.0350
  0.3960
  2.1802
 -0.6564
  0.1246
  0.3150
 -0.8214
  0.5516
  0.1072
  1.4430
 -0.4585
 -0.7405
 -2.2627
 -1.9521
 -0.6212
 -0.2661
 -0.8953
  2.5564
 -0.3373
 -1.2389
  0.1745
 -0.7743
  0.7624
  1.1186
 -1.0424
  0.6071
  2.2983
 -0.1568
 -0.2525
  0.2871
  2.1048
  0.6880
  0.7486
 -1.1339
  1.5712
  0.4554
  2.0859
  0.9036
  1.1834
  1.7656
 -0.2090
  1.0019
  0.7163
  0.2709
 -0.2506
  0.8893
  2.0626
  0.0921
 -0.0846
  0.2871
 -1.2615
  2.0643
 -0.2493
  1.1729
  1.4805
  1.4131
  2.2189
  1.8079
  2.2696
  0.6160
 -1.5983
 -0.6326
  1.3534
 -0.8582
 -0.0815
 -0.3097
 -0.9566
 -0.2276
 -0.2405
  0.7082
 -1.8349
  1.1804
  0.6494
  0.3812
  2.1248
  1.9831
  1.0739
  0.1025
  0.0049
 -0.2311
  0.9824
  1.7964
 -0.9846
  0.1136
  1.9111
  1.0937
  1.3608
  1.0490
  2.0613
  1.4478
  1.6952
  1.4675
 -0.1030
  2.9059
  1.1954
 -0.2487
  1.0003
  2.1038
  1.5679
  2.4823
 -1.1418
  1.3470
  0.4516
  0.6959
 -0.3821
 -2.6844
 -1.6384
  0.4196
  1.3446
  0.9196
 -0.1827
 -0.1170
 -0.4698
 -1.1287
  1.4713
 -1.7572
  1.0908
  1.6243
 -0.4250
  1.1985
 -0.3019
  0.6648
 -0.2334
  2.9474
  0.6299
  0.3328
  0.7553
  1.5728
 -0.4991
 -0.2874
 -0.0796
  1.2476
  0.6851
 -0.0693
  1.1194
  1.4506
  1.3895
  3.2225
  1.2900
  2.2048
  0.9049
 -0.2550
  0.0207
  0.5384
 -2.0797
  0.1939
  1.0153
  0.2155
 -0.5543
  1.0344
  0.3880
  1.7404
  0.2201
  1.5336
  0.1880
  0.9953
  1.1569
  0.9228
  0.8792
  0.4196
 -1.7004
 -0.5576
  0.3069
 -0.9890
  0.1553
  2.4138
  0.4023
  0.2618
 -0.8701
  1.4772
  0.5117
  0.1192
  0.2058
  1.3298
  1.7365
 -0.8802
 -0.0480
 -0.2624
  0.0831
  1.0553
  0.6970
  1.4467
  1.9340
 -0.5254
  0.7580
  1.0710
 -0.3951
  1.5098
  0.7413
  0.3566
  0.6969
  0.3524
  0.0357
  0.9793
  1.5768
  0.3410
 -0.2764
  1.7445
 -0.1505
  0.6216
 -0.0821
  0.6125
 -1.4823
  1.4641
  0.4071
                                                                                                                ./._simimg.nei                                                                                      000755  000765  000765  00000000312 11541743341 013066  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      simimg.nei                                                                                          000755  000765  000024  00000000103 11541743341 013002  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         -1 1
-1 1
 0.20  0.80  0.20 
 0.80  0.00  0.80 
 0.20  0.80  0.20 
                                                                                                                                                                                                                                                                                                                                                                                                                                                             ./._simimg.cr                                                                                       000755  000765  000765  00000000312 11541743340 012716  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      simimg.cr                                                                                           000755  000765  000024  00000047040 11541743340 012645  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
1
1
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
1
1
1
1
1
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
1
1
1
1
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
2
2
2
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
2
1
1
1
1
1
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
1
1
1
1
1
1
1
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
1
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
1
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
1
1
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
1
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
2
1
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
1
2
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
1
1
1
1
1
1
1
1
1
2
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
1
1
1
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
1
1
2
2
2
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
2
2
2
2
2
2
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
2
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
1
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
2
1
1
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
2
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
1
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
2
2
2
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
2
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
2
2
1
2
2
2
2
2
2
1
1
1
1
1
1
1
1
1
1
1
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
2
1
1
1
1
1
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ./._simimg.cf                                                                                       000755  000765  000765  00000000312 11541743337 012710  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      simimg.cf                                                                                           000755  000765  000024  00000047205 11541743337 012642  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 
2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 
2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 
2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 
2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 
2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 
2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 
1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 
2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 
2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 1 1 1 1 
1 1 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2 2 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 1 1 1 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 2 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 2 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 
2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 
1 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 1 2 2 2 2 
2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 1 1 1 2 2 2 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 2 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 
2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 
2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 
1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 

                                                                                                                                                                                                                                                                                                                                                                                           ./._simimg.mf                                                                                       000755  000765  000765  00000000312 11541743341 012715  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      simimg.mf                                                                                           000755  000765  000024  00000000442 11541743341 012637  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         Command line : 

  ./nem_exe simimg 2 

Criteria U=NEM, D=Hathaway, L=mixture, M=markov ps-like, error

  -1872.97    -19456.8    -15440    -21236.1   nan

Beta (fixed)
  1.0000
Mu (1), Pk, and sigma (1) of the 2 classes

     -0.067     0.5     0.983476 
       1.06     0.5     0.983476 
                                                                                                                                                                                                                              ./._nem_user.txt                                                                                    000755  000765  000765  00000000312 11541743335 013465  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      nem_user.txt                                                                                        000755  000765  000024  00000053775 11541743335 013430  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                                                 ==========
                        NEM Manual
                        ==========

%PrintHelpGeneral

Goal
====

    This program computes a partition of a given set of objects 
    described by one or several numeric variables and by their 
    spatial relationships, using the 'Neighborhood EM' algorithm 
    (NEM). This algorithm is derived from the EM algorithm applied
    to a hidden Markov random field model. Its new feature consists in
    taking into account some spatial interdependance between the 
    objects.

    It may be used for:

    - unsupervised segmentation of color or gray-level images 
      (points = pixel values, geographic position = pixel coordinates) ;

    - clustering of spatial data like socio-economical activities of
      neighbouring counties, etc.

    The algorithm takes as input an objects-variables table, 
    and a specification of the neighborhood relationship between 
    the objects. It produces as output a fuzzy or a hard partition 
    of the objects.  The main algorithm is described in the following
    paper:

Ambroise, C., Dang, V.M. and Govaert, G. (1997). Clustering of spatial
  data by the EM algorithm, in A.~Soares, J.~G\'omez-Hernandez and
  R.~Froidevaux, eds, `geoENV I - Geostatistics for Environmental
  Applications', Vol. 9 of `Quantitative Geology and Geostatistics', 
  Kluwer Academic Publisher, pp.~493--504.


Changing default parameters
===========================

    The behaviour of this clustering algorithm can be adjusted in many
    ways to fit a particular problem. The main possibilities
    are described below.  

    - The assumed degree of spatial interdependance is controlled by
    the value of the 'beta' coefficient (option '-b beta_value'). The higher
    it is, the smoother the partition will look in the geographic
    space, but the less it will fit to the data. The default value (1.0)
    seems to work well in most image segmentation problems where the
    patches are supposed to be spatially smooth. For beta's lowest
    value (0.0), the algorithm is the same as Dempster et al's EM
    algorithm (1974), and does a 'spatially blind' segmentation.

    - The algorithm 'Neighborhood Classification EM' (NCEM) can be used
    instead of NEM (option '-a ncem'). The principle of NCEM consists in
    'hardening' the classification matrix at each iteration (C-step
    after the E-step). Practically, NCEM converges faster than NEM,
    but gives a poorer segmentation on data containing a high level of
    noise.  

    - In some applications, the class may be already known for a part of
    the sample. Such a knowledge can be taken into account by the
    Neighborhood EM algorithm, and may improve considerably the resulting
    classification (option '-s l file.ck').  

    - Incomplete observations are taken into account in the
    probabilistic model and the program.  It is simply assumed that
    the missingness occurs at random, i.e. it does not depend on the
    missing value itself nor on the unobserved class.


Notice
======

    The program is only provided to make it easier to test the
    behavior of the Neighborhood EM clustering algorithm and compare
    it with other algorithms.    Although I have tried to write and 
    test the program as carefully as possible, it is not guaranteed 
    to be error-free.  Please contact me if you have
    any questions or problems in using it.  Please also mention its
    origin if you use it for a published work. Finally I would be
    interested to know for what kind of problem you have found
    this program to be of use.

Van M� Dang
Van.Mo.Dang@utc.fr
http://www.hds.utc.fr/~mdang


%end PrintHelpGeneral

%PrintHelpOptions

Command Syntax
==============

 Usage :    nem_exe   file  K  [ option1 option2 ... ]
 ------- 

 Arguments :
 -----------
   file       base name of input files ___.str and ___.dat
   K          number of classes

 Options :   [ default  ] { possible values }
 --------- 

 
  -a algo   [ nem      ]   { nem ncem gem }
     Algorithm to compute classification at E-step of each iteration : 
      ncem = crisp classification by ICM procedure
      nem  = fuzzy classification by mean field approximation
      gem  = fuzzy classification by Gibbs sampling (Monte-Carlo simulations)

  -b beta   [ 1        ]   (0.0 <= beta <= 4)
     Coefficient of spatial smoothing to apply. This matches the 
     Potts random field strength of interaction with 4-neighbor contexts.
     Notice that b = 0 is equivalent to EM for a mixture model.

  -c wh thr [clas  0.04]   { none clas crit } and (>0)

       none     = no convergence test, i.e. do all specified iterations.

       clas thr = stop the iterations when the largest difference
                  between previous and current classification matrix
                  is <= threshold. A threshold 0.04 is usually optimal.

       crit thr = stop the iterations when 
                  | (current_crit - last_crit)/current_crit | < threshold.
                  A threshold of 0.001 is usually best.  The test uses 
                  the criterion selected by the option -C.

  -f format [ hard     ]   { hard fuzzy }
     Format of output partition. Hard = N integers having values from
     1 to nk : for each observation, give the number of the class where
     it has highest grade of membership. Fuzzy = N x K reals between
     0 and 1 : for each observation, give its grade of membership in
     each class.

  -i itmax  [ 100      ]   (>= 0)
     Maximum number of NEM iterations.

  -l dolog  [ n        ]   { y n }
     Produce a log file or not to see the results of each iteration.

  -m f p d  [norm p_ s__]  { norm lapl } { p_ pk } { s__ sk_ s_d skd }

     Mixture model assumption to use

      norm/lapl/bern : normal, Laplace or Bernoulli distributions
                       Bernoulli distributions are for binary data

      p_/pk     : clusters have equal / varying proportions

      s__/...   : variance model
                    s__ : same variance in all clusters and variables
                    sk_ : one variance per cluster, same in all variables
                    s_d : one variance per variable, same in all clusters
                    skd : one variance per cluster and variable

                  the variables are assumed independent within a cluster

  -n neigh  [ 4        ]   { 4 f }
     Neighborhood specification to use in the case of an image.
     4 = default 4-nearest neighbor system. f = specify neighborhood
     window in file.nei.

  -o fout   [ file     ]
     Output files basename ___.cf, ___.mf and ___.log. Default is 
     to use input file basename. Specify '-' to output the
     classification to standard output; useful to pipe the result to an 
     'nem_exe -s f -' session.

  -s init   [ s 1      ]   { s <v> | f <ini.uf> | r <n> | l <file> } 
     Initialization mode.
     -s s <v>   
        Sort observations by variable <v>, then divide them
        in K quantiles of equal size to get initial partition.

     -s f <ini.uf>
        Read initial fuzzy classification from file <ini.uf>. 
        '-s f -' reads from standard input.

     -s r <n>
        Start <n> times from random parameters (means chosen at
        random among the observations), then keep result with highest
        criterion.

     -s l <file>
        Use partially known labels given in <file> to compute initial 
        parameters. Those labels remain fixed throughout the
        clustering process.

     -s mi <para> 
     -s mf <para>
        Use specified parameters at beginning (mi) or throughout the
        clustering process (mf). Parameters syntax :
        p_1 ... p_{K-1}   m_11 .. m_1D m21 ... m_KD   s_11 ... s_KD.
        The s_kd are the standard errors for normal distributions, or
        the scale parameters for the Laplace distributions.

  -t tie    [ random   ]   { random first }
     How to choose the class with highest probability when several
     classes have same maximum probability in MAP
     classification. 'random' draws uniformly between ex-aequo
     classes, 'first' chooses class with lowest index.
 

  -B bmod   [ fix      ]   { fix psgrad heu_d heu_l }
     Procedure to estimate beta automatically :
      fix   = no estimation of beta, use beta given by option '-b'
      psgrad = pseudo-likelihood gradient ascent
      heu_d = heuristic using drop of fuzzy within cluster inertia
      heu_l = heuristic using drop of mixture likelihood

  -C crit   [ U        ]   local maximum criterion { U M D L }
     Criterion used to select the best local solution from random starts :
      U = fuzzy spatial clustering criterion
      M = fuzzy pseudo-likelihood
      D = fuzzy within cluster inertia
      L = likelihood of mixture parameters

  -G nit conv step rand [  1 0.001 0.0 0 ]  
     Parameters of beta gradient estimation
      nit = number of gradient iterations
      conv = threshold to test convergence (|g'|<conv*N is tested)
      step = > 0 for fixed step, <= 0 for Newton step = 1/g''
      rand = in -s r  init mode, initial beta random (1) or fixed by -b (0)

  -H bstep bmax ddrop dloss lloss [ 0.10 2.0 0.8 0.5 0.02 ]  
     Parameters of beta D and L heuristics :
      bstep = step of beta increase
      bmax  = maximal value of beta to test
      ddrop = threshold of allowed D drop (higher = less detection)
      dloss = threshold of allowed D loss (higher = less detection)
      lloss = threshold of allowed L loss (higher = less detection)

  -I eiter  [ 1        ]   (>= 1)
     Number of internal E-step iterations (nem and ncem algorithms),
     i.e. number of sweeps through whole dataset to compute
     classification at each iteration.  For gem algorithm, indicates
     number of sweeps through the dataset to compute the average
     frequency of class occurrence -> a large value is recommended
     for the gem algorithm (typically 50).

  -M miss   [ replace  ]    { replace ignore }
     How to deal with missing data.  Replace by expected value (EM) or
     ignore when computing mean (maximize fuzzy clustering criterion).

  -O order  [ direct   ]   order of site visit { direct random }
     Order in which to classify the observations at E-step.
      direct = given order 1..N
      random = random permutation of 1..N

  -S seed   [ <time>  ]    (integer)
     Specify a seed for the random number generator. Default uses
     current system clock.

  -T test   [ n       ]    { y n }
     Print some debugging information or not.

  -U update [ seq      ]   { seq para }
     Update the class of the sites in a sequential or parallel manner.
     'seq' works best. 'para' is more grounded, because it is EM with
     mean field approximation, but it requires a low spatial smoothing.

You may also just type arguments : 

  -v                      versions information
  -h help_topic           longer help - help topics are
     general
     options
     examples
     filein
     fileout
     versions
 
 
%end PrintHelpOptions

%PrintHelpExamples

    Examples :
    ----------
 nem_exe  myimg  3 -b 0.5 -s r 10 -o Res/myimg3r_05 >&! Res/myimg3r_05.out &

    This unix command clusters data set myimg (files myimg.dat, myimg.str)
    into 3 classes ; spatial coefficient is 0.5 ; do 10 random starts ;
    save results in files Res/myimg3r_05.* (___.cf ___.log ___.mf) .
    Last part of the command is unix-specific : it saves screen output
    to file Res/myimg3r_05.out and executes the program in background.

%end PrintHelpExamples


%PrintHelpFileIn

Input Files
===========

    Input files are in ASCII format.
    2 or 3 input files are required : file.str   file.dat   [ file.nei ]
    2 optional input files :          file.u0    file.ck

 1) file.str
 -----------
    This gives the structure of the data : 
    type of spatial repartition (image, spatial or non-spatial),
    number of objects and variables. This file may start with
    comment lines to describe the dataset (lines beginning with #).

    - Ex 1 : Color image (each pixel is described by 3 variables) 
             of 200 lines (height) and 300 columns (width)
    # RGB biomedical coloscopic image. Look for 3 or 4 classes.
    I 200 300 3


    - Ex 2 : Spatial dataset of 3000 objects described by 4 variables
    # Economic data on counties (region of Centre). About 5 classes.
    S 3000 4


    - Ex 3 : Non-spatial dataset of 3000 objects described by 4 variables
    # Companies characteristics for risk assessment. About 6 classes.
    N 3000 4


 2) file.dat
 -----------
    Contains the objects-variables table (only the non-spatial
    variables).  Missing data are specified by NaN.


    If the dataset is an image, the pixels must be listed 
    line-by-line first, i.e. : 

    x(1,1) x(1,2) ... x(1,nc)  
    x(2,1) x(2,2) ... x(2,nc)
    ...
    x(nl,1) x(nl,2) ... x(nl,nc)

    where nl = number of lines, nc = number of columns, and 
    x(i,j) = values of pixel at line i and column j (a set of 3
    numbers for a color image).

    - Ex 1 : Color image
    50 100 120
    51  99 122
    ...

    - Ex 2 : Spatial data (4 variables)
    0.31 200 41 1200 
    0.28 202 43 1180
    ...


 3) file.nei
 -----------
    Specifies the spatial relationships between the objects. 
    * For an image, this file is optional and allows to specify
    a particular neighborhood system ; if no file is specified,
    default neighborhood system taken is 4-nearest neighbours.
    * For other spatial data, this file is required.

    Format is different for an image and other spatial data
    (see examples below).

    - Ex 1 : 4 nearest-neighbours in an image (default)
    -1 1           /* at most 1 pixel on the left and 1 on the right */
    -1 1           /* at most 1 pixel up and 1 down */
    0 1 0          
    1 0 1          /* 4 equally weighted neighbors : */
    0 1 0	   /* up,left,right,down */            

    - Ex 2 : other spatial data, unweighted neighborhood graph
    0              /* 0 = no weight specified (use default weight = 1.0) */
    1   3   2 5 7  /* object 1 has 3 neighbors : objects 2, 5 and 7 */
    2   2   1 3
    ...


    - Ex 2 : other spatial data, weighted neighborhood graph
    1                           /* 1 = weights are specified */
    1   3   2 5 7  0.5 0.6 0.8  /* object 1 has 3 neighbors : 2, 5 and 7 */
    2   2   1 3    0.3 0.5
    ...


 4.a) file.ck  (option -s l file.ck)
 ------------
    Gives the class of the points for which the label is already
    known. Expected format is N + 1 integers (N being the total
    number of objects in the sample) :
    - the first number is the number of classes K ;
    - the N following integers, in the same order as the
    objects in file.dat, indicate the class to which the
    corresponding object belongs.  If the value is 0 or greater than the
    number of classes, then the object is considered to have no known
    label (its label will be computed by the algorithm). 

    In the current implementation, each class must contain at least
    one observation with known label. Those pre-labeled observations
    are used to initialize the centers of the clusters.

    5     /* number of classes */
    0     /* object 1 has no known label */
    5     /* object 2 belongs to class 5 */
    1     /* object 3 belongs to class 1 */
    ...


 4.b) file.u0  (option -s f)
 ------------
    Gives an initial fuzzy classification to start the algorithm. 

    0.9 0.1  /* object 1 : initial membership = 0.9/0.1 in class 1/2 */
    0.3 0.7  /* object 2 : initial membership = 0.3/0.7 in class 1/2 */


%end PrintHelpFileIn



%PrintHelpFileOut

Output Files
============

    Output files are in ASCII format.
    2 files are output :              file.cf    file.mf
    1 optional file is output :       file.log


 1) file.cf (or file.uf for fuzzy partition)
 -------------------------------------------
    Gives the partition found by the algorithm.

    - Ex 1 : Image segmented in 2 'hard' classes (linewise order, as the 
           input data file)
    1        /* Object 1 belongs to class 1 */
    1        /* Object 2 belongs to class 1 */
    2        /* Object 3 belongs to class 2 */
    1
    2
    2
    ...

    - Ex 2 : Data segmented in 3 fuzzy classes
    0.2  0.7  0.1    /* Object 1 : class 1 with probability 0.2 , etc. */
    0.8  0.05 0.015
    ...


 2) file.mf
 ----------
    
    This file gives the values of criteria optimized by the algorithm,
    and the parameters of the mixture calculated by the
    algorithm (means, proportions, scale parameter). Example: 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~ START EXAMPLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Command line : 
 
  nem_exe intro 2 -b 1.00 -f fuzzy -R intro.cr -o intro_b10 -s r 10 -C M 
 
Criteria U=NEM, D=Hathaway, L=mixture, M=markov ps-like, error

  56.3807    -2671.73    -1788.76    -2793.36   0.335938
 
Beta (fixed)
  1.0000
Mu (4), Pk, and sigma (4) of the 2 classes
 
   0.947  0.00456  0.0199  0.998   0.5    1.15226  1.15226  1.15226  1.15226 
   0.0614 1.03     0.864  -0.163   0.5    1.15226  1.15226  1.15226  1.15226 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END EXAMPLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



 3) file.log
 -----------
    Details each iteration results (optimized criterion and 
    class parameters).


%end PrintHelpFileOut


%PrintHelpVersions

History of modifications
========================

Version 0.00  (01.02.1996)
------------
First version released on WEB. Initialization by sorting variable.


Version 1.00  (31.05.1996)  
------------
Random initializations. Image default neighborhood system. Long help.


Version 1.01  (27.06.1996)  
------------
Added an -a ncem option to implement the crisp clustering version of NEM.

Each observation is now updated sequentially in turn. The previous 
parallel updating would produce chessboard-like images for high betas.


Version 1.02  (17.10.1996)
------------
Added the possibility to take into account a partial knowledge
of the classification into the clustering procedure 
(option '-s l <file.ck>').


Version 1.03  (02.10.1997)  
------------
If a partial knowledge of the classification is available, the
intial cluster centers are computed from the observations with
known labels.  

The log file is made optional (option '-l y').  The ___.mf file now
also contains the estimated cluster proportions, volumes and
covariance matrices and the values of the final criteria.


Version 1.04  (11.01.1998)  
------------

Two heuristics are implemented to estimate the spatial smoothness
coefficient beta.  The first heuristic is based on detecting a sharp
drop in the fuzzy within-cluster inertia D, or a sufficient decrease
from its maximum value, when beta is slowly increased.  The second
heuristic is based on detecting a sufficient decrease of the
log-likelihood L of the mixture parameters from its maximum value.
The heuristics may be invoked with '-B heu_d' or '-B heu_l'. Their
default parameters may be changed with '-H ...'. 

The final partition may now be printed to standard ouput instead of to a
file (option '-o -').  The result can thus be redirected as an 
initial partition to another nem_exe session's input.

The initial partition may now be read from any file with option '-s f
foo.bar' (previously inputbasename.u0 was used).  In particular, the
partition may be read from standard input (option '-s f -').  This allows
to read the initial partition through a pipe from the result of a previous 
nem_exe session.

At each iteration, the fuzzy classification at the E-step may now be
computed by either of two methods :
- Neighborhood EM's fixed point technique (default, '-a nem')
- Gibbsian EM's Gibbs sampler technique ('-a gem').

A longer explanation is given for the options and successive versions
(option '-h helptopic').

Other options have been added, to test the effect of alternative
parameter values.  Those options usually do not change considerably
the default behaviour of the algorithm :

- At the E-step of each iteration, the classification may now be
  updated in a random order instead of 1..N (option '-O random').
- The seed of the random number generator may now be given (option '-S seed').
- The number of E-step internal iterations may be changed (option '-I eiter').
- The convergence threshold may also be changed (option '-c cvthres').
- Another criterion than U may be chosen to select best result ('-C crit').
- Compute the classification error in two-class case ('-R refclass').

A few internal changes were also made to make the program portable to
the djgpp gcc compiler for MS-DOS (srand48 replaced by srandom).


Version 1.05  (09-APR-1998)  
------------

This version mainly adds the capability to deal with missing data.
This means that some of the N observation vectors may be incompletely
observed.  In the 'name.dat' file, use NaN to indicate unobserved
components of an observation vector.  Two slightly different
techniques are provided to deal with missing data.  The first and
default behaviour (invoked using switch '-M replace') roughly consists
in replacing any missing component with its expected value.  This
technique implements the EM procedure and finds parameters maximizing
the likelihood.  The alternative behaviour (invoked using switch '-M
ignore') consists in ignoring missing components.  This means that the
means are computed using only observed components.  This alternative
technique finds a classification matrix and parameters which maximize
the fuzzy classifying log-likelihood.  It appears to converge a bit
faster than the 'replace' mode.

The iteration count is displayed in a more economic way now (all the
iteration numbers were displayed separately).

A new criterion is computed, the fuzzy pseudo-likelihood, named M.
Using this criterion in order to choose the best result 
may prove less sensitive to the value of beta than using U.

In the random start strategy, two initialization tactics have been
made more sensible.  The initial volumes are computed as whole volume
/ number of classes (the whole volume was used previously).  The means
are redrawn until all drawn means are different --- this avoids the
problem of artificially merging together two classes.

A few internal changes were also made in order to allow direct call to
the program from as a Matlab function.  This allowed to detect and
remove a few potential bugs that had gone unnoticed (a file not
closed, use of memory just after freeing it).


%end PrintHelpVersions
   ./._makefile                                                                                        000755  000765  000765  00000000312 11541743317 012607  0                                                                                                    ustar 00tom                             tom                             000000  000000                                                                                                                                                                             Mac OS X            	   2   �      �                                      ATTR       �   �   .                  �   .  com.dropbox.attributes   x��V*�/άP�R�VJ�HM.-IL�IsSK�%�@[[���Z @                                                                                                                                                                                                                                                                                                                      makefile                                                                                            000755  000765  000024  00000015747 11541743317 012547  0                                                                                                    ustar 00tom                             staff                           000000  000000                                                                                                                                                                         #
#	--- Main  groups of targets to be updated ---
#

# V1.04-a  04-10-1997  MD  add lint
# V1.04-b  18-11-1997  MD  replace cc with gcc, add -lemu option
# V1.04-c  11-01-1998  MD  archive name as variable
# V1.05-a  25-01-1998  MD  add exemain.c -> nem_exe
# V1.05-b  25-01-1998  MD  add mexmain.c and mex target -> nem_exe.mexsol
# V1.05-c  26-01-1998  MD  add exememo.c -> nem_exe and mexmemo.c -> .mexsol
# V1.06-a  17-06-1998  MD  change nem_noi -> nem_mod
# V1.06-b  23-06-1998  MD  flag -Wall
# V1.06-c  23-06-1998  MD  add target t_nem_arg
# V1.06-d  30-06-1998  MD  add target t_nem_mod
# V1.07-a  29-02-1999  MD  archive name : new version

CFLAGS = -Wall -D__GO32__

# Flags for MS-DOS compilation using DJGPP's gcc
#CFLAGS = -lemu -Wall

# Select C compiler here
CC = gcc -I .
#CC = cc

# Warning messages


# Name of uncompressed archive file to create
ARC = nem107.tar

#	Optimized application executables (default make)
exe :		txt2hlp nem_exe geo2nei randord err2 tcpu

#	Matlab executable
mex :		nem_exe.mexsol

#	Debuggable executables
dbg :		g_nem_exe g_geo2nei g_randord

#	Archive
arc :		$(ARC).Z

#	Everything
all :		remobj mex exe dbg arc

#	Test version of nem_exe
test :		nem_tmp

#	Tight syntax check
check :
	lint nem_exe.c nem_arg.c nem_alg.c nem_nei.c \
	nem_mod.c nem_rnd.c nem_ver.c nem_hlp.c \
	lib_io.c exemain.c exememo.c
	lint err2.c
	lint tcpu.c

#	Remove all object files
remobj :
	\rm *.o

#
#	--- Single targets ---
#

#	Archive commands to update source file archiving
$(ARC).Z :	nem_exe.c nem_arg.c nem_alg.c nem_nei.c \
		nem_mod.c nem_rnd.c nem_ver.c nem_hlp.c \
		nem_typ.h nem_arg.h nem_alg.h nem_nei.h \
		nem_mod.h nem_rnd.h nem_ver.h nem_hlp.h \
		lib_io.c  lib_io.h genmemo.h exemain.c exememo.c \
		mexmain.c mexmemo.c mainfunc.h \
		geo2nei.c randord.c \
		err2.c tcpu.c txt2hlp.c \
		simimg.str simimg.dat simimg.nei simimg.cr \
		simimg.cf simimg.mf \
		makefile
	tar cvf $(ARC) nem_*.[ch] lib_*.[ch] exe*.c mex*.c gen*h \
	err2.c geo2nei.c randord.c txt2hlp.c tcpu.c mainfunc.h \
	simimg.str simimg.dat simimg.nei simimg.cr simimg.cf simimg.mf \
	nem_user.txt makefile
	gzip -vf $(ARC)
	\cp -p $(ARC).gz /home_f/mdang/public_html/Progs/$(ARC).gz


#	Linkage commands to update the executables
nem_exe.mexsol : nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
		nem_ver.o nem_hlp.o lib_io.o mexmain.c mexmemo.c genmemo.h
	cmex nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
	nem_ver.o nem_hlp.o lib_io.o  mexmain.c mexmemo.c \
	-o nem_exe -lm -O

nem_tmp : 	nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
		nem_ver.o nem_hlp.o lib_io.o exemain.o exememo.o 
	$(CC) nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
	nem_ver.o nem_hlp.o lib_io.o  exemain.o exememo.o \
	-o nem_tmp -lm -O $(CFLAGS)

nem_exe : 	nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
		nem_ver.o nem_hlp.o lib_io.o exemain.o exememo.o 
	$(CC) nem_exe.o nem_arg.o nem_alg.o nem_nei.o nem_mod.o nem_rnd.o \
	nem_ver.o nem_hlp.o lib_io.o  exemain.o exememo.o \
	-o nem_exe -lm -O $(CFLAGS)

g_nem_exe : 	g_nem_exe.o g_nem_arg.o g_nem_alg.o g_nem_nei.o \
		g_nem_mod.o g_nem_rnd.o g_nem_ver.o g_nem_hlp.o g_lib_io.o \
                g_exemain.o g_exememo.o
	$(CC) g_nem_exe.o g_nem_arg.o g_nem_alg.o g_nem_nei.o \
	g_nem_mod.o g_nem_rnd.o g_nem_ver.o g_nem_hlp.o g_lib_io.o \
        g_exemain.o g_exememo.o \
	-o g_nem_exe -lm -g $(CFLAGS)


geo2nei :	geo2nei.o lib_io.o
	$(CC) geo2nei.o lib_io.o \
	-o geo2nei  -lm  -O $(CFLAGS)

g_geo2nei :	g_geo2nei.o g_lib_io.o
	$(CC) g_geo2nei.o g_lib_io.o \
	-o g_geo2nei  -lm  -O $(CFLAGS)

randord :	randord.o nem_rnd.o
	$(CC) randord.o nem_rnd.o \
	-o randord  -lm  -O $(CFLAGS)

g_randord :	g_randord.o g_nem_rnd.o
	$(CC) g_randord.o g_nem_rnd.o \
	-o g_randord  -lm  -g $(CFLAGS)

tcpu :	tcpu.c
	$(CC) tcpu.c -o tcpu  -lm  -O $(CFLAGS)


#	Linkage commands to update the test executables

#V1.06-c
t_nem_arg :	g_nem_arg.o g_lib_io.o g_nem_hlp.o g_exememo.o g_nem_ver.o
	$(CC) g_nem_arg.o g_lib_io.o g_nem_hlp.o g_exememo.o g_nem_ver.o \
	-o $@  -lm  -g $(CFLAGS)

#V1.06-d
t_nem_mod.out :	nem_mod.c t_nem_mod.c t_nem_mod.sh
	$(CC) g_exememo.o t_nem_mod.c -o t_nem_mod -lm  -g $(CFLAGS)
	@t_nem_mod.sh > tmp.out
	@echo "\ndiff tmp.out t_nem_mod.out :\n"
	if ( diff tmp.out t_nem_mod.out ) ; \
	then \
          echo "test of nem_mod.c OK" ; touch t_nem_mod.out ; \
	else \
          echo "\07\nWarning \07: output of test of nem_mod.c changed" ; \
#	  echo "Change t_nem_mod.out ?" ;  ; \
#	  if "$a" = "y"; then \mv tmp.out t_nem_mod.out; fi ; \
	fi
	@\rm tmp.out
	@\rm t_nem_mod

t_nem_mod :	nem_mod.c t_nem_mod.c g_exememo.o 
	$(CC) g_exememo.o t_nem_mod.c \
	-o $@  -lm  -g $(CFLAGS)
#
#	--- Secondary targets (object files) ---
#

#	Compilation commands to update the object files
nem_exe.o :	nem_exe.c nem_typ.h nem_arg.h nem_alg.h nem_rnd.h lib_io.h
	$(CC) -c nem_exe.c -O $(CFLAGS)

g_nem_exe.o :	nem_exe.c nem_typ.h nem_arg.h nem_alg.h nem_rnd.h lib_io.h
	$(CC) -c nem_exe.c -g $(CFLAGS) -o g_nem_exe.o


nem_arg.o :	nem_arg.c nem_typ.h nem_arg.h
	$(CC) -c nem_arg.c -O $(CFLAGS)

g_nem_arg.o :	nem_arg.c nem_typ.h nem_arg.h
	$(CC) -c nem_arg.c -g $(CFLAGS) -o g_nem_arg.o


nem_alg.o :	nem_alg.c nem_typ.h nem_alg.h nem_nei.h nem_mod.h nem_rnd.h
	$(CC) -c nem_alg.c -O $(CFLAGS)

g_nem_alg.o :	nem_alg.c nem_typ.h nem_alg.h nem_nei.h nem_mod.h nem_rnd.h
	$(CC) -c nem_alg.c -g $(CFLAGS) -o g_nem_alg.o


nem_nei.o :	nem_nei.c nem_typ.h nem_nei.h
	$(CC) -c nem_nei.c -O $(CFLAGS)

g_nem_nei.o :	nem_nei.c nem_typ.h nem_nei.h
	$(CC) -c nem_nei.c -g $(CFLAGS) -o g_nem_nei.o


nem_mod.o :	nem_mod.c nem_typ.h nem_mod.h
	$(CC) -c nem_mod.c -O $(CFLAGS)

g_nem_mod.o :	nem_mod.c nem_typ.h nem_mod.h
	$(CC) -c nem_mod.c -g $(CFLAGS) -o g_nem_mod.o


nem_hlp.o :	nem_hlp.c nem_hlp.h
	$(CC) -c nem_hlp.c -O $(CFLAGS)

g_nem_hlp.o :	nem_hlp.c nem_hlp.h
	$(CC) -c nem_hlp.c -g $(CFLAGS) -o g_nem_hlp.o


nem_ver.o :	nem_ver.c nem_typ.h nem_ver.h
	$(CC) -c nem_ver.c -O $(CFLAGS)

g_nem_ver.o :	nem_ver.c nem_typ.h nem_ver.h
	$(CC) -c nem_ver.c -g $(CFLAGS) -o g_nem_ver.o


nem_rnd.o :	nem_rnd.c nem_rnd.h
	$(CC) -c nem_rnd.c -O $(CFLAGS)

g_nem_rnd.o :	nem_rnd.c nem_rnd.h
	$(CC) -c nem_rnd.c -g $(CFLAGS) -o g_nem_rnd.o


lib_io.o :	lib_io.c lib_io.h
	$(CC) -c lib_io.c -O $(CFLAGS)

g_lib_io.o :	lib_io.c lib_io.h
	$(CC) -c lib_io.c -g $(CFLAGS) -o g_lib_io.o

exemain.o :	exemain.c
	$(CC) -c exemain.c -O $(CFLAGS)

g_exemain.o :	exemain.c
	$(CC) -c exemain.c -g $(CFLAGS) -o g_exemain.o

exememo.o :	exememo.c genmemo.h
	$(CC) -c exememo.c -O $(CFLAGS)

g_exememo.o :	exememo.c genmemo.h
	$(CC) -c exememo.c -g $(CFLAGS) -o g_exememo.o

geo2nei.o :	geo2nei.c lib_io.h
	$(CC) -c geo2nei.c -O $(CFLAGS)

g_geo2nei.o :	geo2nei.c lib_io.h
	$(CC) -c geo2nei.c -g $(CFLAGS) -o g_geo2nei.o


randord.o :	randord.c nem_rnd.h
	$(CC) -c randord.c -O $(CFLAGS)

g_randord.o :	randord.c nem_rnd.h
	$(CC) -c randord.c -g $(CFLAGS) -o g_randord.o

#
#	--- Auxiliary targets (automated source files) ---
#

nem_hlp.c :	nem_user.txt
	./txt2hlp nem_user.txt nem_hlp y \
	"Project NEM : display of long help" > /dev/null
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         