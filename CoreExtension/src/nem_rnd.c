/*\

    NEM_RND.C

    Programme NEM (Neighborhood EM) : routines de tirage aleatoire

    Van Mo DANG       Avril 96

Vers-mod  Date         Who  Description

1.04-a    10-OCT-1997  MD   add RandomPermutationAlgo()
1.04-b    05-NOV-1997  MD   use random/srandom instead of lrand48/srand48
1.04-c    12-JAN-1997  MD   add RandomReal()
1.08-a    05-JUL-1999  MD   use rand/srand, not random/srandom (not portable)
\*/


#include <sys/types.h>   /* time_t */
#include <time.h>        /* time() */
#include <stdlib.h>      /* srand(), rand() */

#include "nem_rnd.h"


void  RandomSeedByTime( void ) 
{
    time_t   t ;

    t = time( 0 ) ;

    srand( (unsigned) t ) ;
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

    nrandom = rand( ) ;

    /* result = ( (int) ( nrandom % span ) ) + Mini  ; V1.08-a*/
    result = (int) ((double)nrandom / ((double)RAND_MAX + 1) * span) + Mini ;

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

    nrandom = rand( ) ;

    result = ( (float) nrandom / RAND_MAX ) * span + Mini  ;

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

int main( int argc, char *argv[] )
{
  int       n ;
  int       i ;

  if ( argc < 4 )
    {
      printf( "Au moins 3 args\n" );
      return 2 ;
    }
  n = atoi( argv[ 1 ] ) ;

 {  int       mini, maxi ;
  mini = atoi( argv[ 2 ] ) ;
  maxi = atoi( argv[ 3 ] ) ;

  RandomSeedByTime( ) ;

  for ( i = 0 ; i < n ; i ++ )
    {
      int r = RandomInteger( mini , maxi  ) ;

      Rprintf(  "%8d  ", r ) ;
    }
 }

 {  float       mini, maxi ;
  mini = atof( argv[ 2 ] ) ;
  maxi = atof( argv[ 3 ] ) ;

  RandomSeedByTime( ) ;

  for ( i = 0 ; i < n ; i ++ )
    {
      float r = RandomFloat( mini , maxi  ) ;

      Rprintf(  "%8f  ", r ) ;
    }
 }

  Rprintf(  "\n" ) ;
  return 0 ;
}
*/
