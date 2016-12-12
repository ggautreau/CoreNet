void  RandomSeedByTime( void ) ;
int   RandomInteger( int Mini, int Maxi ) ;
float   RandomFloat( float Mini, float Maxi ) ;      /*V1.05-a*/
void  RandomPermutationAlgo( int* TabV , int Nb ) ;  /*V1.04-a*/

/* external library random functions */ /*V1.08-a suppressed for portability */
#if 0
#if !defined(__GO32__) /* djgpp has its own random */
int srandom( unsigned seed );
long random();
#endif
#endif
