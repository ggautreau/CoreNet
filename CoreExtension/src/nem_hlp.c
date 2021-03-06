/*\

    nem_hlp.c

    Project NEM : display of long help

    Thu Jul 15 09:52:40 1999

\*/

#include "nem_hlp.h"  /* Exported prototypes */
#include <R.h>  /* Exported prototypes */

void PrintHelpGeneral( FILE* F )
{
    Rprintf( "\n" ) ;
    Rprintf( "Goal\n" ) ;
    Rprintf( "====\n" ) ;
    Rprintf( "\n" ) ;
    Rprintf( "    This program computes a partition of a given set of objects \n" ) ;
    Rprintf(  "    described by one or several numeric variables and by their \n" ) ;
    Rprintf(  "    spatial relationships, using the 'Neighborhood EM' algorithm \n" ) ;
    Rprintf(  "    (NEM). This algorithm is derived from the EM algorithm applied\n" ) ;
    Rprintf(  "    to a hidden Markov random field model. Its new feature consists in\n" ) ;
    Rprintf(  "    taking into account some spatial interdependance between the \n" ) ;
    Rprintf(  "    objects.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    It may be used for:\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - unsupervised segmentation of color or gray-level images \n" ) ;
    Rprintf(  "      (points = pixel values, geographic position = pixel coordinates) ;\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - clustering of spatial data like socio-economical activities of\n" ) ;
    Rprintf(  "      neighbouring counties, etc.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    The algorithm takes as input an objects-variables table, \n" ) ;
    Rprintf(  "    and a specification of the neighborhood relationship between \n" ) ;
    Rprintf(  "    the objects. It produces as output a fuzzy or a hard partition \n" ) ;
    Rprintf(  "    of the objects.  The main algorithm is described in the following\n" ) ;
    Rprintf(  "    paper:\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Ambroise, C., Dang, V.M. and Govaert, G. (1997). Clustering of spatial\n" ) ;
    Rprintf(  "  data by the EM algorithm, in A.~Soares, J.~G\'omez-Hernandez and\n" ) ;
    Rprintf(  "  R.~Froidevaux, eds, `geoENV I - Geostatistics for Environmental\n" ) ;
    Rprintf(  "  Applications', Vol. 9 of `Quantitative Geology and Geostatistics', \n" ) ;
    Rprintf(  "  Kluwer Academic Publisher, pp.~493--504.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Changing default parameters\n" ) ;
    Rprintf(  "===========================\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    The behaviour of this clustering algorithm can be adjusted in many\n" ) ;
    Rprintf(  "    ways to fit a particular problem. The main possibilities\n" ) ;
    Rprintf(  "    are described below.  \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - The assumed degree of spatial interdependance is controlled by\n" ) ;
    Rprintf(  "    the value of the 'beta' coefficient (option '-b beta_value'). The higher\n" ) ;
    Rprintf(  "    it is, the smoother the partition will look in the geographic\n" ) ;
    Rprintf(  "    space, but the less it will fit to the data. The default value (1.0)\n" ) ;
    Rprintf(  "    seems to work well in most image segmentation problems where the\n" ) ;
    Rprintf(  "    patches are supposed to be spatially smooth. For beta's lowest\n" ) ;
    Rprintf(  "    value (0.0), the algorithm is the same as Dempster et al's EM\n" ) ;
    Rprintf(  "    algorithm (1974), and does a 'spatially blind' segmentation.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - The algorithm 'Neighborhood Classification EM' (NCEM) can be used\n" ) ;
    Rprintf(  "    instead of NEM (option '-a ncem'). The principle of NCEM consists in\n" ) ;
    Rprintf(  "    'hardening' the classification matrix at each iteration (C-step\n" ) ;
    Rprintf(  "    after the E-step). Practically, NCEM converges faster than NEM,\n" ) ;
    Rprintf(  "    but gives a poorer segmentation on data containing a high level of\n" ) ;
    Rprintf(  "    noise.  \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - In some applications, the class may be already known for a part of\n" ) ;
    Rprintf(  "    the sample. Such a knowledge can be taken into account by the\n" ) ;
    Rprintf(  "    Neighborhood EM algorithm, and may improve considerably the resulting\n" ) ;
    Rprintf(  "    classification (option '-s l file.ck').  \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Incomplete observations are taken into account in the\n" ) ;
    Rprintf(  "    probabilistic model and the program.  It is simply assumed that\n" ) ;
    Rprintf(  "    the missingness occurs at random, i.e. it does not depend on the\n" ) ;
    Rprintf(  "    missing value itself nor on the unobserved class.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Notice\n" ) ;
    Rprintf(  "======\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    The program is only provided to make it easier to test the\n" ) ;
    Rprintf(  "    behavior of the Neighborhood EM clustering algorithm and compare\n" ) ;
    Rprintf(  "    it with other algorithms.    Although I have tried to write and \n" ) ;
    Rprintf(  "    test the program as carefully as possible, it is not guaranteed \n" ) ;
    Rprintf(  "    to be error-free.  Please contact me if you have\n" ) ;
    Rprintf(  "    any questions or problems in using it.  Please also mention its\n" ) ;
    Rprintf(  "    origin if you use it for a published work. Finally I would be\n" ) ;
    Rprintf(  "    interested to know for what kind of problem you have found\n" ) ;
    Rprintf(  "    this program to be of use.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Van Mo Dang\n" ) ;
    Rprintf(  "Van.Mo.Dang@utc.fr\n" ) ;
    Rprintf(  "http://www.hds.utc.fr/~mdang\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
} /* end of PrintHelpGeneral() */


void PrintHelpOptions( FILE* F )
{
    Rprintf(  "\n" ) ;
    Rprintf(  "Command Syntax\n" ) ;
    Rprintf(  "==============\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " Usage :    nem_exe   file  K  [ option1 option2 ... ]\n" ) ;
    Rprintf(  " ------- \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " Arguments :\n" ) ;
    Rprintf(  " -----------\n" ) ;
    Rprintf(  "   file       base name of input files ___.str and ___.dat\n" ) ;
    Rprintf(  "   K          number of classes\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " Options :   [ default  ] { possible values }\n" ) ;
    Rprintf(  " --------- \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " \n" ) ;
    Rprintf(  "  -a algo   [ nem      ]   { nem ncem gem }\n" ) ;
    Rprintf(  "     Algorithm to compute classification at E-step of each iteration : \n" ) ;
    Rprintf(  "      ncem = crisp classification by ICM procedure\n" ) ;
    Rprintf(  "      nem  = fuzzy classification by mean field approximation\n" ) ;
    Rprintf(  "      gem  = fuzzy classification by Gibbs sampling (Monte-Carlo simulations)\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -b beta   [ 1        ]   (0.0 <= beta <= 4)\n" ) ;
    Rprintf(  "     Coefficient of spatial smoothing to apply. This matches the \n" ) ;
    Rprintf(  "     Potts random field strength of interaction with 4-neighbor contexts.\n" ) ;
    Rprintf(  "     Notice that b = 0 is equivalent to EM for a mixture model.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -c wh thr [clas  0.04]   { none clas crit } and (>0)\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "       none     = no convergence test, i.e. do all specified iterations.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "       clas thr = stop the iterations when the largest difference\n" ) ;
    Rprintf(  "                  between previous and current classification matrix\n" ) ;
    Rprintf(  "                  is <= threshold. A threshold 0.04 is usually optimal.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "       crit thr = stop the iterations when \n" ) ;
    Rprintf(  "                  | (current_crit - last_crit)/current_crit | < threshold.\n" ) ;
    Rprintf(  "                  A threshold of 0.001 is usually best.  The test uses \n" ) ;
    Rprintf(  "                  the criterion selected by the option -C.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -f format [ hard     ]   { hard fuzzy }\n" ) ;
    Rprintf(  "     Format of output partition. Hard = N integers having values from\n" ) ;
    Rprintf(  "     1 to nk : for each observation, give the number of the class where\n" ) ;
    Rprintf(  "     it has highest grade of membership. Fuzzy = N x K reals between\n" ) ;
    Rprintf(  "     0 and 1 : for each observation, give its grade of membership in\n" ) ;
    Rprintf(  "     each class.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -i itmax  [ 100      ]   (>= 0)\n" ) ;
    Rprintf(  "     Maximum number of NEM iterations.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -l dolog  [ n        ]   { y n }\n" ) ;
    Rprintf(  "     Produce a log file or not to see the results of each iteration.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -m f p d  [norm p_ s__]  { norm lapl } { p_ pk } { s__ sk_ s_d skd }\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "     Mixture model assumption to use\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "      norm/lapl/bern : normal, Laplace or Bernoulli distributions\n" ) ;
    Rprintf(  "                       Bernoulli distributions are for binary data\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "      p_/pk     : clusters have equal / varying proportions\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "      s__/...   : variance model\n" ) ;
    Rprintf(  "                    s__ : same variance in all clusters and variables\n" ) ;
    Rprintf(  "                    sk_ : one variance per cluster, same in all variables\n" ) ;
    Rprintf(  "                    s_d : one variance per variable, same in all clusters\n" ) ;
    Rprintf(  "                    skd : one variance per cluster and variable\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "                  the variables are assumed independent within a cluster\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -n neigh  [ 4        ]   { 4 f }\n" ) ;
    Rprintf(  "     Neighborhood specification to use in the case of an image.\n" ) ;
    Rprintf(  "     4 = default 4-nearest neighbor system. f = specify neighborhood\n" ) ;
    Rprintf(  "     window in file.nei.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -o fout   [ file     ]\n" ) ;
    Rprintf(  "     Output files basename ___.cf, ___.mf and ___.log. Default is \n" ) ;
    Rprintf(  "     to use input file basename. Specify '-' to output the\n" ) ;
    Rprintf(  "     classification to standard output; useful to pipe the result to an \n" ) ;
    Rprintf(  "     'nem_exe -s f -' session.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -s init   [ s 1      ]   { s <v> | f <ini.uf> | r <n> | l <file> } \n" ) ;
    Rprintf(  "     Initialization mode.\n" ) ;
    Rprintf(  "     -s s <v>   \n" ) ;
    Rprintf(  "        Sort observations by variable <v>, then divide them\n" ) ;
    Rprintf(  "        in K quantiles of equal size to get initial partition.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "     -s f <ini.uf>\n" ) ;
    Rprintf(  "        Read initial fuzzy classification from file <ini.uf>. \n" ) ;
    Rprintf(  "        '-s f -' reads from standard input.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "     -s r <n>\n" ) ;
    Rprintf(  "        Start <n> times from random parameters (means chosen at\n" ) ;
    Rprintf(  "        random among the observations), then keep result with highest\n" ) ;
    Rprintf(  "        criterion.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "     -s l <file>\n" ) ;
    Rprintf(  "        Use partially known labels given in <file> to compute initial \n" ) ;
    Rprintf(  "        parameters. Those labels remain fixed throughout the\n" ) ;
    Rprintf(  "        clustering process.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "     -s mi <para> \n" ) ;
    Rprintf(  "     -s mf <para>\n" ) ;
    Rprintf(  "        Use specified parameters at beginning (mi) or throughout the\n" ) ;
    Rprintf(  "        clustering process (mf). Parameters syntax :\n" ) ;
    Rprintf(  "        p_1 ... p_{K-1}   m_11 .. m_1D m21 ... m_KD   s_11 ... s_KD.\n" ) ;
    Rprintf(  "        The s_kd are the standard errors for normal distributions, or\n" ) ;
    Rprintf(  "        the scale parameters for the Laplace distributions.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -t tie    [ random   ]   { random first }\n" ) ;
    Rprintf(  "     How to choose the class with highest probability when several\n" ) ;
    Rprintf(  "     classes have same maximum probability in MAP\n" ) ;
    Rprintf(  "     classification. 'random' draws uniformly between ex-aequo\n" ) ;
    Rprintf(  "     classes, 'first' chooses class with lowest index.\n" ) ;
    Rprintf(  " \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -B bmod   [ fix      ]   { fix psgrad heu_d heu_l }\n" ) ;
    Rprintf(  "     Procedure to estimate beta automatically :\n" ) ;
    Rprintf(  "      fix   = no estimation of beta, use beta given by option '-b'\n" ) ;
    Rprintf(  "      psgrad = pseudo-likelihood gradient ascent\n" ) ;
    Rprintf(  "      heu_d = heuristic using drop of fuzzy within cluster inertia\n" ) ;
    Rprintf(  "      heu_l = heuristic using drop of mixture likelihood\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -C crit   [ U        ]   local maximum criterion { U M D L }\n" ) ;
    Rprintf(  "     Criterion used to select the best local solution from random starts :\n" ) ;
    Rprintf(  "      U = fuzzy spatial clustering criterion\n" ) ;
    Rprintf(  "      M = fuzzy pseudo-likelihood\n" ) ;
    Rprintf(  "      D = fuzzy within cluster inertia\n" ) ;
    Rprintf(  "      L = likelihood of mixture parameters\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -G nit conv step rand [  1 0.001 0.0 0 ]  \n" ) ;
    Rprintf(  "     Parameters of beta gradient estimation\n" ) ;
    Rprintf(  "      nit = number of gradient iterations\n" ) ;
    Rprintf(  "      conv = threshold to test convergence (|g'|<conv*N is tested)\n" ) ;
    Rprintf(  "      step = > 0 for fixed step, <= 0 for Newton step = 1/g''\n" ) ;
    Rprintf(  "      rand = in -s r  init mode, initial beta random (1) or fixed by -b (0)\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -H bstep bmax ddrop dloss lloss [ 0.10 2.0 0.8 0.5 0.02 ]  \n" ) ;
    Rprintf(  "     Parameters of beta D and L heuristics :\n" ) ;
    Rprintf(  "      bstep = step of beta increase\n" ) ;
    Rprintf(  "      bmax  = maximal value of beta to test\n" ) ;
    Rprintf(  "      ddrop = threshold of allowed D drop (higher = less detection)\n" ) ;
    Rprintf(  "      dloss = threshold of allowed D loss (higher = less detection)\n" ) ;
    Rprintf(  "      lloss = threshold of allowed L loss (higher = less detection)\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -I eiter  [ 1        ]   (>= 1)\n" ) ;
    Rprintf(  "     Number of internal E-step iterations (nem and ncem algorithms),\n" ) ;
    Rprintf(  "     i.e. number of sweeps through whole dataset to compute\n" ) ;
    Rprintf(  "     classification at each iteration.  For gem algorithm, indicates\n" ) ;
    Rprintf(  "     number of sweeps through the dataset to compute the average\n" ) ;
    Rprintf(  "     frequency of class occurrence -> a large value is recommended\n" ) ;
    Rprintf(  "     for the gem algorithm (typically 50).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -M miss   [ replace  ]    { replace ignore }\n" ) ;
    Rprintf(  "     How to deal with missing data.  Replace by expected value (EM) or\n" ) ;
    Rprintf(  "     ignore when computing mean (maximize fuzzy clustering criterion).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -O order  [ direct   ]   order of site visit { direct random }\n" ) ;
    Rprintf(  "     Order in which to classify the observations at E-step.\n" ) ;
    Rprintf(  "      direct = given order 1..N\n" ) ;
    Rprintf(  "      random = random permutation of 1..N\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -S seed   [ <time>  ]    (integer)\n" ) ;
    Rprintf(  "     Specify a seed for the random number generator. Default uses\n" ) ;
    Rprintf(  "     current system clock.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -T test   [ n       ]    { y n }\n" ) ;
    Rprintf(  "     Print some debugging information or not.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -U update [ seq      ]   { seq para }\n" ) ;
    Rprintf(  "     Update the class of the sites in a sequential or parallel manner.\n" ) ;
    Rprintf(  "     'seq' works best. 'para' is more grounded, because it is EM with\n" ) ;
    Rprintf(  "     mean field approximation, but it requires a low spatial smoothing.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "You may also just type arguments : \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  -v                      versions information\n" ) ;
    Rprintf(  "  -h help_topic           longer help - help topics are\n" ) ;
    Rprintf(  "     general\n" ) ;
    Rprintf(  "     options\n" ) ;
    Rprintf(  "     examples\n" ) ;
    Rprintf(  "     filein\n" ) ;
    Rprintf(  "     fileout\n" ) ;
    Rprintf(  "     versions\n" ) ;
    Rprintf(  " \n" ) ;
    Rprintf(  " \n" ) ;
} /* end of PrintHelpOptions() */


void PrintHelpExamples( FILE* F )
{
    Rprintf(  "\n" ) ;
    Rprintf(  "    Examples :\n" ) ;
    Rprintf(  "    ----------\n" ) ;
    Rprintf(  " nem_exe  myimg  3 -b 0.5 -s r 10 -o Res/myimg3r_05 >&! Res/myimg3r_05.out &\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    This unix command clusters data set myimg (files myimg.dat, myimg.str)\n" ) ;
    Rprintf(  "    into 3 classes ; spatial coefficient is 0.5 ; do 10 random starts ;\n" ) ;
    Rprintf(  "    save results in files Res/myimg3r_05.* (___.cf ___.log ___.mf) .\n" ) ;
    Rprintf(  "    Last part of the command is unix-specific : it saves screen output\n" ) ;
    Rprintf(  "    to file Res/myimg3r_05.out and executes the program in background.\n" ) ;
    Rprintf(  "\n" ) ;
} /* end of PrintHelpExamples() */


void PrintHelpFileIn( FILE* F )
{
    Rprintf(  "\n" ) ;
    Rprintf(  "Input Files\n" ) ;
    Rprintf(  "===========\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    Input files are in ASCII format.\n" ) ;
    Rprintf(  "    2 or 3 input files are required : file.str   file.dat   [ file.nei ]\n" ) ;
    Rprintf(  "    2 optional input files :          file.u0    file.ck\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " 1) file.str\n" ) ;
    Rprintf(  " -----------\n" ) ;
    Rprintf(  "    This gives the structure of the data : \n" ) ;
    Rprintf(  "    type of spatial repartition (image, spatial or non-spatial),\n" ) ;
    Rprintf(  "    number of objects and variables. This file may start with\n" ) ;
    Rprintf(  "    comment lines to describe the dataset (lines beginning with #).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Ex 1 : Color image (each pixel is described by 3 variables) \n" ) ;
    Rprintf(  "             of 200 lines (height) and 300 columns (width)\n" ) ;
    Rprintf(  "    # RGB biomedical coloscopic image. Look for 3 or 4 classes.\n" ) ;
    Rprintf(  "    I 200 300 3\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Ex 2 : Spatial dataset of 3000 objects described by 4 variables\n" ) ;
    Rprintf(  "    # Economic data on counties (region of Centre). About 5 classes.\n" ) ;
    Rprintf(  "    S 3000 4\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Ex 3 : Non-spatial dataset of 3000 objects described by 4 variables\n" ) ;
    Rprintf(  "    # Companies characteristics for risk assessment. About 6 classes.\n" ) ;
    Rprintf(  "    N 3000 4\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " 2) file.dat\n" ) ;
    Rprintf(  " -----------\n" ) ;
    Rprintf(  "    Contains the objects-variables table (only the non-spatial\n" ) ;
    Rprintf(  "    variables).  Missing data are specified by NaN.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    If the dataset is an image, the pixels must be listed \n" ) ;
    Rprintf(  "    line-by-line first, i.e. : \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    x(1,1) x(1,2) ... x(1,nc)  \n" ) ;
    Rprintf(  "    x(2,1) x(2,2) ... x(2,nc)\n" ) ;
    Rprintf(  "    ...\n" ) ;
    Rprintf(  "    x(nl,1) x(nl,2) ... x(nl,nc)\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    where nl = number of lines, nc = number of columns, and \n" ) ;
    Rprintf(  "    x(i,j) = values of pixel at line i and column j (a set of 3\n" ) ;
    Rprintf(  "    numbers for a color image).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Ex 1 : Color image\n" ) ;
    Rprintf(  "    50 100 120\n" ) ;
    Rprintf(  "    51  99 122\n" ) ;
    Rprintf(  "    ...\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Ex 2 : Spatial data (4 variables)\n" ) ;
    Rprintf(  "    0.31 200 41 1200 \n" ) ;
    Rprintf(  "    0.28 202 43 1180\n" ) ;
    Rprintf(  "    ...\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " 3) file.nei\n" ) ;
    Rprintf(  " -----------\n" ) ;
    Rprintf(  "    Specifies the spatial relationships between the objects. \n" ) ;
    Rprintf(  "    * For an image, this file is optional and allows to specify\n" ) ;
    Rprintf(  "    a particular neighborhood system ; if no file is specified,\n" ) ;
    Rprintf(  "    default neighborhood system taken is 4-nearest neighbours.\n" ) ;
    Rprintf(  "    * For other spatial data, this file is required.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    Format is different for an image and other spatial data\n" ) ;
    Rprintf(  "    (see examples below).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Ex 1 : 4 nearest-neighbours in an image (default)\n" ) ;
    Rprintf(  "    -1 1           /* at most 1 pixel on the left and 1 on the right */\n" ) ;
    Rprintf(  "    -1 1           /* at most 1 pixel up and 1 down */\n" ) ;
    Rprintf(  "    0 1 0          \n" ) ;
    Rprintf(  "    1 0 1          /* 4 equally weighted neighbors : */\n" ) ;
    Rprintf(  "    0 1 0	   /* up,left,right,down */            \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Ex 2 : other spatial data, unweighted neighborhood graph\n" ) ;
    Rprintf(  "    0              /* 0 = no weight specified (use default weight = 1.0) */\n" ) ;
    Rprintf(  "    1   3   2 5 7  /* object 1 has 3 neighbors : objects 2, 5 and 7 */\n" ) ;
    Rprintf(  "    2   2   1 3\n" ) ;
    Rprintf(  "    ...\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Ex 2 : other spatial data, weighted neighborhood graph\n" ) ;
    Rprintf(  "    1                           /* 1 = weights are specified */\n" ) ;
    Rprintf(  "    1   3   2 5 7  0.5 0.6 0.8  /* object 1 has 3 neighbors : 2, 5 and 7 */\n" ) ;
    Rprintf(  "    2   2   1 3    0.3 0.5\n" ) ;
    Rprintf(  "    ...\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " 4.a) file.ck  (option -s l file.ck)\n" ) ;
    Rprintf(  " ------------\n" ) ;
    Rprintf(  "    Gives the class of the points for which the label is already\n" ) ;
    Rprintf(  "    known. Expected format is N + 1 integers (N being the total\n" ) ;
    Rprintf(  "    number of objects in the sample) :\n" ) ;
    Rprintf(  "    - the first number is the number of classes K ;\n" ) ;
    Rprintf(  "    - the N following integers, in the same order as the\n" ) ;
    Rprintf(  "    objects in file.dat, indicate the class to which the\n" ) ;
    Rprintf(  "    corresponding object belongs.  If the value is 0 or greater than the\n" ) ;
    Rprintf(  "    number of classes, then the object is considered to have no known\n" ) ;
    Rprintf(  "    label (its label will be computed by the algorithm). \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    In the current implementation, each class must contain at least\n" ) ;
    Rprintf(  "    one observation with known label. Those pre-labeled observations\n" ) ;
    Rprintf(  "    are used to initialize the centers of the clusters.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    5     /* number of classes */\n" ) ;
    Rprintf(  "    0     /* object 1 has no known label */\n" ) ;
    Rprintf(  "    5     /* object 2 belongs to class 5 */\n" ) ;
    Rprintf(  "    1     /* object 3 belongs to class 1 */\n" ) ;
    Rprintf(  "    ...\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " 4.b) file.u0  (option -s f)\n" ) ;
    Rprintf(  " ------------\n" ) ;
    Rprintf(  "    Gives an initial fuzzy classification to start the algorithm. \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    0.9 0.1  /* object 1 : initial membership = 0.9/0.1 in class 1/2 */\n" ) ;
    Rprintf(  "    0.3 0.7  /* object 2 : initial membership = 0.3/0.7 in class 1/2 */\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
} /* end of PrintHelpFileIn() */


void PrintHelpFileOut( FILE* F )
{
    Rprintf(  "\n" ) ;
    Rprintf(  "Output Files\n" ) ;
    Rprintf(  "============\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    Output files are in ASCII format.\n" ) ;
    Rprintf(  "    2 files are output :              file.cf    file.mf\n" ) ;
    Rprintf(  "    1 optional file is output :       file.log\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " 1) file.cf (or file.uf for fuzzy partition)\n" ) ;
    Rprintf(  " -------------------------------------------\n" ) ;
    Rprintf(  "    Gives the partition found by the algorithm.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Ex 1 : Image segmented in 2 'hard' classes (linewise order, as the \n" ) ;
    Rprintf(  "           input data file)\n" ) ;
    Rprintf(  "    1        /* Object 1 belongs to class 1 */\n" ) ;
    Rprintf(  "    1        /* Object 2 belongs to class 1 */\n" ) ;
    Rprintf(  "    2        /* Object 3 belongs to class 2 */\n" ) ;
    Rprintf(  "    1\n" ) ;
    Rprintf(  "    2\n" ) ;
    Rprintf(  "    2\n" ) ;
    Rprintf(  "    ...\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "    - Ex 2 : Data segmented in 3 fuzzy classes\n" ) ;
    Rprintf(  "    0.2  0.7  0.1    /* Object 1 : class 1 with probability 0.2 , etc. */\n" ) ;
    Rprintf(  "    0.8  0.05 0.015\n" ) ;
    Rprintf(  "    ...\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " 2) file.mf\n" ) ;
    Rprintf(  " ----------\n" ) ;
    Rprintf(  "    \n" ) ;
    Rprintf(  "    This file gives the values of criteria optimized by the algorithm,\n" ) ;
    Rprintf(  "    and the parameters of the mixture calculated by the\n" ) ;
    Rprintf(  "    algorithm (means, proportions, scale parameter). Example: \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ START EXAMPLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" ) ;
    Rprintf(  "Command line : \n" ) ;
    Rprintf(  " \n" ) ;
    Rprintf(  "  nem_exe intro 2 -b 1.00 -f fuzzy -R intro.cr -o intro_b10 -s r 10 -C M \n" ) ;
    Rprintf(  " \n" ) ;
    Rprintf(  "Criteria U=NEM, D=Hathaway, L=mixture, M=markov ps-like, error\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "  56.3807    -2671.73    -1788.76    -2793.36   0.335938\n" ) ;
    Rprintf(  " \n" ) ;
    Rprintf(  "Beta (fixed)\n" ) ;
    Rprintf(  "  1.0000\n" ) ;
    Rprintf(  "Mu (4), Pk, and sigma (4) of the 2 classes\n" ) ;
    Rprintf(  " \n" ) ;
    Rprintf(  "   0.947  0.00456  0.0199  0.998   0.5    1.15226  1.15226  1.15226  1.15226 \n" ) ;
    Rprintf(  "   0.0614 1.03     0.864  -0.163   0.5    1.15226  1.15226  1.15226  1.15226 \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END EXAMPLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  " 3) file.log\n" ) ;
    Rprintf(  " -----------\n" ) ;
    Rprintf(  "    Details each iteration results (optimized criterion and \n" ) ;
    Rprintf(  "    class parameters).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
} /* end of PrintHelpFileOut() */


void PrintHelpVersions( FILE* F )
{
    Rprintf(  "\n" ) ;
    Rprintf(  "History of modifications\n" ) ;
    Rprintf(  "========================\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Version 0.00  (01.02.1996)\n" ) ;
    Rprintf(  "------------\n" ) ;
    Rprintf(  "First version released on WEB. Initialization by sorting variable.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Version 1.00  (31.05.1996)  \n" ) ;
    Rprintf(  "------------\n" ) ;
    Rprintf(  "Random initializations. Image default neighborhood system. Long help.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Version 1.01  (27.06.1996)  \n" ) ;
    Rprintf(  "------------\n" ) ;
    Rprintf(  "Added an -a ncem option to implement the crisp clustering version of NEM.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Each observation is now updated sequentially in turn. The previous \n" ) ;
    Rprintf(  "parallel updating would produce chessboard-like images for high betas.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Version 1.02  (17.10.1996)\n" ) ;
    Rprintf(  "------------\n" ) ;
    Rprintf(  "Added the possibility to take into account a partial knowledge\n" ) ;
    Rprintf(  "of the classification into the clustering procedure \n" ) ;
    Rprintf(  "(option '-s l <file.ck>').\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Version 1.03  (02.10.1997)  \n" ) ;
    Rprintf(  "------------\n" ) ;
    Rprintf(  "If a partial knowledge of the classification is available, the\n" ) ;
    Rprintf(  "intial cluster centers are computed from the observations with\n" ) ;
    Rprintf(  "known labels.  \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "The log file is made optional (option '-l y').  The ___.mf file now\n" ) ;
    Rprintf(  "also contains the estimated cluster proportions, volumes and\n" ) ;
    Rprintf(  "covariance matrices and the values of the final criteria.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Version 1.04  (11.01.1998)  \n" ) ;
    Rprintf(  "------------\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Two heuristics are implemented to estimate the spatial smoothness\n" ) ;
    Rprintf(  "coefficient beta.  The first heuristic is based on detecting a sharp\n" ) ;
    Rprintf(  "drop in the fuzzy within-cluster inertia D, or a sufficient decrease\n" ) ;
    Rprintf(  "from its maximum value, when beta is slowly increased.  The second\n" ) ;
    Rprintf(  "heuristic is based on detecting a sufficient decrease of the\n" ) ;
    Rprintf(  "log-likelihood L of the mixture parameters from its maximum value.\n" ) ;
    Rprintf(  "The heuristics may be invoked with '-B heu_d' or '-B heu_l'. Their\n" ) ;
    Rprintf(  "default parameters may be changed with '-H ...'. \n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "The final partition may now be printed to standard ouput instead of to a\n" ) ;
    Rprintf(  "file (option '-o -').  The result can thus be redirected as an \n" ) ;
    Rprintf(  "initial partition to another nem_exe session's input.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "The initial partition may now be read from any file with option '-s f\n" ) ;
    Rprintf(  "foo.bar' (previously inputbasename.u0 was used).  In particular, the\n" ) ;
    Rprintf(  "partition may be read from standard input (option '-s f -').  This allows\n" ) ;
    Rprintf(  "to read the initial partition through a pipe from the result of a previous \n" ) ;
    Rprintf(  "nem_exe session.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "At each iteration, the fuzzy classification at the E-step may now be\n" ) ;
    Rprintf(  "computed by either of two methods :\n" ) ;
    Rprintf(  "- Neighborhood EM's fixed point technique (default, '-a nem')\n" ) ;
    Rprintf(  "- Gibbsian EM's Gibbs sampler technique ('-a gem').\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "A longer explanation is given for the options and successive versions\n" ) ;
    Rprintf(  "(option '-h helptopic').\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Other options have been added, to test the effect of alternative\n" ) ;
    Rprintf(  "parameter values.  Those options usually do not change considerably\n" ) ;
    Rprintf(  "the default behaviour of the algorithm :\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "- At the E-step of each iteration, the classification may now be\n" ) ;
    Rprintf(  "  updated in a random order instead of 1..N (option '-O random').\n" ) ;
    Rprintf(  "- The seed of the random number generator may now be given (option '-S seed').\n" ) ;
    Rprintf(  "- The number of E-step internal iterations may be changed (option '-I eiter').\n" ) ;
    Rprintf(  "- The convergence threshold may also be changed (option '-c cvthres').\n" ) ;
    Rprintf(  "- Another criterion than U may be chosen to select best result ('-C crit').\n" ) ;
    Rprintf(  "- Compute the classification error in two-class case ('-R refclass').\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "A few internal changes were also made to make the program portable to\n" ) ;
    Rprintf(  "the djgpp gcc compiler for MS-DOS (srand48 replaced by srandom).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Version 1.05  (09-APR-1998)  \n" ) ;
    Rprintf(  "------------\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "This version mainly adds the capability to deal with missing data.\n" ) ;
    Rprintf(  "This means that some of the N observation vectors may be incompletely\n" ) ;
    Rprintf(  "observed.  In the 'name.dat' file, use NaN to indicate unobserved\n" ) ;
    Rprintf(  "components of an observation vector.  Two slightly different\n" ) ;
    Rprintf(  "techniques are provided to deal with missing data.  The first and\n" ) ;
    Rprintf(  "default behaviour (invoked using switch '-M replace') roughly consists\n" ) ;
    Rprintf(  "in replacing any missing component with its expected value.  This\n" ) ;
    Rprintf(  "technique implements the EM procedure and finds parameters maximizing\n" ) ;
    Rprintf(  "the likelihood.  The alternative behaviour (invoked using switch '-M\n" ) ;
    Rprintf(  "ignore') consists in ignoring missing components.  This means that the\n" ) ;
    Rprintf(  "means are computed using only observed components.  This alternative\n" ) ;
    Rprintf(  "technique finds a classification matrix and parameters which maximize\n" ) ;
    Rprintf(  "the fuzzy classifying log-likelihood.  It appears to converge a bit\n" ) ;
    Rprintf(  "faster than the 'replace' mode.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "The iteration count is displayed in a more economic way now (all the\n" ) ;
    Rprintf(  "iteration numbers were displayed separately).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "A new criterion is computed, the fuzzy pseudo-likelihood, named M.\n" ) ;
    Rprintf(  "Using this criterion in order to choose the best result \n" ) ;
    Rprintf(  "may prove less sensitive to the value of beta than using U.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "In the random start strategy, two initialization tactics have been\n" ) ;
    Rprintf(  "made more sensible.  The initial volumes are computed as whole volume\n" ) ;
    Rprintf(  "/ number of classes (the whole volume was used previously).  The means\n" ) ;
    Rprintf(  "are redrawn until all drawn means are different --- this avoids the\n" ) ;
    Rprintf(  "problem of artificially merging together two classes.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "A few internal changes were also made in order to allow direct call to\n" ) ;
    Rprintf(  "the program from as a Matlab function.  This allowed to detect and\n" ) ;
    Rprintf(  "remove a few potential bugs that had gone unnoticed (a file not\n" ) ;
    Rprintf(  "closed, use of memory just after freeing it).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Version 1.06  (26-FEB-1999)\n" ) ;
    Rprintf(  "------------\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "The beta parameter can now be estimated in a statistically more\n" ) ;
    Rprintf(  "justified way, by maximizing the fuzzy pseudo-likelihood at each\n" ) ;
    Rprintf(  "iteration.  This computation is actually equivalent to the M-step of\n" ) ;
    Rprintf(  "Zhang's (1992) mean field EM algorithm (our algorithm differs from\n" ) ;
    Rprintf(  "Zhang's through the E-step, by allowing a sequential update of the\n" ) ;
    Rprintf(  "classification instead of Zhang's parallel fixed point mechanism).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "The Laplace family of models can now be used instead of the normal\n" ) ;
    Rprintf(  "(Gaussian) family.  This may make the algorithm more robust to\n" ) ;
    Rprintf(  "outliers.  Indeed, the center of a Laplace distribution is estimated\n" ) ;
    Rprintf(  "as the median of the class, instead of the mean for normal\n" ) ;
    Rprintf(  "distributions.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Other changes in functionality:\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "- the parameters of the mixture of distributions (proportions,\n" ) ;
    Rprintf(  "  centers, dispersions) can be given by the user on the command line.\n" ) ;
    Rprintf(  "  Those given parameters can be used, whether as fixed values to be\n" ) ;
    Rprintf(  "  used throughout the algorithm ('-s mf <para>'), or as \n" ) ;
    Rprintf(  "  initial starting parameters ('-s mi <para>') ;\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "- by option '-U' the sites can be updated whether:\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "   - in a parallel way, i.e. the new classification matrix is computed\n" ) ;
    Rprintf(  "     as a whole, based on the whole old classification matrix ;\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "   - in a sequential way, i.e. each (eventually fuzzy) classification \n" ) ;
    Rprintf(  "   vector is updated based on the other classification vectors, and \n" ) ;
    Rprintf(  "   this updated classification vector is used to update the following\n" ) ;
    Rprintf(  "   classification vectors.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "- option '-c' is modified in order to allow different convergence tests : \n" ) ;
    Rprintf(  "   - 'none' : no convergence test, do the specified number iterations\n" ) ;
    Rprintf(  "   - 'clas <maxdif>' (default) : stop if the maximum change in the \n" ) ;
    Rprintf(  "     classification matrix over the last two iterations is less than '<maxdif>'\n" ) ;
    Rprintf(  "   - 'crit <dif>' : stop if the difference in the selected criterion\n" ) ;
    Rprintf(  "   over the last two iterations is less than '<dif>'\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "- option '-C' allows to select different criteria to choose between\n" ) ;
    Rprintf(  "  different local optima of the algorithm\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "A few internal changes were also made in order to accomodate for the\n" ) ;
    Rprintf(  "new functionalities.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Version 1.07  (08-APR-1999)\n" ) ;
    Rprintf(  "------------\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "This version adds the Bernoulli family of distributions which is\n" ) ;
    Rprintf(  "suited to binary data (or qualitative data encoded as binary).\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Reference :\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Dang, M. V. and Govaert, G. (1998). Fuzzy clustering of spatial binary\n" ) ;
    Rprintf(  "data.  Kybernetika, 34(4), pp.~393--398.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "Version 1.08  (15-JUL-1999)\n" ) ;
    Rprintf(  "------------\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "- Internal changes aim to reduce the computation time.\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "- Also internally, changed srandom/random to srand/rand to enhance\n" ) ;
    Rprintf(  "  the portability of the program\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
    Rprintf(  "\n" ) ;
} /* end of PrintHelpVersions() */


