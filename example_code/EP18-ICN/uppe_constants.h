#ifndef __uppe_constants_h__
#define __uppe_constants_h__

/* some fundamental constants */

#define cnst_c         ( 299792458.0        )
#define cnst_mu0       ( M_PI*4.0e-07       )
#define cnst_epsilon0  ( 1.0/(cnst_c*cnst_c*cnst_mu0) )
#define cnst_ht        ( 1.05457266e-34     )
#define cnst_e         ( 1.60217733e-19     )

#define cnst_me        ( 9.10938970e-31     )

/* Atomic units               */

#define cnst_au_Eh     ( 4.35974417e-18     )
#define cnst_au_a0     ( 0.5291772108e-10   )
#define cnst_au_time   ( 2.418884326505e-17 )
#define cnst_au_EF     ( 5.14220642e+11     )



/* some unit conversions      */

#define om2eV          ( cnst_ht / cnst_e  )
#define eV2om          ( cnst_e  / cnst_ht )
#define om2nm          ( 2.0*M_PI*cnst_c*1.0e9   )
#define nm2om          ( 2.0*M_PI*cnst_c*1.0e9   )
#define eV2nm          ( 2.0*M_PI*cnst_c*1.0e+9/eV2om )

#endif // __uppe_constants_h__
