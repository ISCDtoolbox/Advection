#include "advect.h"
#include "ad_calls.h"


int AD_advect(ADst *adst) {
  int    ier;

  if ( adst->info.dim == 2)
		ier = advec1_2d(adst);
	else
		ier = advec1_3d(adst);

  return(ier);
}