#ifndef _ADVECT_H
#define _ADVECT_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "ad_calls.h"

#define AD_VER   "3.0a"
#define AD_REL   "Feb. 19, 2016"
#define AD_CPY   "(C) Copyright 2007- , ICS-SU"

#define AD_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define AD_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define AD_MIN3(a,b,c) ( (a) < (b) ? ((a)<(c) ? (a) : (c)) : ((b)<(c) ? (b) : (c)) )
#define AD_MAX3(a,b,c) ( (a) > (b) ? ((a)>(c) ? (a) : (c)) : ((b)>(c) ? (b) : (c)) )

#define AD_EPS     1.e-6
#define AD_EPS1    1.e-4
#define AD_EPS2    1.e-12
#define AD_EPSD    1.e-200
#define AD_TGV     1.e+30
#define AD_DTM     1.e-4
#define AD_LONMAX  1024


/* data structures */
typedef struct {
  double    c[3];
  int       s,ref;
  char      flag;
} Point;
typedef Point * pPoint;

typedef struct {
  int       v[3],adj[3],ref,mark;
} Tria;
typedef Tria * pTria;

typedef struct {
  int       v[4],adj[4],ref,mark;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
  int       dim,ver;
  int       np,nt,ne,mark;
  char      verb,nocfl,noex;
  mytime    ctim[TIMEMAX];
} Info;

typedef struct {
  int       mark;
  char     *name;
  Point    *point;
  Tria     *tria;
  Tetra    *tetra;
} Mesh;
typedef Mesh * pMesh;

typedef struct {
  double  *u,*chi,*new,dt,hmin,umax;
  char    *namein,*nameout,*namechi;
  /* ri = initialized rho, rv = rho at each vertex */
} Sol;
typedef Sol * pSol;

struct _ADst {
  Mesh    mesh;
  Sol     sol;
  Info    info;
};


/* prototypes */
int   loadMesh(ADst *adst);
int   loadSol(ADst *adst);
int   loadChi(ADst *adst);
int   saveChi(ADst *adst);
int   hashel_2d(ADst *adst);
int   hashel_3d(ADst *adst);
int   advec1_2d(ADst *adst);
int   advec1_3d(ADst *adst);

int   AD_advect(ADst *adst);


#endif
