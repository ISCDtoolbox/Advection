#include "advect.h"

/* Tangential projection of vector field v onto the plane of triangle k */
static int tanproj(ADst *adst,int k,double *v) {
  pTria     pt;
  pPoint    p0,p1,p2;
  double    n[3],det,ps,norm;
  
  pt = &adst->mesh.tria[k];
  p0 = &adst->mesh.point[pt->v[0]];
  p1 = &adst->mesh.point[pt->v[1]];
  p2 = &adst->mesh.point[pt->v[2]];
  
  n[0] = (p1->c[1]-p0->c[1])*(p2->c[2]-p0->c[2]) - (p1->c[2]-p0->c[2])*(p2->c[1]-p0->c[1]);
  n[1] = (p1->c[2]-p0->c[2])*(p2->c[0]-p0->c[0]) - (p1->c[0]-p0->c[0])*(p2->c[2]-p0->c[2]);
  n[2] = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);
  
  norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if ( norm > AD_EPS2 ) norm = 1.0 / norm;
  
  n[0] *= norm;
  n[1] *= norm;
  n[2] *= norm;
  ps = v[0]*n[0] + v[1]*n[1] + v[2]*n[2];
  
  v[0] -= ps*n[0];
  v[1] -= ps*n[1];
  v[2] -= ps*n[2];

  return(1);
}

/* Calculate the (triangle) ball of point i in triangle k */
static int boulet_s(pMesh mesh,int start,char i,int *list) {
  pTria       pt;
  int         ilist,k;
  char        i0,i1,i2,voy;
  
  ilist = 0;
  k = start;
  i0 = i;
  
  do {
    pt = &mesh->tria[k];
    i1 = (i0+1) % 3;
    
    list[ilist] = 3*k+i0;
    ilist++;
    if ( ilist > AD_LONMAX-2 ) return(-ilist);
    
    k = pt->adj[i1] / 3;
    voy = pt->adj[i1] % 3;
    i0 = (voy+1) % 3;
  }
  while ( k && k != start );
  
  if ( k == start )  return(ilist);
  
  /* Reverse loop */
  k = start;
  pt = &mesh->tria[k];
  i0 = i;
  i2 = (i0+2) % 3;
  k  = pt->adj[i2] / 3;
  voy = pt->adj[i2] % 3;
  i0 = (voy+2) % 3;

  while ( k ) {
    list[ilist] = 3*k+i0;
    ilist++;
    if ( ilist > AD_LONMAX-2 ) return(-ilist);
    
    pt = &mesh->tria[k];
    i2 = (i0+2) % 3;
    k = pt->adj[i2] / 3;
    voy = pt->adj[i2] % 3;
    i0 = (voy+2) % 3;
  }
  
  return(ilist);
}

/* */
static int travel_s(ADst *adst,double *cb,int *iel,double *dt) {
  pTria     pt;
  pPoint    p[3];
  double    *u0,*u1,*u2,u[3];
  int       nz,k;
  char      i,i0,i1,i2,ia;
  
  nz = 0;
  for (i=0; i<3; i++)
    if ( cb[i] < AD_EPS ) {
      nz++;
      i0 = i;
    }
  else
    ia = i;
  
  /* Advection starts from vertex i0 of pt */
  if ( nz == 2 ) {
    
    // take boulet_s
    
  }
  /* Advection starts from inside edge ia */
  if ( nz == 1 ) {
    // Check direction
  }
  /* Advection starts from inside pt */
  k = *iel;
  pt = &adst->mesh.tria[k];
  p[0] = &adst->mesh.point[pt->v[0]];
  p[1] = &adst->mesh.point[pt->v[1]];
  p[2] = &adst->mesh.point[pt->v[2]];
  
  /* u = P1 velocity at the point with barycentric coordinates cb */
  u[0] = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0];
  u[1] = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1];
  u[2] = cb[0]*u0[2] + cb[1]*u1[2] + cb[2]*u2[2];
  tanproj(adst,k,u);
  
  if ( u[0]*u[0]+u[1]*u[1]+u[2]*u[2] < AD_EPSD )  return(0);
  
  return(1);
}

/* Find next point in characteristic backtracking by the 4th order Runge-Kutta method */
static int nxtptR_s(ADst *adst,int *iel,double *c,double *cb,double step,double *v) {
  return(1);
}

/* Solve advection */
int advec1_s(ADst *adst) {
  pTria      pt;
  double     dt,tol;
  int        nstep,j,k,kprv,ip;
  char       i;

  /* allocate structure for solution */
  if ( !adst->sol.new ) {
    adst->sol.new = (double *)malloc((adst->info.np+1)*sizeof(double));
    assert(adst->sol.new);
    memcpy(adst->sol.new,adst->sol.chi,(adst->info.np+1)*sizeof(double));
  }
  
  if ( !adst->info.nocfl ) {
    printf(" *** Unavailable option -- Abort.\n");
    exit(0);
  }
  
  dt    = adst->sol.dt;
  nstep = 100;
  tol   = adst->sol.dt / nstep;
  
  if ( adst->info.verb != '0' ) {
    fprintf(stdout,"    Time stepping: %g\n",dt);
    fprintf(stdout,"    Solving: "); fflush(stdout);
  }
  
  /* Backtrack characteristic line emerging from all points */
  /* for (k=1; k<=adst->info.nt; k++) {
    pt = &adst->mesh.tria[k];
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      ppt = &adst->mesh.point[ip];
      if ( ppt->flag == 1 ) continue;
      kprv = k;
      
    }
  }*/
  
  exit(0);

  return(1);
}
