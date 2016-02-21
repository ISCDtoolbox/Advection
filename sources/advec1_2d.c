#include "advect.h"


/* barycentric coordinates of point c[] in iel=p0,p1,p2 */
static inline int bar_2d(pPoint pp[3],double *c,int iel,double *cb) {
  double    det;
  char      i,i1,i2;

  det = (pp[1]->c[0] - pp[0]->c[0]) * (pp[2]->c[1] - pp[0]->c[1])
      - (pp[1]->c[1] - pp[0]->c[1]) * (pp[2]->c[0] - pp[0]->c[0]);
  if ( fabs(det) < AD_EPSD )  return(0);
  det = 1.0 / det;

  /* barycentric coordinate, by use of Cramer's formulas */
  for (i=0; i<2; i++) {
    i1 = (i+1) % 3;
    i2 = (i+2) % 3;
    cb[i] = (pp[i1]->c[0]-c[0])*(pp[i2]->c[1]-c[1]) - (pp[i1]->c[1]-c[1])*(pp[i2]->c[0]-c[0]);
    cb[i] *= det;
  }
  cb[2] = 1.0 - cb[0] - cb[1];

  return(1);
}


/* return velocity P1_interpolation in v in element iv[3] */
static inline double vecint_2d(double *un,int *iv,double *cb,double *v) {
  double   *u0,*u1,*u2,dd;

  u0 = &un[2*(iv[0]-1)+1];
  u1 = &un[2*(iv[1]-1)+1];
  u2 = &un[2*(iv[2]-1)+1];

  /* P1 interpolate of the speed */   
  v[0] = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0];  
  v[1] = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1];
  dd   = sqrt(v[0]*v[0] + v[1]*v[1]);

  return(dd);
}


/* find element containing c, starting from nsdep, return baryc. coord */
static int locelt_2d(pMesh mesh,int nsd,double *c,double *cb) {
  pTria     pt;
  pPoint    p0,p1,p2;
  double    ax,ay,bx,by,cx,cy,eps;
  double    aire1,aire2,aire3,dd; 
  int       i,sgn,nsf,nsp;

  nsf  = nsd;
  nsp  = 0;
  mesh->mark = ++mesh->mark;
  while ( nsf > 0 ) {
    pt = &mesh->tria[nsf];
    if ( pt->mark == mesh->mark )  return(nsp);
    pt->mark = mesh->mark;

    /* area of triangle */
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    bx = p2->c[0] - p0->c[0];
    by = p2->c[1] - p0->c[1];
    dd = ax*by - ay*bx;
    sgn = dd > 0.0 ? 1 : -1;
    eps = sgn == 1 ? -AD_EPS*dd : AD_EPS*dd;

    /* barycentric */
    bx = p1->c[0] - c[0];
    by = p1->c[1] - c[1];
    cx = p2->c[0] - c[0];
    cy = p2->c[1] - c[1];

    /* p in half-plane lambda_0 > 0 */
    aire1 = sgn*(bx*cy - by*cx);
    if ( aire1 < eps ) {
      nsp = nsf;
      nsf = pt->adj[0] / 3;
      if ( !nsf ) {
        cb[0] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    ax = p0->c[0] - c[0];
    ay = p0->c[1] - c[1];
    aire2 = sgn*(cx*ay - cy*ax);
    if ( aire2 < eps ) {
      nsp = nsf;
      nsf = pt->adj[1] / 3;
      if ( !nsf ) {
        cb[1] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    aire3 = sgn*dd - aire1 - aire2;
    if ( aire3 < eps ) {
      nsp = nsf;
      nsf = pt->adj[2] / 3;
      if ( !nsf ) {
        cb[2] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    aire1 = AD_MAX(aire1,0.0);
    aire2 = AD_MAX(aire2,0.0);
    aire3 = AD_MAX(aire3,0.0);
    dd    = aire1 + aire2 + aire3;
    if ( dd > AD_EPSD ) {
      dd = 1.0 / dd;
      cb[0] = aire1 * dd;
      cb[1] = aire2 * dd;
      cb[2] = aire3 * dd;
    }
    return(nsf);
  }

  /* no need for exhaustive search */
  return(nsp);
}


/* computes the characteristic line emerging from point with barycentric coordinates cb
 in triangle iel and follows the characteristic on a length dt at most. Most often has to 
 cross the boundary of the current triangle, thus stores in it the new triangle, in cb
 the barycentric coordinates in the new triangle of the crossing point, and updates dt 
 with the remaining time to follow characteristic line */
static int travel_2d(ADst *adst,double *cb,int *iel,double *dt) {
  pTria       pt;
  pPoint      p[3];
  double     *u0,*u1,*u2,m[3],ddt,tol,ux,uy,c[2],cb1[3];
  int         k,*adj;
  char        i,i0,i1,i2;

  tol = *dt;
  k   = *iel;
  pt  = &adst->mesh.tria[k];

  p[0] = &adst->mesh.point[pt->v[0]];
  p[1] = &adst->mesh.point[pt->v[1]];
  p[2] = &adst->mesh.point[pt->v[2]];

  /* the velocity at each vertex of iel */
  u0 = &adst->sol.u[2*(pt->v[0]-1)+1];
  u1 = &adst->sol.u[2*(pt->v[1]-1)+1];
  u2 = &adst->sol.u[2*(pt->v[2]-1)+1];

  /* u = P1 velocity at the point of barycentric coordinates cb */
  ux = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0];
  uy = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1];
  if ( ux*ux+uy*uy < AD_EPSD )  return(0);

  /* endpoint of the characteristic line starting from cb */
  /* segment is of norm sqrt(ux*ux+uy*uy), and is to be followed for time dt, in the - direction */
  c[0] = (cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0]) - tol*ux;   
  c[1] = (cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1]) - tol*uy;

  /* barycentric coordinate of uu in the current triangle */
  if ( !bar_2d(p,c,k,cb1) )  return(0);

  /* check if endpoint is in k */
  for (i=0; i<3; i++)
    if ( cb1[i] < 0.0 )  break;
  if ( i == 3 ) {
    memcpy(cb,cb1,3*sizeof(double));
    *dt -= tol;
    return(*dt > 0.0);
  }

  /* argument of the smallest value among 3 */
  ddt = AD_TGV;
  i0  = -1;
  for (i=0; i<3; i++) {
    m[i] = cb[i] - cb1[i];
    if ( m[i] > 0.0 ) {
      if ( tol*cb[i]/m[i] < ddt ) {
        ddt = tol*cb[i]/m[i];
        i0  = i;
      } 
    }
  }

  /* case when advection stays in same triangle */
  if ( ddt > tol ) {
    memcpy(cb,cb1,3*sizeof(double));
    *dt -= tol;
  }
  /* advection goes out of triangle: advect a minimum value */
  if ( ddt < AD_EPS*tol ) {
    c[0] = cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0] - AD_EPS*ux;
    c[1] = cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1] - AD_EPS*uy;
    /* find the new triangle */
    k = locelt_2d(&adst->mesh,k,c,cb); 
    if ( !k )  return(0);
    *iel = k;
    *dt -= AD_EPS;
  }
  else {
    /* barycentric coordinates of the exit point */   
    for (i=0; i<3; i++)
      cb1[i] = cb[i] - ddt * m[i] / tol;
    *dt -= ddt;

    /* find output triangle */
    if ( !pt->adj[i0] ) {
      memcpy(cb,cb1,3*sizeof(double));
      return(0);
    }
    else {
      i1 = (i0+1) % 3;
      i2 = (i0+2) % 3;
      *iel = pt->adj[i0] / 3;
      i0   = pt->adj[i0] % 3; 
      cb[i0]         = 0.0;
      cb[(i0+1) % 3] = cb1[i2];
      cb[(i0+2) % 3] = cb1[i1];
    }
  }

  return(*dt > 0.0);
}


/* 4th order Runge-Kutta: for backtracking characteristic line, step is <0 */
/* v = initial speed at point c */
static int nxtptR_2d(ADst *adst,int *iel,double *c,double *cb,double step,double *v) {
  double  h6,v1[2],v2[2],v3[2];
  double  xp1[2],xp2[2],xp3[2],cc[2];
  int     k;

  /* first integration point = middle point*/
  k = *iel;
  xp1[0] = c[0] - 0.5*step*v[0];
  xp1[1] = c[1] - 0.5*step*v[1];
  k = locelt_2d(&adst->mesh,k,xp1,cb);
  if ( k < 1 )   return(k);
  vecint_2d(adst->sol.u,adst->mesh.tria[k].v,cb,v1);  

  xp2[0] = c[0] - 0.5*step*v1[0];
  xp2[1] = c[1] - 0.5*step*v1[1];
  k = locelt_2d(&adst->mesh,k,xp2,cb);
  if ( k < 1 )  return(k);
  vecint_2d(adst->sol.u,adst->mesh.tria[k].v,cb,v2);

  xp3[0] = c[0] - step*v2[0];
  xp3[1] = c[1] - step*v2[1];
  k = locelt_2d(&adst->mesh,k,xp3,cb);
  if ( k < 1 )   return(k);
  vecint_2d(adst->sol.u,adst->mesh.tria[k].v,cb,v3);

  /* final RK4 step */  
  h6 = step / 6.0;
  cc[0] = c[0] - h6 * (v[0] + 2*(v1[0] + v2[0]) + v3[0]);
  cc[1] = c[1] - h6 * (v[1] + 2*(v1[1] + v2[1]) + v3[1]);
  k = locelt_2d(&adst->mesh,k,cc,cb);
  if ( k < 1 )   return(k);

  /* update */  
  c[0] = cc[0];
  c[1] = cc[1];
  vecint_2d(adst->sol.u,adst->mesh.tria[k].v,cb,v);

  *iel = k;
  return(1);
}


/* Euler scheme for backtracking characteristic line */
static int nxtptE_2d(ADst *adst,int *iel,double *c,double *cb,double step,double *v) {
  double  cc[2];
  int     k;

  k = *iel;
  cc[0] = c[0] - step*v[0];
  cc[1] = c[1] - step*v[1];
  k = locelt_2d(&adst->mesh,k,cc,cb);
  if ( k < 1 )   return(k);

  /* update */  
  c[0] = cc[0];
  c[1] = cc[1];
  vecint_2d(adst->sol.u,adst->mesh.tria[k].v,cb,v);
  *iel = k;

  return(1);
}


/* solve advection, solution in rv */
int advec1_2d(ADst *adst) {
  pTria    pt,pt1;
  pPoint   ppt;
  double   cb[3],v[2],c[2],dt,dte,norm,tol,step;  
  int      j,k,ip,iel,kprv,nt,nstep;
  char     i;

  if ( !adst->sol.new ) {
    adst->sol.new = (double *)malloc((adst->info.np+1)*sizeof(double));
    assert(adst->sol.new);
    memcpy(adst->sol.new,adst->sol.chi,(adst->info.np+1)*sizeof(double));
  }

  /* check mesh size and velocity */
  if ( adst->sol.dt < 0.0 ) {
    adst->sol.dt = adst->sol.hmin / adst->sol.umax;
  }

  dt    = adst->sol.dt;
  tol   = adst->sol.dt / 2;
  nstep = (int)(dt/tol);
  step  = dt / nstep;

  if ( adst->info.verb != '0' ) {
    fprintf(stdout,"    Time stepping: %g\n",dt);
    fprintf(stdout,"    Solving: "); fflush(stdout);
  }

  ++adst->mesh.mark;
  nt = 0;
  for (k=1; k<=adst->info.nt; k++) {
    pt = &adst->mesh.tria[k];
    for (i=0; i<3; i++) {
      ip  = pt->v[i];
      ppt = &adst->mesh.point[ip];
      /* check if already processed */
      if ( ppt->flag == 1 )  continue;
      ppt->flag = 1;
      kprv = k;

      /* coordinate and velocity at starting point */
      v[0] = adst->sol.u[2*(ip-1)+1];
      v[1] = adst->sol.u[2*(ip-1)+2];
      norm = sqrt(v[0]*v[0] + v[1]*v[1]);
      if ( norm < AD_EPSD )  continue;

      /* barycentric coordinates of point p in triangle k */
      memset(cb,0,3*sizeof(double));
			cb[i] = 1.0;

      /* next point = foot of the characteristic line */
      c[0] = ppt->c[0];
      c[1] = ppt->c[1];
      dte  = dt;
      for (iel=k,j=0; j<nstep; j++) {
        kprv = iel;
        if ( nxtptR_2d(adst,&iel,c,cb,step,v) < 1 )  break;
        dte -= step;
      }
      if ( j < nstep ) {
        iel = kprv;
        while ( travel_2d(adst,cb,&iel,&dte) );
      }
      /* check if characteristic remains inside domain */
      if ( dte > AD_EPS ) {
        if ( iel < 0 )  iel = kprv;
        iel = locelt_2d(&adst->mesh,iel,c,cb);
				if ( iel < 1 )  iel = kprv;
      }
      /* interpolate value at foot  */
      if ( iel == 0 )  return(0);
      pt1 = &adst->mesh.tria[iel];
      adst->sol.new[ip] = cb[0]*adst->sol.chi[pt1->v[0]] \
                        + cb[1]*adst->sol.chi[pt1->v[1]] \
                        + cb[2]*adst->sol.chi[pt1->v[2]];
      nt++;
    }
  }

  if ( adst->info.verb != '0' )
    fprintf(stdout,"%d characteristics\n",nt);

  return(1);
}