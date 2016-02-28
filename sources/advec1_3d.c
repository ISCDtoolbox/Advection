#include "advect.h"


static double det_3d(double *cp0,double *cp1,double *cp2,double *cp3) {
  double  det,x01,x02,x03,y01,y02,y03,z01,z02,z03; 

  x01 = cp1[0] - cp0[0];  y01 = cp1[1] - cp0[1];  z01 = cp1[2] - cp0[2];
  x02 = cp2[0] - cp0[0];  y02 = cp2[1] - cp0[1];  z02 = cp2[2] - cp0[2];
  x03 = cp3[0] - cp0[0];  y03 = cp3[1] - cp0[1];  z03 = cp3[2] - cp0[2];

  det = x01*(y02*z03-y03*z02) + y01*(z02*x03-z03*x02) + z01*(x02*y03-x03*y02);
  return(det);  
}

/* barycentric coordinates of point c[] in iel=p0,p1,p2,p3 */
static inline int bar_3d(pPoint pp[4],double *c,int iel,double *cb) { 
  double    det,cp0[3],cp1[3],cp2[3],cp3[3];
  char      i;

  for (i=0; i<3; i++) {
    cp0[i] = pp[0]->c[i];
    cp1[i] = pp[1]->c[i];
    cp2[i] = pp[2]->c[i];
    cp3[i] = pp[3]->c[i];
  }
  det = det_3d(cp0,cp1,cp2,cp3);
  if ( fabs(det) < AD_EPSD )  return(0);
  det = 1.0 / det;

  /* barycentric coordinate */
  cb[0] =  det_3d(c,cp1,cp2,cp3) * det;
  cb[1] = -det_3d(c,cp2,cp3,cp0) * det;
  cb[2] =  det_3d(c,cp3,cp0,cp1) * det;
  cb[3] = 1.0 - cb[0] - cb[1] - cb[2];

  return(1);
}

/* return velocity interpolation in v in element iv[3] */
static double vecint_3d(double *un,int *iv,double *cb,double *v) {
  double   *u0,*u1,*u2,*u3,dd;

  u0 = &un[3*(iv[0]-1)+1];
  u1 = &un[3*(iv[1]-1)+1];
  u2 = &un[3*(iv[2]-1)+1];
  u3 = &un[3*(iv[3]-1)+1];

  /* P1 interpolate of the speed */   
  v[0] = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0] + cb[3]*u3[0];
  v[1] = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1] + cb[3]*u3[1];
  v[2] = cb[0]*u0[2] + cb[1]*u1[2] + cb[2]*u2[2] + cb[3]*u3[2];
  dd   = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  return(dd);
}


/* find element containing c, starting from nsdep, return baryc. coord */
static int locelt_3d(pMesh mesh,int nsd,double *c,double *cb) {
  pTetra   pt;
  pPoint   p0,p1,p2,p3;
  double   bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
  double   eps,vto,vol1,vol2,vol3,vol4,dd; 
  int      i,nsf,nsp;

  nsf = nsd;
  nsp = nsd;
  ++mesh->mark;
  while ( nsf > 0 ) {
    pt = &mesh->tetra[nsf];
    if ( pt->mark == mesh->mark )  return(nsp);
    pt->mark = mesh->mark;
    
    /* measure of element */
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    p3 = &mesh->point[pt->v[3]];

    /* barycentric and volume */
    bx  = p1->c[0] - p0->c[0];
    by  = p1->c[1] - p0->c[1];
    bz  = p1->c[2] - p0->c[2];
    cx  = p2->c[0] - p0->c[0];
    cy  = p2->c[1] - p0->c[1];
    cz  = p2->c[2] - p0->c[2];
    dx  = p3->c[0] - p0->c[0];
    dy  = p3->c[1] - p0->c[1];
    dz  = p3->c[2] - p0->c[2];

    /* test volume */
    vx  = cy*dz - cz*dy;
    vy  = cz*dx - cx*dz;
    vz  = cx*dy - cy*dx;
    vto = bx*vx + by*vy + bz*vz;
		eps = AD_EPS*vto;

    /* barycentric */
    apx = c[0] - p0->c[0];
    apy = c[1] - p0->c[1];
    apz = c[2] - p0->c[2];

    /* p in 2 */
    vol2  = apx*vx + apy*vy + apz*vz;
    if ( vol2 < eps ) {
      nsp = nsf;
      nsf = pt->adj[1] / 4;
      if ( !nsf ) {
        cb[1] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    /* p in 3 */
    vx  = by*apz - bz*apy;
    vy  = bz*apx - bx*apz;
    vz  = bx*apy - by*apx;
    vol3 = dx*vx + dy*vy + dz*vz;
    if ( vol3 < eps ) {
			nsp = nsf;
      nsf = pt->adj[2] / 4;
      if ( !nsf ) {
        cb[2] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    /* p in 4 */
    vol4 = -cx*vx - cy*vy - cz*vz;
    if ( vol4 < eps ) {
			nsp = nsf;
      nsf = pt->adj[3] / 4;
      if ( !nsf ) {
        cb[3] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    /* p in 1 */
    vol1 = vto - vol2 - vol3 - vol4;
    if ( vol1 < eps ) {
			nsp = nsf;
      nsf = pt->adj[0] / 4;
      if ( !nsf ) {
        cb[0] = 0.0;
        nsf = nsp;
      }
      else
        continue;
    }
    dd = fabs(vol1+vol2+vol3+vol4);
    if ( dd > AD_EPSD ) {
	    dd = 1.0 / dd;
      cb[0] = fabs(vol1) * dd;
      cb[1] = fabs(vol2) * dd;
      cb[2] = fabs(vol3) * dd;
      cb[3] = fabs(vol4) * dd;
    }
    return(nsf);
  }

  return(nsp);
}


/* computes the characteristic line emerging from point with barycentric coordinates cb
 in tetra iel and follows the characteristic on a length dt at most. Most often has to 
 cross the boundary of the current tetra, thus stores in it the new tetra, in cb
 the barycentric coordinates in the new tetra of the crossing point, and updates dt 
 with the remaining time to follow characteristic line */
static int travel_3d(ADst *adst,double *cb,int *iel,double *dt) {
  pTetra      pt;
  pPoint      p[4];
  double     *u0,*u1,*u2,*u3,m[4],dd,ddt,tol,ux,uy,uz,c[3],cb1[4];
  int         k;
  char        i,i0,i1,i2,i3;

  tol  = *dt;
  k    = *iel;
  pt   = &adst->mesh.tetra[k];
  
  p[0] = &adst->mesh.point[pt->v[0]];
  p[1] = &adst->mesh.point[pt->v[1]];
  p[2] = &adst->mesh.point[pt->v[2]];
  p[3] = &adst->mesh.point[pt->v[3]];

  /* the speed at each vertex of iel */
  u0 = &adst->sol.u[3*(pt->v[0]-1)+1];
  u1 = &adst->sol.u[3*(pt->v[1]-1)+1];
  u2 = &adst->sol.u[3*(pt->v[2]-1)+1];
  u3 = &adst->sol.u[3*(pt->v[3]-1)+1];

  /* u = P1 velocity at the point of barycentric coordinates cb */
  ux = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0] + cb[3]*u3[0];
  uy = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1] + cb[3]*u3[1];
  uz = cb[0]*u0[2] + cb[1]*u1[2] + cb[2]*u2[2] + cb[3]*u3[2];  
  if ( ux*ux+uy*uy+uz*uz < AD_EPSD )  return(0);
  
  /* endpoint of the characteristic line starting from cb */
  /* segment is of norm sqrt(ux*ux+uy*uy+uz*uz), and is to be followed for time dt, in the - direction */ 
  c[0] = (cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0] + cb[3]*p[3]->c[0]) - ux;
  c[1] = (cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1] + cb[3]*p[3]->c[1]) - uy;
  c[2] = (cb[0]*p[0]->c[2] + cb[1]*p[1]->c[2] + cb[2]*p[2]->c[2] + cb[3]*p[3]->c[2]) - uz;

  /* barycentric coordinate of uu in the current tetra */
  if ( !bar_3d(p,c,k,cb1) )  return(0);

  /* check if endpoint is in k */
  for (i=0; i<4; i++)
    if ( cb1[i] < 0.0 )  break;
  if ( i == 4 ) {
    memcpy(cb,cb1,4*sizeof(double));
    *dt -= tol;
    return(*dt > 0.0);
  }

  /* ddt = least value of cb[i]/m[i] = (weighted) distance to lambda_i0 = 0 */
  /* ddt = interval of time elapsed before going out of the current element */
  /* condition m[i] > 0 ensures going out of element in "marching direction"*/  
  ddt = AD_TGV;
  i0  = -1;
  for (i=0; i<4; i++) {
    m[i] = cb[i] - cb1[i];
    if ( m[i] > 0.0 ) {
      if ( tol*cb[i]/m[i] < ddt ) {
        ddt = tol*cb[i]/m[i];
        i0  = i;
      }
    }
  }

  /* incomplete advection remain in element */
  if ( ddt > tol ) {
    for (i=0; i<4; i++)  cb[i] = cb[i] - tol*m[i];
    *dt -= tol;
  }
  /* advection goes out of element: advect a minimum value */
  if ( ddt < AD_EPS*tol ) {
    c[0] = cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0] - AD_EPS*ux;
    c[1] = cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1] - AD_EPS*uy;
    c[2] = cb[0]*p[0]->c[2] + cb[1]*p[1]->c[2] + cb[2]*p[2]->c[2] - AD_EPS*uz;
    /* find the new element */
    k = locelt_3d(&adst->mesh,k,c,cb); 
    if ( k < 1 )  return(0);
    *iel = k;
    *dt -= AD_EPS;
  }
  else {
    i1 = (i0+1) % 4;
    i2 = (i1+1) % 4;
    i3 = (i2+1) % 4;

    /* barycentric coordinates of the exit point */ 
    cb1[i0] = 0.0;
    cb1[i1] = cb[i1] - ddt * m[i1];
    cb1[i2] = cb[i2] - ddt * m[i2];
    cb1[i3] = 1.0 - cb1[i1] - cb1[i2];
    if ( cb1[i1] < AD_EPS2 ) {
	    cb1[i2] -= (AD_EPS2-cb1[i1]) / 2;
      cb1[i3] -= (AD_EPS2-cb1[i1]) / 2;
      cb1[i1]  =  AD_EPS2;
    }
    else if ( cb1[i1] > 1.0-AD_EPS2 ) {
	    cb1[i2] += (cb1[i1]-1.0+AD_EPS2) / 2;
      cb1[i3] += (cb1[i1]-1.0+AD_EPS2) / 2;
      cb1[i1]  = 1.0-AD_EPS2;
    }
    if ( cb1[i2] < AD_EPS2 ) {
	    cb1[i3] -= (AD_EPS2-cb1[i2]) / 2;
      cb1[i1] -= (AD_EPS2-cb1[i2]) / 2;
      cb1[i2]  = AD_EPS2;
    }
    else if ( cb1[i2] > 1.0-AD_EPS2 ) {
	    cb1[i3] += (cb1[i2]-1.0+AD_EPS2) / 2;
      cb1[i1] += (cb1[i2]-1.0+AD_EPS2) / 2;
      cb1[i2]  = 1.0-AD_EPS2;
    }
    *dt -= ddt;

    /* coordinates of the point of barycentric coordinates cb1 */
    c[0] = cb1[0]*p[0]->c[0] + cb1[1]*p[1]->c[0] + cb1[2]*p[2]->c[0] + cb1[3]*p[3]->c[0];
    c[1] = cb1[0]*p[0]->c[1] + cb1[1]*p[1]->c[1] + cb1[2]*p[2]->c[1] + cb1[3]*p[3]->c[1];
    c[2] = cb1[0]*p[0]->c[2] + cb1[1]*p[1]->c[2] + cb1[2]*p[2]->c[2] + cb1[3]*p[3]->c[2];

    if ( !pt->adj[i0] ) {
      memcpy(cb,cb1,4*sizeof(double));
      return(0);
    }
    else {
      *iel = pt->adj[i0] / 4;
      if ( !bar_3d(p,c,*iel,cb) )  return(0);
    }
  }

  return(*dt > 0.0);
}


/* 4th order Runge-Kutta: for backtracking characteristic line, step is <0 */ 
/* v = initial speed at point c */
static int nxtptR_3d(ADst *adst,int *iel,double *c,double *cb,double step,double *v) {
  double  h6,v1[3],v2[3],v3[3];
  double  xp1[3],xp2[3],xp3[3],cc[3];
  int     k;

  k = *iel;
  xp1[0] = c[0] - 0.5*step*v[0];
  xp1[1] = c[1] - 0.5*step*v[1];
  xp1[2] = c[2] - 0.5*step*v[2];
  k = locelt_3d(&adst->mesh,k,xp1,cb);
  if ( k < 1 )   return(k);
  vecint_3d(adst->sol.u,adst->mesh.tetra[k].v,cb,v1);

  xp2[0] = c[0] - 0.5*step*v1[0];
  xp2[1] = c[1] - 0.5*step*v1[1];
  xp2[2] = c[2] - 0.5*step*v1[2];
  k = locelt_3d(&adst->mesh,k,xp2,cb);
  if ( k < 1 )   return(k);
  vecint_3d(adst->sol.u,adst->mesh.tetra[k].v,cb,v2);

  xp3[0] = c[0] - step*v2[0];
  xp3[1] = c[1] - step*v2[1];
  xp3[2] = c[2] - step*v2[2];
  k = locelt_3d(&adst->mesh,k,xp3,cb);
  if ( k < 1 )   return(k);
  vecint_3d(adst->sol.u,adst->mesh.tetra[k].v,cb,v3);

  /* final RK4 step */  
  h6 = step / 6.0;
  cc[0] = c[0] - h6 * (v[0] + 2.0*(v1[0] + v2[0]) + v3[0]);
  cc[1] = c[1] - h6 * (v[1] + 2.0*(v1[1] + v2[1]) + v3[1]);
  cc[2] = c[2] - h6 * (v[2] + 2.0*(v1[2] + v2[2]) + v3[2]);
  k = locelt_3d(&adst->mesh,k,cc,cb);
  if ( k < 1 )   return(k);

  /* update */
  c[0] = cc[0];
  c[1] = cc[1];
  c[2] = cc[2];
  vecint_3d(adst->sol.u,adst->mesh.tetra[k].v,cb,v);

  *iel = k;
	return(1);
}


/* Euler scheme for backtracking characteristic line */
static int nxtptE_3d(ADst *adst,int *iel,double *c,double *cb,double step,double *v) {
  double  cc[3];
  int     k;

  k = *iel;
  cc[0] = c[0] - step*v[0];
  cc[1] = c[1] - step*v[1];
  cc[2] = c[2] - step*v[2];
  k = locelt_3d(&adst->mesh,k,cc,cb);
  if ( k < 1 )   return(k);

  /* update */  
  c[0] = cc[0];
  c[1] = cc[1];
  c[2] = cc[2];
  vecint_3d(adst->sol.u,adst->mesh.tetra[k].v,cb,v);
  *iel = k;

  return(1);
}


static void savedt(double dt) {
  FILE   *out;
  
  out = fopen(".dt","w");
  fprintf(out,"%g\n",dt);
  fclose(out);
}


/* solve advection, solution in rv */
int advec1_3d(ADst *adst) {
  pTetra   pt,pt1;
  pPoint   ppt;
  double   cb[4],v[3],c[3],dt,dte,norm,tol,step;  
  int      j,k,ip,iel,kprv,nt,nstep;
  char     i;

  if ( !adst->sol.new ) {
    adst->sol.new = (double *)malloc((adst->info.np+1)*sizeof(double));
    assert(adst->sol.new);
    memcpy(adst->sol.new,adst->sol.chi,(adst->info.np+1)*sizeof(double));
  }

  /* check mesh size and velocity */
  dt = adst->sol.hmin / adst->sol.umax;
  if ( adst->sol.dt < 0.0 ) {
    adst->sol.dt = dt;
    savedt(dt);
  }
  else if ( dt < adst->sol.dt/10.0 ) {
    adst->sol.dt = dt;
    savedt(dt);
  }
  
  /* adjust iterations */
  dt    = adst->sol.dt;
  tol   = adst->sol.dt / 10;
  nstep = (int)(dt/tol);
  step  = dt / nstep;

  if ( adst->info.verb != '0' ) {
    fprintf(stdout,"    Time stepping: %g\n",dt);
    fprintf(stdout,"    Solving: "); fflush(stdout);
  }

  ++adst->mesh.mark;
  nt = 0;
  for (k=1; k<=adst->info.ne; k++) {
    pt = &adst->mesh.tetra[k];
    for (i=0; i<4; i++) {
      ip  = pt->v[i];
      ppt = &adst->mesh.point[ip];
      /* check if already processed */
      if ( ppt->flag == 1 )  continue;
      ppt->flag = 1;
      kprv = k;

      /* coordinate and velocity at starting point */
      v[0] = adst->sol.u[3*(ip-1)+1];
      v[1] = adst->sol.u[3*(ip-1)+2];
      v[2] = adst->sol.u[3*(ip-1)+3];
      norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      if ( norm < AD_EPSD )  continue;

      /* barycentric coordinates of point p in triangle k */
      memset(cb,0,4*sizeof(double));
			cb[i] = 1.0;

      /* next point = foot of the characteristic line */
      c[0] = ppt->c[0];
      c[1] = ppt->c[1];
      c[2] = ppt->c[2];
      dte  = dt;
      for (iel=k,j=0; j<nstep; j++) {
        kprv = iel;
        if ( nxtptR_3d(adst,&iel,c,cb,step,v) < 1 )  break;
        dte -= step;
      }
      if ( j < nstep ) {
        iel = kprv;
        while ( travel_3d(adst,cb,&iel,&dte) );
      }
      /* check if characteristic remains inside domain */
      if ( dte > AD_EPS ) {
        if ( iel < 0 )  iel = kprv;
        iel = locelt_3d(&adst->mesh,iel,c,cb);
				if ( iel < 1 )  iel = kprv;
      }
      /* interpolate value at foot  */
      if ( iel == 0 )  return(0);
      pt1 = &adst->mesh.tetra[iel];
      adst->sol.new[ip] = cb[0]*adst->sol.chi[pt1->v[0]] \
                        + cb[1]*adst->sol.chi[pt1->v[1]] \
                        + cb[2]*adst->sol.chi[pt1->v[2]] \
                        + cb[3]*adst->sol.chi[pt1->v[3]];
      nt++;
    }
  }

  if ( adst->info.verb != '0' )
    fprintf(stdout,"%d characteristics\n",nt);

  return(1);
}