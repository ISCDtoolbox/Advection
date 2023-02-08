#include "advect.h"

extern int ddb;

static double det_3vec(double *v1,double *v2,double *v3) {
  double  det;
  
  det = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1]) - v1[1]*(v2[0]*v3[2]-v2[2]*v3[0]) + v1[2]*(v2[0]*v3[1]-v2[1]*v3[0]);
  
  return(det);
}

/* Barycentric coordinates of the orthogonal projection of point c[] in iel=p0,p1,p2 */
static inline int bar_s(pPoint p[3],double *c,double *cb) {
  double    norm2,n[3],t1[3],t2[3],tc[3];
  
  t1[0] = p[1]->c[0]-p[0]->c[0] ;   t2[0] = p[2]->c[0]-p[0]->c[0] ; tc[0] = c[0] - p[0]->c[0];
  t1[1] = p[1]->c[1]-p[0]->c[1] ;   t2[1] = p[2]->c[1]-p[0]->c[1] ; tc[1] = c[1] - p[0]->c[1];
  t1[2] = p[1]->c[2]-p[0]->c[2] ;   t2[2] = p[2]->c[2]-p[0]->c[2] ; tc[2] = c[2] - p[0]->c[2];

  /* Non normalized normal vector to triangle */
  n[0] = t1[1]*t2[2] - t1[2]*t2[1];
  n[1] = t1[2]*t2[0] - t1[0]*t2[2];
  n[2] = t1[0]*t2[1] - t1[1]*t2[0];

  norm2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( norm2 < AD_EPSD ) return(0);

  norm2 = 1.0 / norm2;

  /* Barycentric coordinate */
  cb[1] = det_3vec(tc,t2,n);
  cb[2] = det_3vec(t1,tc,n);
  cb[1] *= norm2;
  cb[2] *= norm2;

  cb[0] = 1.0 - cb[1] - cb[2];

  return(1);
}

/* Calculate the (point) ball of point i in triangle k */
static int boulep_s(pMesh mesh,int start,char i,int *list) {
  pTria       pt;
  int         ilist,k;
  char        i0,i1,i2,voy;
  
  ilist = 0;
  k = start;
  i0 = i;
  
  do {
    pt = &mesh->tria[k];
    i1 = (i0+1)%3;
    
    list[ilist] = pt->v[i1];
    ilist++;
    if ( ilist > AD_LONMAX-2 ) return(-ilist);
    k = pt->adj[i1] / 3;
    voy = pt->adj[i1] % 3;
    i0 = (voy+1)%3;
  }
  while ( k && k != start );
  
  if ( k == start )  return(ilist);
  
  /* Add last point to list */
  i2 = (i1+1)%3;
  list[ilist] = pt->v[i2];
  ilist++;
  if ( ilist > AD_LONMAX-2 ) return(-ilist);
  
  /* Reverse loop */
  k = start;
  pt = &mesh->tria[k];
  i0 = i;
  i2 = (i0+2)%3;
  k  = pt->adj[i2] / 3;
  i1 = pt->adj[i2] % 3;
  i2 = (i1+1)%3;
  
  while ( k ) {
    pt = &mesh->tria[k];
    list[ilist] = pt->v[i1];
    ilist++;
    if ( ilist > AD_LONMAX-2 ) return(-ilist);
    
    k  = pt->adj[i2] / 3;
    i1 = pt->adj[i2] % 3;
    i2 = (i1+1)%3;
  }
  
  return(ilist);
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


/* Unit normal vector of triangle k */
static int nortri(pMesh mesh,int k,double *n) {
  pTria    pt;
  pPoint   p0,p1,p2;
  double   norm;
  
  pt = &mesh->tria[k];
  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];
  
  n[0] = (p1->c[1]-p0->c[1])*(p2->c[2]-p0->c[2]) - (p1->c[2]-p0->c[2])*(p2->c[1]-p0->c[1]);
  n[1] = (p1->c[2]-p0->c[2])*(p2->c[0]-p0->c[0]) - (p1->c[0]-p0->c[0])*(p2->c[2]-p0->c[2]);
  n[2] = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);
  
  norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if ( norm > AD_EPS2 ) norm = 1.0 / norm;
  
  n[0] *= norm;
  n[1] *= norm;
  n[2] *= norm;
  
  return(1);
}

/* Unit normal vector at vertex i in triangle k */
static int norver(pMesh mesh,int k,int i,double *n) {
  pTria    pt;
  pPoint   p0,p1,p2;
  double   norm;
  int      l,iel,ilist,list[AD_LONMAX];
  
  memset(n,0,3*sizeof(double));
  
  ilist = boulet_s(mesh,k,i,list);
  
  /* Average of the unit normal vector to the triangles surrounding i (weighted by the triangle areas) */
  for (l=0; l<ilist; l++) {
    iel = list[l] / 3;
    pt  = &mesh->tria[iel];
    
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    
    n[0] += (p1->c[1]-p0->c[1])*(p2->c[2]-p0->c[2]) - (p1->c[2]-p0->c[2])*(p2->c[1]-p0->c[1]);
    n[1] += (p1->c[2]-p0->c[2])*(p2->c[0]-p0->c[0]) - (p1->c[0]-p0->c[0])*(p2->c[2]-p0->c[2]);
    n[2] += (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);
  }
  
  norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  if ( norm > AD_EPS2 ) norm = 1.0 / norm;
  
  n[0] *= norm;
  n[1] *= norm;
  n[2] *= norm;
  
  return(1);
}

/* Calculate vector v in normal plane to n1 such that vw = orthogonal projection of v onto normal plane to n2 */
static int invproj(double *n1,double *n2,double *w,double *v) {
  double   n1n2,n1w,dd;
  
  n1n2 = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
  n1w  = n1[0]*w[0] + n1[1]*w[1] + n1[2]*w[2];

  if ( fabs(n1n2) < AD_EPS ) return(0);
  dd   = n1w / n1n2;

  v[0] = w[0] + dd*n2[0];
  v[1] = w[1] + dd*n2[1];
  v[2] = w[2] + dd*n2[2];

  return(1);
}

/* Tangential projection of vector field v onto the plane of triangle k */
static int tanproj(pMesh mesh,int k,double *v) {
  pTria     pt;
  pPoint    p0,p1,p2;
  double    n[3],det,ps,norm;
  
  pt = &mesh->tria[k];
  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];
  
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

/* Perform one step of mesh travel to backtrack characteristic line of u,
   stemming fro bcs cb in iel, for a total time dt.
   Return 1 when successful step: iel, cb are updated, and fraction of time dt updated
   Return 0 when either advection complete,
                        or an external boundary is met: iel and cb then contain the last triangle, with the cbs of the exit point */
static int travel_s(ADst *adst,double *cb,int *iel,double *dt) {
  pTria     pt;
  pPoint    p[3];
  double    *u0,*u1,*u2,tol,ps1,ps2,pps,psmax,t1[3],t2[3],ta[3],n[3],u[3],c[3],m[3],cb1[3];
  double    det,ddt,ddt1,ddt2,norm;
  int       nz,k,kmax,jel,kel,ilist,list[AD_LONMAX];
  char      i,i0,i1,i2,iim,imax,j0,j1,j2,ia,je;

  tol = *dt;
  pt = &adst->mesh.tria[*iel];

  nz = 0;
  for (i=0; i<3; i++) {
    if ( cb[i] < AD_EPS ) {
      nz++;
      ia = i;
    }
    else
      i0 = i;
  }
    
  /* Advection starts from vertex i0 of pt */
  if ( nz == 2 ) {
    if ( ddb ) printf("Cas 1 \n");

    u0 = &adst->sol.u[3*(pt->v[i0]-1)+1];
    ilist = boulet_s(&adst->mesh,*iel,i0,list);
    
    /* Identify initial triangle for advection by u */
    psmax = -AD_TGV;
    for (k=0; k<ilist; k++) {
      jel = list[k] / 3;
      j0  = list[k] % 3;
      j1  = (j0+1) % 3;
      j2  = (j0+2) % 3;

      pt = &adst->mesh.tria[jel];
      p[0] = &adst->mesh.point[pt->v[0]];
      p[1] = &adst->mesh.point[pt->v[1]];
      p[2] = &adst->mesh.point[pt->v[2]];
      
      /* Endpoint of the tangent line to the characteristic */
      c[0]  = p[j0]->c[0] - tol*u0[0];
      c[1]  = p[j0]->c[1] - tol*u0[1];
      c[2]  = p[j0]->c[2] - tol*u0[2];
      
      if ( !bar_s(p,c,cb1) ) continue;
      
      /* -u0 goes outside pt, on the side of edge j2 */
      if ( cb1[j1] < cb1[j2] ) {
        pps = cb1[j1];
        iim = j2;
      }
  
      /* -u0 goes outside pt, on the side of edge i1 (=t2) */
      else {
        pps = cb1[j2];
        iim = j1;
      }

      if ( pps > 0.0 ) break;
      
      /* If no triangle of the ball contains -u0, kmax = that closest to contain it (imax = index of defect edge) */
      if ( pps > psmax ) {
        psmax = pps;
        kmax = k;
        imax = iim;
      }
    }

    /* Vector -u0 falls into triangle jel: exit through edge j0 */
    if ( k < ilist ) {
      if ( pt->mark == adst->mesh.mark )  return(0);
    
      cb[j0] = 1.0;
      cb[j1] = 0.0;
      cb[j2] = 0.0;
            
      /* ddt = time needed to exit triangle through edge j0 */
      for (i=0; i<3; i++)
        m[i] = cb[i] - cb1[i];
      
      if ( m[j0] < AD_EPSD  ) return(0);
    
      ddt = tol * 1.0 / m[j0];
      
      /* Case when advection stays in the same triangle */
      if ( ddt > tol ) {
        *iel = jel;
        memcpy(cb,cb1,3*sizeof(double));
        *dt -= tol;
        return(*dt>0.0);
      }
      else {
        /* Barycentric coordinates of the exit point */
        for (i=0; i<3; i++)
          cb1[i] = cb[i] - ddt * m[i] / tol;
        *dt -= ddt;
        
        kel = pt->adj[j0] / 3;
    
        /* External boundary met */
        if ( !kel ) {
          memcpy(cb,cb1,3*sizeof(double));
          return(0);
        }
        /* Pass to the next triangle */
        else {
          pt->mark = adst->mesh.mark;
          *iel = kel;
          i0   = pt->adj[j0] % 3;
          cb[i0] = 0.0;
          cb[(i0+1)%3] = cb1[j2];
          cb[(i0+2)%3] = cb1[j1];
          return(*dt>0.0);
        }
      }
    }
    /* No triangle contains (projection of) u0: snapping to edge imax of kmax */
    else {
      return(0);
      jel = list[kmax] / 3;
      j0  = list[kmax] % 3;
      je  = imax;
                  
      pt = &adst->mesh.tria[jel];
      c[0]  = p[j0]->c[0] - tol*u0[0];
      c[1]  = p[j0]->c[1] - tol*u0[1];
      c[2]  = p[j0]->c[2] - tol*u0[2];
      
      if ( !bar_s(p,c,cb1) ) return(0);
      pt->mark = adst->mesh.mark;
            
      /* Endpoint falls in current triangle */
      if ( cb1[j0]  > 0.0 ) {
        cb1[je] = 0.0;
        cb1[j0]  = 1.0 - cb1[(j0+1)%3] - cb1[(j0+2)%3];
        *iel = jel;
        memcpy(cb,cb1,3*sizeof(double));
        *dt -= tol;
        return(*dt>0.0);
      }
      
      /* Otherwise, stop at third vertex */
      m[j0]  = 1.0 - cb1[j0];
      if ( m[j0] < AD_EPSD  ) return(0);
      
      ddt = tol * 1.0 / m[j0];
      
      *iel = jel;
      for (i=0; i<3; i++) cb[i] = 1.0;
      cb[j0] = 0.0;
      cb[je]  = 0.0;
      
      *dt -= ddt;
      return(*dt > 0.0);
    }
  }

  /* Travel starts from inside edge ia = j1j2 */
  if ( nz == 1 ) {
    if ( ddb )  printf("Cas 2 \n");
    i0 = ia;
    i1 = (ia+1) % 3;
    i2 = (ia+2) % 3;

    p[0] = &adst->mesh.point[pt->v[0]];
    p[1] = &adst->mesh.point[pt->v[1]];
    p[2] = &adst->mesh.point[pt->v[2]];

    /* Interpolation of velocity field */
    u1 = &adst->sol.u[3*(pt->v[i1]-1)+1];
    u2 = &adst->sol.u[3*(pt->v[i2]-1)+1];
    
    u[0] = cb[i1]*u1[0] + cb[i2]*u2[0];
    u[1] = cb[i1]*u1[1] + cb[i2]*u2[1];
    u[2] = cb[i1]*u1[2] + cb[i2]*u2[2];

    /* Endpoint of the tangent line to the characteristic */
    c[0] = cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0] - tol*u[0];
    c[1] = cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1] - tol*u[1];
    c[2] = cb[0]*p[0]->c[2] + cb[1]*p[1]->c[2] + cb[2]*p[2]->c[2] - tol*u[2];
    
    if ( !bar_s(p,c,cb1) ) return(0);
    
    /* Stay in the same triangle */
    if ( cb1[i0] > 0.0 ) {
      if ( ddb ) printf("Cas 2 - 1 \n");

      if ( pt->mark == adst->mesh.mark ) {
        // printf("Ou bien ici ?\n");
        return(0);
      }
      /* check if endpoint is in k */
      for (i=0; i<3; i++)
        if ( cb1[i] < 0.0 )  break;
      if ( i == 3 ) {
        memcpy(cb,cb1,3*sizeof(double));
        *dt -= tol;
        return(*dt>0.0);
      }

      /* Identify exit edge (j1 or j2) */
      ddt1 = AD_TGV;
      ddt2 = AD_TGV;
      for (i=0; i<3; i++)
        m[i] = cb[i] - cb1[i];
      if ( m[i1] > 0.0 ) ddt1 = tol * cb[i1] / m[i1];
      if ( m[i2] > 0.0 ) ddt2 = tol * cb[i2] / m[i2];
      
      pt->mark = adst->mesh.mark;
      
      /* Exit through i1 */
      if ( ddt1 < ddt2 ) {
        for (i=0; i<3; i++)
          cb1[i] = cb[i] - ddt1 * m[i] / tol;
        *dt -= ddt1;
        jel = pt->adj[i1] / 3;
        
        if ( !jel ) {
          memcpy(cb,cb1,3*sizeof(double));
          return(0);
        }
        else {
          j0  = pt->adj[i1] % 3;
          cb[j0] = 0.0;
          cb[(j0+1)%3] = cb1[i0];
          cb[(j0+2)%3] = cb1[i2];
          *iel = jel;
          return(*dt>0.0);
        }
      }
      /* Exit through i2 */
      else {
        for (i=0; i<3; i++)
          cb1[i] = cb[i] - ddt2 * m[i] / tol;
        *dt -= ddt2;
        jel = pt->adj[i2] / 3;
        
        if ( !jel ) {
          memcpy(cb,cb1,3*sizeof(double));
          return(0);
        }
        else {
          j0  = pt->adj[i2] % 3;
          cb[j0] = 0.0;
          cb[(j0+1)%3] = cb1[i1];
          cb[(j0+2)%3] = cb1[i0];
          *iel = jel;
          return(*dt>0.0);
        }
      }
    }
    /* Change triangle */
    else {
      if ( ddb ) printf("Cas 2 - 2 \n");
      jel = pt->adj[ia] / 3;
      j0  = pt->adj[ia] % 3;
      j1  = (j0+1) % 3;
      j2  = (j0+2) % 3;
      if ( !jel ) return(0);

      cb1[j0] = cb[i0];
      cb1[j1] = cb[i2];
      cb1[j2] = cb[i1];
      memcpy(cb,cb1,3*sizeof(double));
      
      pt = &adst->mesh.tria[jel];
      if ( pt->mark == adst->mesh.mark ) return(0);
  
      p[0] = &adst->mesh.point[pt->v[0]];
      p[1] = &adst->mesh.point[pt->v[1]];
      p[2] = &adst->mesh.point[pt->v[2]];

      if ( !bar_s(p,c,cb1) ) return(0);
            
      /* check if endpoint is in jel */
      for (i=0; i<3; i++)
        if ( cb1[i] < 0.0 )  break;
      if ( i == 3 ) {
        memcpy(cb,cb1,3*sizeof(double));
        *dt -= tol;
        return(*dt > 0.0);
      }
            
      /* Identify exit edge (j1 or j2) */
      ddt1 = AD_TGV;
      ddt2 = AD_TGV;
      for (i=0; i<3; i++)
        m[i] = cb[i] - cb1[i];
      
      if ( m[j1] > 0.0 ) ddt1 = tol * cb[j1] / m[j1];
      if ( m[j2] > 0.0 ) ddt2 = tol * cb[j2] / m[j2];
            
      pt->mark = adst->mesh.mark;
      
      /* Exit through j1 */
      if ( ddt1 < ddt2 ) {
        for (i=0; i<3; i++)
          cb1[i] = cb[i] - ddt1 * m[i] / tol;
        *dt -= ddt1;
        kel = pt->adj[j1] / 3;
                
        if ( !kel ) {
          memcpy(cb,cb1,3*sizeof(double));
          return(0);
        }
        else {
          i0  = pt->adj[j1] % 3;
          cb[i0] = 0.0;
          cb[(i0+1)%3] = cb1[j0];
          cb[(i0+2)%3] = cb1[j2];
          *iel = kel;
          return(*dt>0.0);
        }
      }
      /* Exit through j2 */
      else {
        for (i=0; i<3; i++)
          cb1[i] = cb[i] - ddt2 * m[i] / tol;
        *dt -= ddt2;
        kel = pt->adj[j2] / 3;
        
        if ( !kel ) {
          memcpy(cb,cb1,3*sizeof(double));
          return(0);
        }
        else {
          i0  = pt->adj[j2] % 3;
          cb[i0] = 0.0;
          cb[(i0+1)%3] = cb1[j1];
          cb[(i0+2)%3] = cb1[j0];
          *iel = kel;
          return(*dt>0.0);
        }
      }
    }
  }

  /* Default: travel starts from inside pt */
  if ( ddb ) printf("Cas 3 \n");

  k = *iel;
  pt = &adst->mesh.tria[k];
  p[0] = &adst->mesh.point[pt->v[0]];
  p[1] = &adst->mesh.point[pt->v[1]];
  p[2] = &adst->mesh.point[pt->v[2]];

  /* velocity at each vertex of iel */
  u0 = &adst->sol.u[3*(pt->v[0]-1)+1];
  u1 = &adst->sol.u[3*(pt->v[1]-1)+1];
  u2 = &adst->sol.u[3*(pt->v[2]-1)+1];

  /* u = (projected) P1 velocity at the point with barycentric coordinates cb */
  u[0] = cb[0]*u0[0] + cb[1]*u1[0] + cb[2]*u2[0];
  u[1] = cb[0]*u0[1] + cb[1]*u1[1] + cb[2]*u2[1];
  u[2] = cb[0]*u0[2] + cb[1]*u1[2] + cb[2]*u2[2];
  tanproj(&adst->mesh,k,u); // pas forcement necessaire

  if ( u[0]*u[0] + u[1]*u[1] + u[2]*u[2] < AD_EPSD ) return(0);
  
  /* Endpoint of characteristic line generated by u */
  c[0] = cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0] - tol*u[0];
  c[1] = cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1] - tol*u[1];
  c[2] = cb[0]*p[0]->c[2] + cb[1]*p[1]->c[2] + cb[2]*p[2]->c[2] - tol*u[2];
  
  /* Barycentric coordinates of c in current triangle  */
  if ( !bar_s(p,c,cb1) ) return(0);
  
  /* Case where endpoint is in k */
  for (i=0; i<3; i++)
    if ( cb1[i] < 0.0 )  break;
  if ( i == 3 ) {
    memcpy(cb,cb1,3*sizeof(double));
    pt->mark = adst->mesh.mark;
    *dt -= tol;
    return(*dt > 0.0);
  }
  
  /* ddt = needed time to exit triangle ; i0 = exit edge */
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
  
  /* Barycentric coordinates of the exit point */
  for (i=0; i<3; i++)
    cb1[i] = cb[i] - ddt * m[i] / tol;
  *dt -= ddt;
  pt->mark = adst->mesh.mark;
  
  /* Find output triangle and bcs of the output point */
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
  
  return(*dt > 0.0);
}

/* Perform geodesic backtracking in one single triangle, from point cb in iel, with velocity v, for a time dt */
static int geod_step(ADst *adst,int *iel,double *cb,double *dt, double *v) {
  pTria     pt;
  pPoint    p[3];
  double    pps,psmax,dd,normv,normvt,ddt,ddt1,ddt2,tol,c[3],cb1[3],m[3],vt[3];
  int       nz,k,kmax,jel,kel,ilist,list[AD_LONMAX];
  char      i,iim,imax,ia,i0,i1,i2,j0,j1,j2;
  
  tol = *dt;
  ++adst->mesh.mark;

  pt = &adst->mesh.tria[*iel];
  normv = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  nz = 0;
  for (i=0; i<3; i++) {
    if ( cb[i] < AD_EPS ) {
      nz++;
      ia = i;
    }
    else
      i0 = i;
  }
      
  /* Case 1: Advection starts from vertex i0 of pt */
  if ( nz == 2 ) {
    if  ( ddb ) printf("Cas 2\n");
    ilist = boulet_s(&adst->mesh,*iel,i0,list);
        
    /* Identify initial triangle for advection by v */
    psmax = -AD_TGV;
    for (k=0; k<ilist; k++) {
      jel = list[k] / 3;
      j0  = list[k] % 3;
      j1  = (j0+1) % 3;
      j2  = (j0+2) % 3;

      pt = &adst->mesh.tria[jel];
      p[0] = &adst->mesh.point[pt->v[0]];
      p[1] = &adst->mesh.point[pt->v[1]];
      p[2] = &adst->mesh.point[pt->v[2]];
      
      /* Normalized tangential projection of v */
      memcpy(vt,v,3*sizeof(double));
      tanproj(&adst->mesh,jel,vt);
      normvt = sqrt(vt[0]*vt[0]+vt[1]*vt[1]+vt[2]*vt[2]);
      dd = ( normvt < AD_EPS ) ? normv / (normvt+AD_EPS) : normv / normvt;
      vt[0] *= dd;
      vt[1] *= dd;
      vt[2] *= dd;
      
      /* Endpoint of the tangent line to the characteristic */
      c[0]  = p[j0]->c[0] - tol*vt[0];
      c[1]  = p[j0]->c[1] - tol*vt[1];
      c[2]  = p[j0]->c[2] - tol*vt[2];
      
      if ( !bar_s(p,c,cb1) ) continue;
      
      /* -vt goes outside pt, on the side of edge j2 */
      if ( cb1[j1] < cb1[j2] ) {
        pps = cb1[j1];
        iim = j2;
      }
  
      /* -vt goes outside pt, on the side of edge i1 (=t2) */
      else {
        pps = cb1[j2];
        iim = j1;
      }

      if ( pps > 0.0 ) break;
      
      /* If no triangle of the ball contains -vt, kmax = that closest to contain it (imax = index of defect edge) */
      if ( pps > psmax ) {
        psmax = pps;
        kmax = k;
        imax = iim;
      }
    }
    
    /* No triangle contains (projection of) -vt */
    if ( k == ilist ) return(0);

    /* Vector -vt passes into triangle jel: exit through edge j0 */
    if ( pt->mark == adst->mesh.mark )  return(0);
    cb[j0] = 1.0;
    cb[j1] = 0.0;
    cb[j2] = 0.0;
            
    /* ddt = time needed to exit triangle through edge j0 */
    for (i=0; i<3; i++)
      m[i] = cb[i] - cb1[i];
      
    if ( m[j0] < AD_EPSD  ) return(0);
    ddt = tol * 1.0 / m[j0];
      
    /* Case when advection stays in the same triangle */
    if ( ddt > tol ) {
      *iel = jel;
      memcpy(cb,cb1,3*sizeof(double));
      *dt -= tol;
      return(*dt>0.0);
    }
    else {
      /* Barycentric coordinates of the exit point */
      for (i=0; i<3; i++)
        cb1[i] = cb[i] - ddt * m[i] / tol;
      *dt -= ddt;
        
      kel = pt->adj[j0] / 3;
    
      /* External boundary met */
      if ( !kel ) {
        memcpy(cb,cb1,3*sizeof(double));
        return(0);
      }
      /* Pass to the next triangle */
      else {
        pt->mark = adst->mesh.mark;
        *iel = kel;
        i0   = pt->adj[j0] % 3;
        cb[i0] = 0.0;
        cb[(i0+1)%3] = cb1[j2];
        cb[(i0+2)%3] = cb1[j1];
        return(*dt>0.0);
      }
    }
  }
  
  /* Case 2: Travel starts from inside edge ia = j1j2 */
  if ( nz == 1 ) {
    if  ( ddb ) printf("Cas 1 %d\n",*iel);

    k = *iel;
    pt = &adst->mesh.tria[k];

    i0 = ia;
    i1 = (ia+1) % 3;
    i2 = (ia+2) % 3;

    p[0] = &adst->mesh.point[pt->v[0]];
    p[1] = &adst->mesh.point[pt->v[1]];
    p[2] = &adst->mesh.point[pt->v[2]];

    /* Normalized projection of velocity field in tria iel */
    memcpy(vt,v,3*sizeof(double));
    tanproj(&adst->mesh,k,vt);
    normvt = sqrt(vt[0]*vt[0]+vt[1]*vt[1]+vt[2]*vt[2]);
    dd = ( normvt < AD_EPS ) ? normv / (normvt+AD_EPS) : normv / normvt;
    vt[0] *= dd;
    vt[1] *= dd;
    vt[2] *= dd;
    
    /* Endpoint of the tangent line to the characteristic */
    c[0] = cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0] - tol*vt[0];
    c[1] = cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1] - tol*vt[1];
    c[2] = cb[0]*p[0]->c[2] + cb[1]*p[1]->c[2] + cb[2]*p[2]->c[2] - tol*vt[2];
    
    if ( !bar_s(p,c,cb1) ) return(0);
    
    /* Stay in the same triangle */
    if ( cb1[i0] > 0.0 ) {
      if  ( ddb ) printf("Cas 1 - same \n");

      if ( pt->mark == adst->mesh.mark ) return(0);

      /* Check if endpoint is in elt */
      for (i=0; i<3; i++)
        if ( cb1[i] < 0.0 )  break;
      if ( i == 3 ) {
        memcpy(cb,cb1,3*sizeof(double));
        *dt -= tol;
        return(*dt>0.0);
      }

      /* Identify exit edge (j1 or j2) */
      ddt1 = AD_TGV;
      ddt2 = AD_TGV;
      for (i=0; i<3; i++)
        m[i] = cb[i] - cb1[i];
      if ( m[i1] > 0.0 ) ddt1 = tol * cb[i1] / m[i1];
      if ( m[i2] > 0.0 ) ddt2 = tol * cb[i2] / m[i2];
      
      pt->mark = adst->mesh.mark;
      
      /* Exit through i1 */
      if ( ddt1 < ddt2 ) {
        for (i=0; i<3; i++)
          cb1[i] = cb[i] - ddt1 * m[i] / tol;
        *dt -= ddt1;
        jel = pt->adj[i1] / 3;
        
        if ( !jel ) {
          memcpy(cb,cb1,3*sizeof(double));
          return(0);
        }
        else {
          j0  = pt->adj[i1] % 3;
          cb[j0] = 0.0;
          cb[(j0+1)%3] = cb1[i0];
          cb[(j0+2)%3] = cb1[i2];
          *iel = jel;
          return(*dt>0.0);
        }
      }
      /* Exit through i2 */
      else {
        for (i=0; i<3; i++)
          cb1[i] = cb[i] - ddt2 * m[i] / tol;
        *dt -= ddt2;
        jel = pt->adj[i2] / 3;
        
        if ( !jel ) {
          memcpy(cb,cb1,3*sizeof(double));
          return(0);
        }
        else {
          j0  = pt->adj[i2] % 3;
          cb[j0] = 0.0;
          cb[(j0+1)%3] = cb1[i1];
          cb[(j0+2)%3] = cb1[i0];
          *iel = jel;
          return(*dt>0.0);
        }
      }
    }
    
    /* Change triangle */
    else {
      if  ( ddb ) printf("Cas 1 - change \n");

      jel = pt->adj[ia] / 3;
      j0  = pt->adj[ia] % 3;
      j1  = (j0+1) % 3;
      j2  = (j0+2) % 3;
      if ( !jel ) return(0);

      cb1[j0] = cb[i0];
      cb1[j1] = cb[i2];
      cb1[j2] = cb[i1];
      memcpy(cb,cb1,3*sizeof(double));
      
      pt = &adst->mesh.tria[jel];
      if ( pt->mark == adst->mesh.mark ) return(0);
  
      p[0] = &adst->mesh.point[pt->v[0]];
      p[1] = &adst->mesh.point[pt->v[1]];
      p[2] = &adst->mesh.point[pt->v[2]];

      if ( !bar_s(p,c,cb1) ) return(0);
            
      /* check if endpoint is in jel */
      for (i=0; i<3; i++)
        if ( cb1[i] < 0.0 )  break;
      if ( i == 3 ) {
        memcpy(cb,cb1,3*sizeof(double));
        *dt -= tol;
        return(*dt > 0.0);
      }
            
      /* Identify exit edge (j1 or j2) */
      ddt1 = AD_TGV;
      ddt2 = AD_TGV;
      for (i=0; i<3; i++)
        m[i] = cb[i] - cb1[i];
      
      if ( m[j1] > 0.0 ) ddt1 = tol * cb[j1] / m[j1];
      if ( m[j2] > 0.0 ) ddt2 = tol * cb[j2] / m[j2];
            
      pt->mark = adst->mesh.mark;
      
      /* Exit through j1 */
      if ( ddt1 < ddt2 ) {
        for (i=0; i<3; i++)
          cb1[i] = cb[i] - ddt1 * m[i] / tol;
        *dt -= ddt1;
        kel = pt->adj[j1] / 3;
                
        if ( !kel ) {
          memcpy(cb,cb1,3*sizeof(double));
          return(0);
        }
        else {
          i0  = pt->adj[j1] % 3;
          cb[i0] = 0.0;
          cb[(i0+1)%3] = cb1[j0];
          cb[(i0+2)%3] = cb1[j2];
          *iel = kel;
          return(*dt>0.0);
        }
      }
      /* Exit through j2 */
      else {
        for (i=0; i<3; i++)
          cb1[i] = cb[i] - ddt2 * m[i] / tol;
        *dt -= ddt2;
        kel = pt->adj[j2] / 3;
        
        if ( !kel ) {
          memcpy(cb,cb1,3*sizeof(double));
          return(0);
        }
        else {
          i0  = pt->adj[j2] % 3;
          cb[i0] = 0.0;
          cb[(i0+1)%3] = cb1[j1];
          cb[(i0+2)%3] = cb1[j0];
          *iel = kel;
          return(*dt>0.0);
        }
      }
    }
  }
  
  /* Case 3: travel starts from inside pt */
  if  ( ddb ) printf("Cas 0\n");

  k = *iel;
  pt = &adst->mesh.tria[k];
  p[0] = &adst->mesh.point[pt->v[0]];
  p[1] = &adst->mesh.point[pt->v[1]];
  p[2] = &adst->mesh.point[pt->v[2]];

  /* Normalized projection of velocity field in tria k */
  memcpy(vt,v,3*sizeof(double));
  tanproj(&adst->mesh,k,vt);
  normvt = sqrt(vt[0]*vt[0]+vt[1]*vt[1]+vt[2]*vt[2]);
  dd = ( normvt < AD_EPS ) ? normv / (normvt+AD_EPS) : normv / normvt;
  vt[0] *= dd;
  vt[1] *= dd;
  vt[2] *= dd;
  
  /* Endpoint of characteristic line generated by u */
  c[0] = cb[0]*p[0]->c[0] + cb[1]*p[1]->c[0] + cb[2]*p[2]->c[0] - tol*vt[0];
  c[1] = cb[0]*p[0]->c[1] + cb[1]*p[1]->c[1] + cb[2]*p[2]->c[1] - tol*vt[1];
  c[2] = cb[0]*p[0]->c[2] + cb[1]*p[1]->c[2] + cb[2]*p[2]->c[2] - tol*vt[2];
  
  /* Barycentric coordinates of c in current triangle  */
  if ( !bar_s(p,c,cb1) ) return(0);
  
  /* Case where endpoint is in k */
  for (i=0; i<3; i++)
    if ( cb1[i] < 0.0 )  break;
  if ( i == 3 ) {
    memcpy(cb,cb1,3*sizeof(double));
    pt->mark = adst->mesh.mark;
    *dt -= tol;
    return(*dt > 0.0);
  }
  
  /* ddt = needed time to exit triangle ; i0 = exit edge */
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
  
  /* Barycentric coordinates of the exit point */
  for (i=0; i<3; i++)
    cb1[i] = cb[i] - ddt * m[i] / tol;
  *dt -= ddt;
  pt->mark = adst->mesh.mark;
  
  /* Find output triangle and bcs of the output point */
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
  
  return(*dt > 0.0);
  
}

/* Backtrack geodesic emerging from point cb in iel, with velocity v, for a time dt */
static int geod(ADst *adst,int *iel,double *cb,double dt, double *v) {
  double   dte;
  
  dte = dt;
  
  /* Backtracking geodesic until completion */
  while ( geod_step(adst,iel,cb,&dte,v) );

  /* Incomplete backtracking */
  if ( dte > AD_EPS ) return(0);
  
  return(1);
}

/* Find next point in characteristic backtracking by the 4th order Runge-Kutta method */
static int nxtptR_s(ADst *adst,int *iel,double *cb,double step) {
  pTria     pt;
  pPoint    p[3];
  double    *u0,*u1,*u2,cb1[3],v0[3],v1[3],v2[3],v3[3];
  double    h2,dd,v[3],w[3],n1[3],n[3];
  int       k,nz;
  char      i,ia,i0,i1;
  
  nz = 0;
  h2 = 0.5*step;
    
  /* Calculate normal vector n1 at starting point */
  k = *iel;
  pt = &adst->mesh.tria[k];
  
  for (i=0; i<3; i++) {
    if ( cb[i] < AD_EPS ) {
      nz++;
      ia = i;
    }
    else
      i0 = i;
  }
  
  if ( nz == 2 ) {
    if ( !norver(&adst->mesh,k,i0,n1) ) return(0);
  }
  else if ( nz == 1 ) {
    i0 = (ia+1) % 3;
    i1 = (i0+1) % 3;
    if ( !norver(&adst->mesh,k,i0,n1) ) return(0);
    if ( !norver(&adst->mesh,k,i1,n) ) return(0);
    n1[0] += n[0];
    n1[1] += n[1];
    n1[2] += n[2];
    dd = sqrt( n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2] );
    if ( dd < AD_EPS ) return(0);
    dd = 1.0 / dd;
    
    n1[0] *= dd;
    n1[1] *= dd;
    n1[2] *= dd;
  }
  else {
    if ( !nortri(&adst->mesh,k,n1) ) return(0);
  }
    
  /* First integration point = starting point */
  memcpy(cb1,cb,3*sizeof(double));

  p[0] = &adst->mesh.point[pt->v[0]];
  p[1] = &adst->mesh.point[pt->v[1]];
  p[2] = &adst->mesh.point[pt->v[2]];

  /* velocity at each vertex of k */
  u0 = &adst->sol.u[3*(pt->v[0]-1)+1];
  u1 = &adst->sol.u[3*(pt->v[1]-1)+1];
  u2 = &adst->sol.u[3*(pt->v[2]-1)+1];

  /* v0 = P1 velocity at the point with barycentric coordinates cb */
  v0[0] = cb1[0]*u0[0] + cb1[1]*u1[0] + cb1[2]*u2[0];
  v0[1] = cb1[0]*u0[1] + cb1[1]*u1[1] + cb1[2]*u2[1];
  v0[2] = cb1[0]*u0[2] + cb1[1]*u1[2] + cb1[2]*u2[2];

  dd = v0[0]*n1[0] + v0[1]*n1[1] + v0[2]*n1[2];
  v0[0] -= dd * n1[0];
  v0[1] -= dd * n1[1];
  v0[2] -= dd * n1[2];
  
  /* Second integration point */
  if ( !geod(adst,&k,cb1,h2,v0) ) return(0);
  
  pt = &adst->mesh.tria[k];
  p[0] = &adst->mesh.point[pt->v[0]];
  p[1] = &adst->mesh.point[pt->v[1]];
  p[2] = &adst->mesh.point[pt->v[2]];

  /* velocity at each vertex of k */
  u0 = &adst->sol.u[3*(pt->v[0]-1)+1];
  u1 = &adst->sol.u[3*(pt->v[1]-1)+1];
  u2 = &adst->sol.u[3*(pt->v[2]-1)+1];

  /* w = P1 velocity at the point with barycentric coordinates cb */
  w[0] = cb1[0]*u0[0] + cb1[1]*u1[0] + cb1[2]*u2[0];
  w[1] = cb1[0]*u0[1] + cb1[1]*u1[1] + cb1[2]*u2[1];
  w[2] = cb1[0]*u0[2] + cb1[1]*u1[2] + cb1[2]*u2[2];
  
  /* Transport w onto initial tangent plane to obtain v1 */
  if ( !nortri(&adst->mesh,k,n) ) return(0);
  if ( !invproj(n1,n,w,v1) ) return(0);
  
  /* Third integration point */
  k = *iel;
  memcpy(cb1,cb,3*sizeof(double));

  if ( !geod(adst,&k,cb1,h2,v1) ) return(0);
  
  pt = &adst->mesh.tria[k];
  p[0] = &adst->mesh.point[pt->v[0]];
  p[1] = &adst->mesh.point[pt->v[1]];
  p[2] = &adst->mesh.point[pt->v[2]];

  /* velocity at each vertex of k */
  u0 = &adst->sol.u[3*(pt->v[0]-1)+1];
  u1 = &adst->sol.u[3*(pt->v[1]-1)+1];
  u2 = &adst->sol.u[3*(pt->v[2]-1)+1];

  /* w = P1 velocity at the point with barycentric coordinates cb */
  w[0] = cb1[0]*u0[0] + cb1[1]*u1[0] + cb1[2]*u2[0];
  w[1] = cb1[0]*u0[1] + cb1[1]*u1[1] + cb1[2]*u2[1];
  w[2] = cb1[0]*u0[2] + cb1[1]*u1[2] + cb1[2]*u2[2];
  
  /* Transport w onto initial tangent plane to obtain v2 */
  if ( !nortri(&adst->mesh,k,n) ) return(0);
  if ( !invproj(n1,n,w,v2) ) return(0);
  
  /* Last integration point */
  k = *iel;
  memcpy(cb1,cb,3*sizeof(double));
  
  if ( !geod(adst,&k,cb1,step,v2) ) return(0);
  
  pt = &adst->mesh.tria[k];
  p[0] = &adst->mesh.point[pt->v[0]];
  p[1] = &adst->mesh.point[pt->v[1]];
  p[2] = &adst->mesh.point[pt->v[2]];

  /* velocity at each vertex of k */
  u0 = &adst->sol.u[3*(pt->v[0]-1)+1];
  u1 = &adst->sol.u[3*(pt->v[1]-1)+1];
  u2 = &adst->sol.u[3*(pt->v[2]-1)+1];

  /* w = P1 velocity at the point with barycentric coordinates cb */
  w[0] = cb1[0]*u0[0] + cb1[1]*u1[0] + cb1[2]*u2[0];
  w[1] = cb1[0]*u0[1] + cb1[1]*u1[1] + cb1[2]*u2[1];
  w[2] = cb1[0]*u0[2] + cb1[1]*u1[2] + cb1[2]*u2[2];
  
  /* Transport w onto initial tangent plane to obtain v3 */
  if ( !nortri(&adst->mesh,k,n) ) return(0);
  if ( !invproj(n1,n,w,v3) ) return(0);

  /* Sum transported velocity fields */
  v[0] =  ( v0[0] + 2.0 * (v1[0] + v2[0]) + v3[0] ) / 6.0;
  v[1] =  ( v0[1] + 2.0 * (v1[1] + v2[1]) + v3[1] ) / 6.0;
  v[2] =  ( v0[2] + 2.0 * (v1[2] + v2[2]) + v3[2] ) / 6.0;
  
  /* Backtrack geodesic with final velocity */
  if ( !geod(adst,iel,cb,step,v) ) return(0);

  return(1);
}

/* Solve advection */
int advec1_s(ADst *adst) {
  pTria      pt,pt1;
  pPoint     ppt,p0;
  double     dt,dte,norm,v0,v1,vb,tol,*u0,cb[3];
  int        nstep,nchar,ncor,nb,nc,j,k,kprv,l,iel,ip,ip0,ilist,list[AD_LONMAX];
  char       i;

  /* allocate structure for solution */
  if ( !adst->sol.new ) {
    if ( adst->info.zip)
      adst->sol.new = (double *)malloc((adst->info.npi+1)*sizeof(double));
    else
      adst->sol.new = (double *)malloc((adst->info.np+1)*sizeof(double));
    
    assert(adst->sol.new);
    memcpy(adst->sol.new,adst->sol.chi,(adst->info.np+1)*sizeof(double));
  }
  
  if ( !adst->info.nocfl ) {
    printf(" *** Unavailable option -- Abort.\n");
    exit(0);
  }
  
  /* Reset mark field for triangles and flag field for points */
  for (k=1; k<=adst->info.nt; k++)
    adst->mesh.tria[k].mark = 0;
  
  for (k=1; k<=adst->info.np; k++)
    adst->mesh.point[k].flag = 0;
    
  dt    = adst->sol.dt;
  nstep = 100;
  tol   = adst->sol.dt / nstep;
  nchar = 0;
  ncor  = 0;
  
  if ( adst->info.verb != '0' ) {
    fprintf(stdout,"    Time stepping: %g\n",dt);
    fprintf(stdout,"    Solving: "); fflush(stdout);
  }
  
  for (k=1; k<=adst->info.nt; k++) {
    pt = &adst->mesh.tria[k];
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      ppt = &adst->mesh.point[ip];
                  
      /* Check if already processed */
      if ( ppt->flag == 1 )  continue;
            
      /* Case of very small velocity: no motion */
      u0 = &adst->sol.u[3*(ip-1)+1];
      norm = sqrt( u0[0]*u0[0] + u0[1]*u0[1] + u0[2]*u0[2] );
      if ( norm < AD_EPSD ) {
        adst->sol.new[ip] = adst->sol.chi[ip];
        ppt->flag = 1;
        nchar++;
        continue;
      }
            
      /* Use travel to get to the foot of the characteristic */
      /* iel = k;
      memset(cb,0,3*sizeof(double));
      cb[i] = 1.0;
      dte = dt;
      ++adst->mesh.mark;
      
      while ( travel_s(adst,cb,&iel,&dte) );
      if ( iel == 0 ) return(0); */
      
      /* Go to the foot of the characteristic using RK4 scheme */
      iel = k;
      memset(cb,0,3*sizeof(double));
      cb[i] = 1.0;
      dte = dt;
      
      for (j=0; j<nstep; j++) {
        if ( !nxtptR_s(adst,&iel,cb,tol) ) break;
        dte -= tol;
      }

      if ( j < nstep ) continue;
      
      /* Interpolate value at the foot of characteristic */
      pt1 = &adst->mesh.tria[iel];
      adst->sol.new[ip] = cb[0]*adst->sol.chi[pt1->v[0]] \
                        + cb[1]*adst->sol.chi[pt1->v[1]] \
                        + cb[2]*adst->sol.chi[pt1->v[2]];
            
      /* Advection complete */
      if ( dte < AD_EPS ) {
        ppt->flag = 1;
        nchar++;
      }
      /* Extrapolation of characteristic curve when an external boundary is met */
      else if ( !adst->info.noex ) {
        pt1 = &adst->mesh.tria[iel];
        /* TODO !!*/
        /* Exit through an edge (one cb is 0) */
        
        
        /* Exit through a vertex (two cbs are 0) */

        // ppt->flag = 1;
        // nc++;
      }
    }
  }
    
  /* Post processing; interpolate sol.new at the (few) points where the previous procedure failed, according to the change in one of the neighbours */
  
  do {
    nc = 0;
    for (k=1; k<=adst->info.nt; k++) {
      pt = &adst->mesh.tria[k];
      for (i=0; i<3; i++) {
        ip = pt->v[i];
        ppt = &adst->mesh.point[ip];
                
        if ( ppt->flag ) continue;
        
        nb = 0;
        vb = 0.0;
        ilist = boulep_s(&adst->mesh,k,i,list);
        for (l=0; l<ilist; l++) {
          ip0 = list[l];
          p0 = &adst->mesh.point[ip0];
          if ( p0->flag ) {
            v0 = adst->sol.chi[ip0];
            v1 = adst->sol.new[ip0];
            vb += (v1-v0);
            nb++;
          }
        }
                
        if ( nb ) {
          adst->sol.new[ip] = adst->sol.chi[ip] + vb/nb;
          ppt->flag = 1;
          nc++;
        }
      }
    }
    ncor += nc;
  }
  while ( nc > 0 );
  
  /* Update points which have not yet received correction */
  for (k=1; k<=adst->info.np; k++) {
    ppt = &adst->mesh.point[k];
    if ( ppt->flag ) continue;
    adst->sol.new[k] = adst->sol.chi[k];
    ppt->flag = 1;
  }
  
  if ( adst->info.verb != '0' )
    fprintf(stdout,"%d characteristics, %d corrections\n",nchar,ncor);

  return(1);
}
