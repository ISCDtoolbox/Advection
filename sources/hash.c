#include "advect.h"


/* data structures */
typedef struct {
  int    min,max,nxt,elt,ind;
} Cell;

typedef struct {
  Cell  *cell;
  int    nmax,hsiz,hnxt;
} Htab;


static int hcode_3d(Tetra *tetra,Htab *ht,int a,int b,int c,int k,int i) {
  Cell     *pc;
  pTetra    pt,pt1;
  int       adj,sum,min,max;

  sum = a+b+c;
  if ( sum >= ht->nmax )  return(0);

  /* check if edge ab stored */
  sum = sum % ht->hsiz;
  pc  = &ht->cell[sum];
  min = AD_MIN3(a,b,c);
  max = AD_MAX3(a,b,c);
  if ( !pc->min ) {
    pc->min = min;
    pc->max = max;
    pc->elt = k;
    pc->ind = i;
    return(1);
  }

  /* analyze linked list */
  pt  = &tetra[k];
  do {
    pt1 = &tetra[pc->elt];
    if ( pc->min == min && pc->max == max ) {
      adj = pt1->adj[pc->ind];
      if ( !adj ) {
        pt->adj[i]        = 4*pc->elt+pc->ind;
        pt1->adj[pc->ind] = 4*k+i;
      }
      return(1);
    }
    else if ( !pc->nxt ) {
      pc->nxt = ht->hnxt;
      pc      = &ht->cell[ht->hnxt];
      if ( !pc )  return(0);
      pc->min  = min;
      pc->max  = max;
      pc->elt  = k;
      pc->ind  = i;
      ht->hnxt = pc->nxt;
      pc->nxt  = 0;

      /* check for size overflow */
      if ( !ht->hnxt )  return(0);
      return(1);
    }
    pc = &ht->cell[pc->nxt];
  } while (1);

  return(0);
}


static int hcode_2d(Tria *tria,Htab *ht,int a,int b,int k,int i) {
  Cell     *pc;
  pTria     pt,pt1;
  int       abmin,adj,sum;

  sum = a+b;
  if ( sum >= ht->nmax )  return(0);

  /* check if edge ab stored */
  pc    = &ht->cell[sum];
  abmin = AD_MIN(a,b);
  if ( !pc->min ) {
    pc->min = abmin;
    pc->elt = k;
    pc->ind = i;
    return(1);
  }

  /* analyze linked list */
  pt  = &tria[k];
  do {
    pt1 = &tria[pc->elt];
    if ( pc->min == abmin ) {
      adj = pt1->adj[pc->ind];
      if ( !adj ) {
        pt->adj[i]        = 3*pc->elt+pc->ind;
        pt1->adj[pc->ind] = 3*k+i;
      }
      return(1);
    }
    else if ( !pc->nxt ) {
      pc->nxt = ht->hnxt;
      pc      = &ht->cell[ht->hnxt];
      if ( !pc )  return(0);
      pc->min  = abmin;
      pc->elt  = k;
      pc->ind  = i;
      ht->hnxt = pc->nxt;
      pc->nxt  = 0;

      /* check for size overflow */
      if ( !ht->hnxt ) return(0);
    
      return(1);
    }
    pc = &ht->cell[pc->nxt];
  } while (1);
  
  return(0);  
}

/* build adjacency table */
int hashel_3d(ADst *adst) {
  Htab     ht;
  pTetra   pt,pt1;
  pPoint   ppt;
  int      k,nt;
  char     i,i1,i2,i3;

  if ( adst->info.verb != '0' )  fprintf(stdout,"    Adjacency table: ");

  /* alloc hash */
  ht.nmax = (int)(12.71 * adst->info.np);
  ht.cell = (Cell*)calloc(ht.nmax+2,sizeof(Cell));
  assert(ht.cell);

  ht.hsiz = adst->info.np;
  ht.hnxt = ht.hsiz;
  for (k=ht.hsiz; k<ht.nmax; k++)
    ht.cell[k].nxt = k+1;

  /* update */
  nt = 0;
  for (k=1; k<=adst->info.ne; k++) {
    pt = &adst->mesh.tetra[k];
    for (i=0; i<4; i++) {
      i1 = (i+1) % 4;
      i2 = (i+2) % 4;
      i3 = (i+3) % 4;
      if ( !hcode_3d(adst->mesh.tetra,&ht,pt->v[i1],pt->v[i2],pt->v[i3],k,i) )  return(0);
      nt++;
    }
  }

  /* add seed with point */
  for (k=1; k<=adst->info.ne; k++) {
    pt   = &adst->mesh.tetra[k];
    for (i=0; i<4; i++) {
      if ( !pt->adj[i] )  adst->mesh.point[pt->v[i]].s = k;
    }
  }
  for (k=1; k<=adst->info.ne; k++) {
    pt = &adst->mesh.tetra[k];
    for (i=0; i<4; i++) {
      ppt = &adst->mesh.point[pt->v[i]];
      if ( !ppt->s )  ppt->s = k; 
    }
  }
  free(ht.cell);

  if ( adst->info.verb != '0' )  fprintf(stdout," %d updated\n",nt);

  return(1);
}


/* build adjacency table */
int hashel_2d(ADst *adst) {
  Htab     ht;
  pTria    pt;
  pPoint   ppt;
  int      k,na;
  char     i,i1,i2;

  if ( adst->info.verb != '0' )  fprintf(stdout,"    Adjacency table: ");

  /* alloc hash */
  ht.nmax = (int)(5.0 * adst->info.np);
  ht.cell = (Cell*)calloc(ht.nmax+2,sizeof(Cell));
  assert(ht.cell);

  ht.hsiz = 2 * adst->info.np;
  ht.hnxt = ht.hsiz;
  for (k=ht.hsiz; k<ht.nmax; k++)
    ht.cell[k].nxt = k+1;
  
  /* update */
  na = 0;
  for (k=1; k<=adst->info.nt; k++) {
    pt = &adst->mesh.tria[k];
    for (i=0; i<3; i++) {
      i1 = (i+1) % 3;
      i2 = (i+2) % 3;
      if ( !hcode_2d(adst->mesh.tria,&ht,pt->v[i1],pt->v[i2],k,i) ) return(0);
      na++;
    }
  }

  /* add seed with point */
  for (k=1; k<=adst->info.nt; k++) {
    pt   = &adst->mesh.tria[k];
    for (i=0; i<3; i++) {
      if ( !pt->adj[i] )  adst->mesh.point[pt->v[(i+1)%3]].s = k;
    }
  }
  for (k=1; k<=adst->info.nt; k++) {
    pt = &adst->mesh.tria[k];
    for (i=0; i<3; i++) {
      ppt = &adst->mesh.point[pt->v[i]];
      if ( !ppt->s )  ppt->s = k; 
    }
  }
  free(ht.cell);

  if ( adst->info.verb != '0' )  fprintf(stdout," %d updated\n",na);


  return(1);  
}

/* Extract points of a 3d mesh which are part of the surface triangulation,
   renumber them so that they are contiguous,
   and pack the values of solution and velocity accordingly */
int pack_s(ADst *adst,int *perm) {
  pTria        ptt;
  pPoint       p0,p1;
  int          k,npc;
  char         i;
  
  /* Reset flag field */
  for (k=1; k<=adst->info.np; k++)
    adst->mesh.point[k].flag = 0;
  
  /* Identify points of the surface triangulation */
  for (k=1; k<=adst->info.nt; k++) {
    ptt = &adst->mesh.tria[k];
    if ( !ptt->v[0] ) continue;
    for (i=0; i<3; i++){
      p0 = &adst->mesh.point[ptt->v[i]];
      p0->flag = 1;
    }
  }
  
  /* Compress points: p1->flag contains the original index */
  npc = 0;
  for (k=1; k<=adst->info.np; k++) {
    p0 = &adst->mesh.point[k];
    if ( p0->flag ) {
      npc++;
      p1 = &adst->mesh.point[npc];
      memcpy(p1,p0,sizeof(Point));
      p1->flag = k;
    }
  }
  
  adst->info.npi = adst->info.np;
  adst->info.np = npc;
  
  /* Compress solution */
  for (k=1; k<=adst->info.np; k++) {
    p0 = &adst->mesh.point[k];
    // printf("k %d et flag %d\n",k,p0->flag);
    adst->sol.chi[k] = adst->sol.chi[p0->flag];
    adst->sol.u[3*(k-1)+1] = adst->sol.chi[3*(p0->flag-1)+1];
    adst->sol.u[3*(k-1)+2] = adst->sol.chi[3*(p0->flag-1)+2];
    adst->sol.u[3*(k-1)+3] = adst->sol.chi[3*(p0->flag-1)+3];
  }
  
  /* Store permutation; perm[k] = original index of (new) point k */
  for (k=1; k<=adst->info.np; k++) {
    p0 = &adst->mesh.point[k];
    perm[k] = p0->flag;
  }
  
  return(1);
}

/* Restore values of sol at the position pointed by perm */
int unpack_s(ADst *adst,int *perm) {
  int k,ki;

  for (k=adst->info.np; k>=1; k--) {
    ki = perm[k];
    adst->sol.chi[ki] = adst->sol.chi[k];
    adst->sol.new[ki] = adst->sol.new[k];
    adst->sol.u[3*(ki-1)+1] = adst->sol.chi[3*(k-1)+1];
    adst->sol.u[3*(ki-1)+2] = adst->sol.chi[3*(k-1)+2];
    adst->sol.u[3*(ki-1)+3] = adst->sol.chi[3*(k-1)+3];
  }
  
  adst->info.np = adst->info.npi;
  return(1);
}
