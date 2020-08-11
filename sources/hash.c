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
