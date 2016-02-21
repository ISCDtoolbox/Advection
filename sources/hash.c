#include "advect.h"


int hashel_3d(ADst *adst) {
  return(1);
}


static int hcode_2d(Mesh *mesh,int a,int b,int k,int i) {
  Htab     *ht;
  Cell     *pc;
  pTria     pt,pt1;
  int       abmin,adj,sum;

  ht = &mesh->hash;
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
  pt  = &mesh->tria[k];
  do {
    pt1 = &mesh->tria[pc->elt];
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
      if ( !ht->hnxt )  return(0);
      return(1);
    }
    pc = &ht->cell[pc->nxt];
  } while (1);

  return(0);  
}


/* build adjacency table */
int hashel_2d(ADst *adst) {
  Htab    *ht;
  pTria    pt;
  pPoint   ppt;
  int      k,na;
  char     i,i1,i2;

  if ( adst->info.verb != '0' )  fprintf(stdout,"    Adjacency table: ");

  /* alloc hash */
  ht = &adst->mesh.hash;
  ht->nmax = (int)(3.71 * adst->info.np);
  ht->cell = (Cell*)calloc(ht->nmax+2,sizeof(Cell));
  assert(ht->cell);

  ht->hsiz = 2 * adst->info.np;
  ht->hnxt = ht->hsiz;
  for (k=ht->hsiz; k<ht->nmax; k++)
    ht->cell[k].nxt = k+1;

  /* update */
  na = 0;
  for (k=1; k<=adst->info.nt; k++) {
    pt = &adst->mesh.tria[k];
    for (i=0; i<3; i++) {
      i1 = (i+1) % 3;
      i2 = (i+2) % 3;
      if ( !hcode_2d(&adst->mesh,pt->v[i1],pt->v[i2],k,i) )  return(0);
      na++;
    }
  }

  /* add seed with point */
  for (k=1; k<=adst->info.nt; k++) {
    pt   = &adst->mesh.tria[k];
    for (i=0; i<3; i++) {
      if ( !pt->adj[i] )  adst->mesh.point[pt->v[1]].s = k;
    }
  }
  for (k=1; k<=adst->info.nt; k++) {
    pt = &adst->mesh.tria[k];
    for (i=0; i<3; i++) {
      ppt = &adst->mesh.point[pt->v[i]];
      if ( !ppt->s )  ppt->s = k; 
    }
  }

  if ( adst->info.verb != '0' )  fprintf(stdout," %d updated\n",na);

  return(1);  
}