#include "advect.h"
#include "libmesh5.h"


int loadMesh(ADst *adst) {
  pPoint     ppt;
  pTetra     pt;
  pTria      pt1;
  double    *a,*b,dd;
  float      fp1,fp2,fp3;
  int        i,i1,k,inm;
  char      *ptr,data[128];
  static int edg[6][2] = {0,1, 0,2, 0,3, 1,2, 1,3, 2,3};

  strcpy(data,adst->mesh.name);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if ( !(inm = GmfOpenMesh(data,GmfRead,&adst->info.ver,&adst->info.dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if ( !(inm = GmfOpenMesh(data,GmfRead,&adst->info.ver,&adst->info.dim)) ) {
        fprintf(stderr," # %s file not found.\n",data);
        return(0);
      }
    }
  }
  else if ( !(inm = GmfOpenMesh(data,GmfRead,&adst->info.ver,&adst->info.dim)) ) {
    fprintf(stderr," # %s file not found.\n",data);
    return(0);
  }

  if ( adst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  adst->info.np = GmfStatKwd(inm,GmfVertices);
  adst->info.nt = GmfStatKwd(inm,GmfTriangles);
  adst->info.ne = GmfStatKwd(inm,GmfTetrahedra);

  if ( !adst->info.np ) {
    if ( adst->info.verb != '0' )  fprintf(stdout,"\n # missing data\n");
    return(0);
  }

  /* memory alloc */
  adst->mesh.point = (Point*)calloc(adst->info.np+1,sizeof(Point));
  assert(adst->mesh.point);
  if ( adst->info.nt > 0 ) {
    adst->mesh.tria  = (Tria*)calloc(adst->info.nt+1,sizeof(Tria));
    assert(adst->mesh.tria);
  }
  if ( adst->info.ne > 0 ) {
    adst->mesh.tetra  = (Tetra*)calloc(adst->info.ne+1,sizeof(Tetra));
    assert(adst->mesh.tetra);
  }

  /* read mesh vertices */
  if ( adst->info.dim == 2 ) {
    GmfGotoKwd(inm,GmfVertices);
    for (k=1; k<=adst->info.np; k++) {
      ppt = &adst->mesh.point[k];
      if ( adst->info.ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ppt->ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->ref);
    }
  }
  else {
    GmfGotoKwd(inm,GmfVertices);
    for (k=1; k<=adst->info.np; k++) {
      ppt = &adst->mesh.point[k];
      if ( adst->info.ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
    }
  } 

  /* read triangles */
  GmfGotoKwd(inm,GmfTriangles);
  adst->sol.hmin = AD_TGV;
  for (k=1; k<=adst->info.nt; k++) {
    pt1 = &adst->mesh.tria[k];
    GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
    if ( adst->info.dim == 2 ) {
      for (i=0; i<3; i++) {
        i1 = (i+1) % 3;
        a  = &adst->mesh.point[pt1->v[i]].c[0];
        b  = &adst->mesh.point[pt1->v[i1]].c[0];
        dd = sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]));
        dd = AD_MAX(AD_EPS,dd);
        adst->sol.hmin = AD_MIN(adst->sol.hmin,dd);
      }
    }
  }
  /* read tetrahedra */
  GmfGotoKwd(inm,GmfTetrahedra);
  for (k=1; k<=adst->info.ne; k++) {
    pt = &adst->mesh.tetra[k];
		pt->mark =0;
    GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&pt->ref);
    for (i=0; i<6; i++) {
      a  = &adst->mesh.point[pt->v[edg[i][0]]].c[0];
      b  = &adst->mesh.point[pt->v[edg[i][1]]].c[0];
      dd = sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) + (b[2]-a[2])*(b[2]-a[2]));
      adst->sol.hmin = AD_MIN(adst->sol.hmin,dd);
    }
  }
  GmfCloseMesh(inm);

  if ( adst->info.verb != '0' ) {
    fprintf(stdout," %d vertices",adst->info.np);
    if ( adst->info.nt )  fprintf(stdout,", %d triangles",adst->info.nt);
    if ( adst->info.ne )  fprintf(stdout,", %d tetrahedra",adst->info.ne);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* load solution (velocity) */
int loadSol(ADst *adst) {
  double       bufd[GmfMaxTyp],dd;
  float        buf[GmfMaxTyp];
  int          i,k,dim,ver,np,inm,type,size,offset,typtab[GmfMaxTyp];
  char        *ptr,data[128];

	if ( !adst->sol.namein )  return(-1);
  strcpy(data,adst->sol.namein);

  /* remove .mesh extension */
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';

  /* look for data file */
  ptr = strstr(data,".sol");
  if ( ptr ) {
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
  }
  else {
    /* first try to read binary file */
    strcat(data,".solb");
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    if ( !inm ) {
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    }
  }
  if ( !inm )  return(-1);

  if ( dim != adst->info.dim )  return(-1);
  np = GmfStatKwd(inm,GmfSolAtVertices,&type,&offset,&typtab);
  if ( !np || typtab[0] != GmfVec || np != adst->info.np )  return(-1);

  if ( adst->info.verb != '0' )  fprintf(stdout,"    %s :",data);

  /* read sol: assume velocity if 1st field */
  GmfGotoKwd(inm,GmfSolAtVertices);
  adst->sol.umax = 0.0;
  if ( ver == GmfFloat ) {
	  for (k=0; k<adst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,&buf);
      dd = 0.0;
      for (i=0; i<dim; i++) {
        adst->sol.u[dim*k+i+1] = buf[i];
        dd += buf[i]*buf[i];
      }
      dd = sqrt(dd);
      adst->sol.umax = AD_MAX(adst->sol.umax,dd);
    }
  }
  else {
	  for (k=0; k<adst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,&bufd);
      dd = 0.0;
      for (i=0; i<dim; i++) {
        adst->sol.u[dim*k+i+1] = bufd[i];
        dd += bufd[i]*bufd[i];
			}
      dd = sqrt(dd);
      adst->sol.umax = AD_MAX(adst->sol.umax,dd);
    }
  }
  GmfCloseMesh(inm);

  if ( adst->info.verb != '0' ) {
    fprintf(stdout," %d data vectors\n",adst->info.np);
  }

  return(1);
}

/* load characteristic function */
int loadChi(ADst *adst) {
  double       bufd[GmfMaxTyp];
  float        buf[GmfMaxTyp];
  int          k,inm,np,ver,dim,type,size,typtab[GmfMaxTyp];
  char        *ptr,data[128];

  if ( !adst->sol.namechi )  return(-1);
  strcpy(data,adst->sol.namechi);
  
  /* remove .mesh extension */
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';

  /* look for data file */
  ptr = strstr(data,".sol");
  if ( ptr ) {
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
  }
  else {
    /* first try to read binary file */
    strcat(data,".chi.solb");
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    if ( !inm ) {
      ptr  = strstr(data,".chi.solb");
      *ptr = '\0';
      strcat(data,".chi.sol");
      inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    }
  }
  if ( !inm )  return(-1);
  
  if ( dim != adst->info.dim )  return(-1);
  np = GmfStatKwd(inm,GmfSolAtVertices,&type,&size,typtab);
  if ( !np || typtab[0] != GmfSca || np != adst->info.np )  return(-1);

  if ( adst->info.verb != '0' )  fprintf(stdout,"    %s :",data);

  /* read characteristic function */
  GmfGotoKwd(inm,GmfSolAtVertices);
  if ( ver == GmfFloat ) {
    for (k=1; k<=adst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,buf);
      adst->sol.chi[k] = buf[0];
    }
	}
  else {
    for (k=1; k<=adst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,bufd);
      adst->sol.chi[k] = bufd[0];
    }
  }
  GmfCloseMesh(inm);

  if ( adst->info.verb != '0' ) {
    fprintf(stdout," %d data scalar\n",adst->info.np);
  }

  return(1);
}


int saveChi(ADst *adst) {
  double       dbuf[GmfMaxTyp];
  int          i,k,outm,type,typtab[GmfMaxTyp];
  char        *ptr,data[128];

  strcpy(data,adst->sol.nameout);  
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".new.chi.solb");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".new.chi.sol");
  }

  adst->info.ver = GmfDouble;
  if ( !(outm = GmfOpenMesh(data,GmfWrite,adst->info.ver,adst->info.dim)) ) {
    fprintf(stderr," # unable to open %s\n",data);
    return(0);
  }
  if ( adst->info.verb != '0' )  fprintf(stdout,"    %s:",data);
  type = 1;
  typtab[0] = GmfSca;

  /* write sol */
  GmfSetKwd(outm,GmfSolAtVertices,adst->info.np,type,typtab);
  for (k=1; k<=adst->info.np; k++) {
    dbuf[0] = adst->sol.new[k];
    GmfSetLin(outm,GmfSolAtVertices,dbuf);
  }
  GmfCloseMesh(outm);

  if ( adst->info.verb != '0' )  fprintf(stdout," %d data vectors\n",adst->info.np);

  return(1);
}
