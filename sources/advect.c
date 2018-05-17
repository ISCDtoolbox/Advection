/*
 * main program file for advect
 * (C) Copyright 2007- , ICS-SU
 *
 * This file is part of advect.
 *
 * advect is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * advect is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with advect.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "advect.h"


static void excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout,"  Abnormal stop\n");  break;
    case SIGBUS:
      fprintf(stdout,"  Code error...\n");  break;
    case SIGFPE:
      fprintf(stdout,"  Floating-point exception\n"); break;
    case SIGILL:
      fprintf(stdout,"  Illegal instruction\n"); break;
    case SIGSEGV:
      fprintf(stdout,"  Segmentation fault\n");  break;
    case SIGTERM:
    case SIGINT:
      fprintf(stdout,"  Program killed\n");  break;
  }
  exit(1);
}

static void usage(char *prog) {
  fprintf(stdout,"usage: %s [+/-v | -h] [-dt step] data_mesh[.mesh] [-c chi[.sol]] [-s data_vel[.sol]] [-o output[.sol]]\n",prog);
  fprintf(stdout,"\nOptions and flags:\n\
  --help       show the syntax and exit.\n\
  --version    show the version and date of release and exit.\n\n\
  -dt step     time step (time units)\n\
  -nocfl       avoid truncation of the advection time period. \n\
  -v           suppress any message (for use with function call).\n\
  +v           increase the verbosity level for output.\n\n\
  source.mesh    name of the mesh file\n\
  chi.sol        characteristic function (scalar)\n\
  data.sol       name of file containing the initial solution or boundary conditions\n\
  output.sol     name of the output file (displacement field)\n");
  exit(1);
}


/* parse command line */
static int parsar(int argc,char *argv[],ADst *adst) {
  int      i;
  char    *ptr,*data;
  
  i = 1;
  while ( i < argc ) {
    if ( (*argv[i] == '-') || (*argv[i] == '+') ) {
      switch(argv[i][1]) {
      case '-':
        if ( !strcmp(argv[i],"--help") )
          usage(argv[0]);
        else if ( !strcmp(argv[i],"--version") ) {
          fprintf(stdout,"%s: version: %s release: %s\n",argv[0],AD_VER,AD_REL);
          exit(1);
        }
        break;
      case 'h':
      case '?':
        usage(argv[0]);
        break;
      case 'c':
        if ( ++i < argc ) {
          adst->sol.namechi = argv[i];
          ptr = strstr(adst->sol.namechi,".sol");
          if ( !ptr )  strcat(adst->sol.namechi,".sol");
        }
        else {
          fprintf(stdout,"%s: missing parameter file\n", argv[0]);
          usage(argv[0]);
        }
        break;  
      case 'd':
        if ( !strcmp(argv[i],"-dt") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            adst->sol.dt = strtod(argv[i],NULL);
          else {
            fprintf(stderr,"%s: missing argument option\n",argv[0]);
            usage(argv[0]);
          }
        }
        else {
          fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
          usage(argv[0]); 
        }
        break;
      case 'i':
        if ( ++i < argc ) {
          adst->mesh.name = argv[i];
          ptr = strstr(adst->mesh.name,".mesh");
          if ( !ptr )  strcat(adst->mesh.name,".mesh");
        }
        else {
          fprintf(stdout,"%s: missing input file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nocfl") )
          adst->info.nocfl = 1;
        else if ( !strcmp(argv[i],"-noex") )
          adst->info.noex = 1;
        break;
      case 'o':
        if ( ++i < argc ) {
          adst->sol.nameout = argv[i];
          ptr = strstr(adst->sol.nameout,".sol");
          if ( !ptr )  strcat(adst->sol.nameout,".sol");
        }
        else {
          fprintf(stdout,"%s: missing data file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 's':
        if ( ++i < argc ) {
          adst->sol.namein = argv[i];
          ptr = strstr(adst->sol.namein,".sol");
          if ( !ptr )  strcat(adst->sol.namein,".sol");
        }
        else {
          fprintf(stdout,"%s: missing data file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'v':
        if ( !strcmp(argv[i],"-v") )
          adst->info.verb = '0';
        else if ( !strcmp(argv[i],"+v") )
          adst->info.verb = '+';
        else {
          fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
          usage(argv[0]);
        }
        break;
      default:
        fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    else {
      if ( adst->mesh.name == NULL ) {
        data = (char*)calloc(strlen(argv[i])+10,sizeof(char));
        strcpy(data,argv[i]);
        ptr = strstr(data,".mesh");
        if ( !ptr )  strcat(data,".mesh");
        adst->mesh.name = data;
      }
      else {
        fprintf(stdout,"%s: illegal option %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* check file names */
  if ( adst->mesh.name == NULL ) {
    if ( adst->info.verb != '0' )  fprintf(stderr,"%s: missing argument\n",argv[0]);
    usage(argv[0]);
  }

  return(1);
}


/* parsing boundary conditions */
static int parsdt(ADst *adst) {
  double   dt;
  FILE    *in;

  /* check for parameter file */
  in = fopen(".dt","r");
  if ( !in )  return(1);
  fscanf(in,"%lf",&dt);
  fclose(in);

  /* check value */
  if ( adst->sol.dt > 0.0 )
    adst->sol.dt = AD_MIN(adst->sol.dt,dt);

  return(1);
}


int main(int argc,char *argv[]) {
  ADst    adst;
  int     ier;
  char    stim[32];

  tminit(adst.info.ctim,TIMEMAX);
  chrono(ON,&adst.info.ctim[0]);

  /* interrupts */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  signal(SIGBUS,excfun);

  /* init structure */
  memset(&adst.mesh,0,sizeof(Mesh));
  memset(&adst,0,sizeof(Sol));
  adst.sol.dt = -1.0;

  /* global parameters */
  adst.info.dim  = 3;
  adst.info.ver  = 1;
  adst.info.verb = '1';
  adst.info.nocfl = 0;
  adst.info.noex = 0;
  
  /* parse command line */
  if ( !parsar(argc,argv,&adst) )  return(1);
  if ( adst.sol.dt < 0.0 && !parsdt(&adst) )  return(1);

  /* loading data */
	chrono(ON,&adst.info.ctim[1]);

  if ( adst.info.verb != '0' ) {
    fprintf(stdout," - ADVECT, Release %s, %s\n   %s\n\n",AD_VER,AD_REL,AD_CPY);
    fprintf(stdout," - LOADING DATA\n");
  }

  /* loading mesh */
	ier = loadMesh(&adst);
  if ( ier <= 0 )  return(1);

  /* allocating memory */
  if ( !adst.sol.u ) {
    adst.sol.u = (double*)calloc(adst.info.dim*adst.info.np+1,sizeof(double));
    assert(adst.sol.u);
  }

  /* loading velocity */
  if ( adst.sol.namein ) {
    ier = loadSol(&adst);
    if ( !ier )  return(1);
  }
  
  /* load characteristic function */
  adst.sol.chi = (double*)calloc(adst.info.dim*adst.info.np+1,sizeof(double));
  assert(adst.sol.chi);
  ier = loadChi(&adst);
  if ( ier <= 0 ) {
    if ( adst.info.verb != '0' )  fprintf(stdout," # missing or wrong file %s",adst.sol.namechi);
    return(1);
  }

  /* build adjacency table */
  ier = adst.info.dim == 2 ? hashel_2d(&adst) : hashel_3d(&adst);
  if ( !ier )  return(1);

  chrono(OFF,&adst.info.ctim[1]);
  printim(adst.info.ctim[1].gdif,stim);
  if ( adst.info.verb != '0' )  fprintf(stdout," - COMPLETED: %s\n",stim);

  /* solve advection */
  chrono(ON,&adst.info.ctim[2]);
  if ( adst.info.verb != '0' )
    fprintf(stdout,"\n ** MODULE ADVECT: %s\n",AD_VER);

  ier = AD_advect(&adst);
  if ( !ier )  return(1);

  chrono(OFF,&adst.info.ctim[2]);
  if ( adst.info.verb != '0' ) {
		printim(adst.info.ctim[2].gdif,stim);
    fprintf(stdout," ** COMPLETED: %s\n\n",stim);
	}

  /* save file */
  if ( adst.info.verb != '0' )  fprintf(stdout," - WRITING DATA\n");
  chrono(ON,&adst.info.ctim[3]);

  if ( !adst.sol.nameout ) {
    adst.sol.nameout = (char *)calloc(128,sizeof(char));
    assert(adst.sol.nameout);
    strcpy(adst.sol.nameout,adst.mesh.name);
  }

  ier = saveChi(&adst);
	if ( !ier )   return(1);
  chrono(OFF,&adst.info.ctim[3]);
  if ( adst.info.verb != '0' ) {
    printim(adst.info.ctim[3].gdif,stim);
    fprintf(stdout," - COMPLETED: %s\n",stim);
  }

  /* free mem */
	free(adst.sol.u);
  free(adst.sol.chi);

  chrono(OFF,&adst.info.ctim[0]);
  if ( adst.info.verb != '0' ) {
    printim(adst.info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s.\n",stim);
  }

  return(0);
}
