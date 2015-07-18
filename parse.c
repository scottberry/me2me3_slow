#include "definitions.h"
/* 
   Functions for parsing command line arguments for gillespie
   algorithm simulations of chromatin-based epigenetics
   ============================================================
   Author: Scott Berry
   Institute: John Innes Centre
   ============================================================
*/

void usage(void)
{
  printf("Usage:\n");
  printf(" -c <control region>\n");
  printf(" -a <alpha>\n");
  printf(" -b <beta>\n");
  printf(" -g <G2 duration>\n");
  printf(" -t <firing threshold>\n");
  printf(" -i <identifier>\n");
  printf(" -p <PRC2 inhibition>\n");
  printf(" -s (seed based on identifier)\n");
  printf(" -r (DNA replication ON)\n");
  printf(" -m (start in K27me3 state)\n");
  printf(" -u (start in unmodified state)\n");
  exit (8);
}

void parseCommandLine(int argc, char *const *argv, chromatin *c, parameters *p) {
  int j;
  char buffer[256]="";

  /* set defaults */
  c->controlSites = c->sites;
  p->G2duration = 0.0; 
  p->alpha = 1.0;
  p->beta = 1.0; 
  p->firingThreshold = 1.0; 
  p->PRC2inhibition = 1.0;
  p->startM = FALSE;
  p->startU = FALSE;
  p->randomSeed = TRUE;
  p->seed = 0;
  strcpy(p->id,"\0");

  /* parse command line args */
  opterr = 0;
  while ((j = getopt (argc, argv, "c:a:b:i:smurg:p:t:")) != -1)
    switch (j)
      {
      case 'c':
        sprintf(buffer,"%s",optarg);
        c->controlSites = atoi(buffer);
        break;
        
      case 'a':
        sprintf(buffer,"%s",optarg);
        p->alpha = atof(buffer);
        break;

      case 'b':
        sprintf(buffer,"%s",optarg);
        p->beta = atof(buffer);
        break;

      case 'i':
        sprintf(p->id,"%s",optarg);
        p->seed = atoi(p->id);
        sprintf(p->id,"_%s",optarg);
        break;

      case 's':
        p->randomSeed = FALSE;
        break;
        
      case 'm':
        p->startM = TRUE;
        break;

      case 'u':
        p->startU = TRUE;
        break;

      case 'r':
        p->DNAreplication = TRUE;
        break;

      case 'g':
        sprintf(buffer,"%s",optarg);
        p->G2duration = atof(buffer);
        break;

      case 'p':
        sprintf(buffer,"%s",optarg);
        p->PRC2inhibition = atof(buffer);
        break;
        
      case 't':
        sprintf(buffer,"%s",optarg);
        p->firingThreshold = atof(buffer);
        break;
        
      default:
        usage();
      }
  return;
}
