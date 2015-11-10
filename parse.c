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
  printf(" -s (seed based on identifier)\n");
  printf(" -r (DNA replication ON)\n");
  printf(" -m (start in K27me3 state)\n");
  printf(" -u (start in unmodified state)\n");
  printf(" -n <translation efficiency (noise)>\n");
  printf(" -h <histone turnover rate (per histone per transcription)>\n");
  exit (8);
}

void parseCommandLine(int argc, char *const *argv, chromatin *c, parameters *p) {
  int j;
  char buffer[256]="";

  /* set defaults */
  c->controlSites = c->sites;
  p->alpha = 1.0;
  p->beta = 1.0; 
  p->firingThreshold = 1.0; 
  p->startM = FALSE;
  p->startU = FALSE;
  p->randomSeed = TRUE;
  p->seed = 0;
  p->stochasticTranslationEfficiency = 1;
  p->transcription_turnover = 0.004;
  strcpy(p->id,"\0");
  sprintf(p->executable,"%s",argv[0]);
  
  /* parse command line args */
  opterr = 0;
  while ((j = getopt (argc, argv, "c:a:b:i:smurt:n:h:")) != -1)
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
        
      case 't':
        sprintf(buffer,"%s",optarg);
        p->firingThreshold = atof(buffer);
        break;

      case 'n':
        sprintf(buffer,"%s",optarg);
        p->stochasticTranslationEfficiency = atoi(buffer);
        break;

      case 'h':
        sprintf(buffer,"%s",optarg);
        p->transcription_turnover = atof(buffer);
        break;
              
      default:
        usage();
      }
  return;
}
