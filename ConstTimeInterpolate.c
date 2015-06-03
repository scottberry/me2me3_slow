#include "definitions.h"

void usage(void)
{
  printf("Usage:\n");
  printf(" -t<timefile>\n");
  printf(" -d<datafile>\n");
  printf(" -s<timestep>\n");
  exit (8);
}

/* Replace a character of a string */
char *str_replace(char *orig, char *rep, char *with) {
    char *result; // the return string
    char *ins;    // the next insert point
    char *tmp;    // varies
    int len_rep;  // length of rep
    int len_with; // length of with
    int len_front; // distance between rep and end of last rep
    int count;    // number of replacements

    if (!orig)
        return NULL;
    if (!rep)
        rep = "";
    len_rep = strlen(rep);
    if (!with)
        with = "";
    len_with = strlen(with);

    ins = orig;
    for (count = 0; (tmp = strstr(ins, rep)); ++count) {
        ins = tmp + len_rep;
    }

    // first time through the loop, all the variable are set correctly
    // from here on,
    //    tmp points to the end of the result string
    //    ins points to the next occurrence of rep in orig
    //    orig points to the remainder of orig after "end of rep"
    tmp = result = malloc(strlen(orig) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    while (count--) {
        ins = strstr(orig, rep);
        len_front = ins - orig;
        tmp = strncpy(tmp, orig, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        orig += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, orig);
    return result;
}

int main(int argc, char *argv[]) {
  FILE *tFilePtr, *dFilePtr, *outFilePtr;
  char tFileName[128]="", dFileName[128]="", buffer[128]="";
  char outFileName[128]="", txt[16]=".txt", int_txt[16]="_int.txt";
  long step = 60, j;
  logical standardOut = FALSE;
  double t_in, level, levelold, t_out;

  /* Parse command line */
  opterr = 0;
  while ((j = getopt (argc, argv, "t:d:s:o")) != -1)
    switch (j)
      {
      case 't':
        sprintf(tFileName,"%s",optarg);
        break;
        
      case 'd':
        sprintf(dFileName,"%s",optarg);
        break;

      case 's':
        sprintf(buffer,"%s",optarg);
        step = atoi(buffer);
        break;

      case 'o':
        standardOut = TRUE;
        break;
        
      default:
        usage();
      }

  if (!strcmp(tFileName,"") || !strcmp(dFileName,"")) {
    printf("Mandatory argument(s) missing\n");
    exit(1);
  }
  
  tFilePtr = fopen(tFileName,"r");
  dFilePtr = fopen(dFileName,"r");

  sprintf(outFileName,"%s",str_replace(dFileName,txt,int_txt));
  outFilePtr = fopen(outFileName,"w");
  if (standardOut == TRUE)
    outFilePtr = stdout;
  
  t_out = 0;
  fscanf(tFilePtr, "%lf", &t_in);
  fscanf(dFilePtr, "%lf", &level);
  fprintf(outFilePtr, "%0.4f %0.4f\n", t_out, level);

  while(!feof(tFilePtr)) {
    levelold = level; 
    fscanf(tFilePtr, "%lf", &t_in);
    fscanf(dFilePtr, "%lf", &level);
    while (t_out < t_in) {
      fprintf(outFilePtr, "%0.4f %0.4f\n", t_out, levelold);
      t_out += step;
    }
  }

  fclose(tFilePtr);
  fclose(dFilePtr);
  fclose(outFilePtr);
  
  return(1);
}
