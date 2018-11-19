/* PCF  */
/* May 13, 2014 */
/* Program to filter fasta files based on size criteria */


#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <sys/types.h>


#define CR 13
#define	LF	10
#define BOLD "\033[1m"
#define NORMAL "\033[0m"



int counter=-1;

int outputs;
int num[]={1,2,3,4};
float gcp;
void show_seq();
char csome[50];
float freq;
long aacode;
float total[100];
float tfreq;
int maxaacnt;
typedef unsigned long ulong;
char seqname[500];

char fatal_str[100];
int fatal_flag;

typedef struct st_dynamic_string
{
  char *str;
  uint length,max_length,alloc_increment;
} DYNAMIC_STRING;

typedef char my_bool;

extern my_bool init_dynamic_string(DYNAMIC_STRING *str, const char *init_str,uint init_alloc, uint alloc_increment);
extern my_bool dynstr_append(DYNAMIC_STRING *str, const char *append);
extern my_bool dynstr_append_mem(DYNAMIC_STRING *str, const char *append, uint length);
extern my_bool dynstr_set(DYNAMIC_STRING *str, const char *init_str);
extern my_bool dynstr_realloc(DYNAMIC_STRING *str, ulong additional_size);
extern void dynstr_free(DYNAMIC_STRING *str);


DYNAMIC_STRING d_query;
DYNAMIC_STRING *query;




typedef struct 
{
DYNAMIC_STRING d_item;
DYNAMIC_STRING *item;
char name[100];
int score;
double pvalue;
int flag;
} HITALIGN;

FILE *ofp;
int file_num=1;
char line[5000];
char sname[500];
char pname[500];
char oname[5000];
int debug=0;
int seq_per_file;
int tflag;
int pflag;
char router;


void useage()
{

	printf ("%s Useage is %s -f sequence file -n number of sequences per file -o output file name  ]%s\n",BOLD,pname,NORMAL);

}




void useage();
int mygetline(char *s, int n, FILE *f);
int main(argc,argv)
int argc;
char *argv[];
{
extern char *optarg;
extern int optind;
char fname[500];

int i,j;
int c,errflg;
char *dashQ, *dashM;
int queryN, matchN;
char *element; 
char selfmatch[200];
int rflag=0;
int srange_flag;
int qrange_flag;


int sum_count;
long sum_pos;
float sum_score;


int ident1,ident2;
double identpcnt,pospcnt;
int pos1,pos2;
int score;
double bits,expect,P;
int qstart, qend;
int sstart, send;
char strand1, strand2;
char *namesplit;
char qname[50];
int status;
int cnt;


int aacnt;
int scale=1000;
float myfreq;
long loc;
FILE *tfp;
float *array;
float sum;
float mean;
int high[12],low[12];
int excess=5;


seq_per_file=0;
sname[0]='\0';
errflg=3;

outputs=0;

	strcpy(pname,argv[0]);	   
        optind=1;
          while ((c = getopt(argc, argv, "Df:n:o:")) != -1)
               switch (c) {
	       case 'D':
	       		debug=1;
			break;
	       case 'f':  /* ref file */
                    strcpy(fname,optarg); errflg--;
                    if (debug) fprintf(stderr,"fasta file = %s\n", fname);
                    break;   		       						                 		    
 	       case 'n':
	       		seq_per_file=atoi(optarg);errflg--; 
			break;
 	       case 'o':
                    strcpy(oname,optarg); errflg--;
                    if (debug) fprintf(stderr,"fasta file = %s\n", fname);
                    break;  	
		    	
               default :
                    errflg++;
               }

          if (errflg)	
	  {
		useage();
                exit (2);
          }
	

       if (strcmp(fname,"-")==0) tfp = stdin;
        else if ((tfp=fopen(fname, "r"))==NULL)
                        {
                        printf ("ERROR: Unable to access file \"%s\" \n",fname);
                        exit (1);               
                        }
maxaacnt=0;
			    
cnt=0;	
aacnt=0;


					
	while (!(mygetline(line,5000,tfp)))
	{
	if (line[0] =='>')
		{
		if (cnt>0) show_seq();
		init_dynamic_string(&d_query, "", 16384, 65536);
		query = &d_query;
		cnt++;
		strncpy(seqname,line,499);
		}
	else	
		{
/*  		printf ("%s\n",line); */	
		dynstr_append_mem(query,line, strlen(line)); 
/* 		printf ("%s\n",d_query.str); */
		}
	} /*  end of while */	
	
	
show_seq();		

fprintf (stderr,"A total of %d sequences split into %d files with ~%d sequences per file\n",outputs,file_num-1,seq_per_file);			 
			 			
	if (ofp!=NULL) fclose (ofp);
			
} /*  end of MAIN */

/* now process the sequence in d_query*/		
void show_seq()
{			
int seqsize;
int i;
char outfile[500];


	
	if (counter==-1 || counter>=seq_per_file)
	{
	counter=0;
	if (ofp!=NULL) fclose (ofp);
	
	sprintf (outfile,"%s_%d.fa",oname,file_num++);
	
	if ((ofp=fopen(outfile, "w"))==NULL)
                        {
                        printf ("ERROR: Unable to open out file file \"%s\" \n",outfile);
                        exit (1);               
                        }	
	}



			
			
	seqsize=strlen(d_query.str);

	{	
	 if (debug) 	fprintf (stderr,"%s\t%d\t%d\t%d\n",seqname+1,seqsize,counter,file_num);
		fprintf (ofp,"%s\n",seqname);

			for (i=0; i<seqsize-70; i=i+70)
				{
 			fprintf (ofp,"%.70s\n",d_query.str+i);
				}
 			fprintf (ofp,"%.70s\n",d_query.str+i); 
	outputs++;
	counter++;
	}
	dynstr_free(query); 	
}




int mygetline(char *s, int n, FILE *f) {
    register int i=0;

    while(1) {
        s[i] = (char)fgetc(f);

        if(s[i] == CR)
            s[i] = fgetc(f);

        if((s[i] == 0x4) || (s[i] == LF) || (i == (n-1))) {
            s[i] = '\0';
            return (feof(f) ? 1 : 0);
        }
        ++i;
    }
}



