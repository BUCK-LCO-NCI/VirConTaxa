/* This program takes a blastn output and added contig and taxonomy infop*/
/* Pcf 09/01/2016 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include <unistd.h>
//#include <mysql.h>

#define CR 13
#define	LF	10
#define BOLD "\033[1m"
#define NORMAL "\033[0m"
#define BUFSIZE 1000



char key[50][50]={"# scontig id", "contig length", "contig reads", "contig coverage","subject id","% identity","alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score", "subject title","tax id","Organism","Taxonomy"};



int key_len=19;

/*
typedef struct
{

int tax_id;
char query_id[100];
char subject_id[100];
double identity;
int align_len;
int mismatches;
int gap_opens;
int  q_start;
int  q_end;
int  s_start;
int  s_end;
double  evalue;
double bit_score;
char subject_title[200];
} FIELDS;

FIELDS field; */


char in_fields[25][8000];
typedef struct
{
char name[500];
int length;
int reads;
float coverage;
} ITEM;
int items;

ITEM *list_array;


ITEM target, *result;

void process_list();
int critter_cmp (const void *v1, const void *v2);
int find_critter ( const char *name,int count);
void print_critter (const ITEM *c);

char str[BUFSIZ]; /* Space to read string */
int found_item; /* Name to look for in table. */

int mygetline(char *s, int n, FILE *f);
int hashcount;
char outname[256];

int max_flag=0;
void no_match();
int unique;

char line[50000];

char *hselement;
char hsline[8000];
char *crelement;
char crline[8000];
char sname[500];
char pname[500];
char rname[500];
char cname[500];
char seqname[500];
int debug;
int b_cutoff;
int nflag;
int icount;
int segs;
int primerflag;
int plotflag;
char scale_flag;
int max;
char ofmt;
void useage()
{

	printf ("%s Useage %s %s\n",BOLD,pname,NORMAL);

}




void useage();
int mygetline(char *s, int n, FILE *f);
int main(argc,argv)
int argc;
char *argv[];
{
extern char *optarg;
extern int optind;
char *tline;
int i;
int c,errflg;
int count;
char *element;
int field=1;
char query[1000];
char *spacer;
int comment_flag=0;
int fld_cnt=0;

char last_id[250];
//    MYSQL *conn;
//    MYSQL_RES *res;
//    MYSQL_ROW row;
//    char *server = "biobase";
//    char *user = "pathuser";
//    char *password = "D@ta1reader"; /* set me first */
//    char *database = "pathodb";
//    conn = mysql_init(NULL);
//
//
// errflg=0;
//
//    /* Connect to database */
//    if (!mysql_real_connect(conn, server,
//          user, password, database, 0, NULL, 0)) {
//       fprintf(stderr, "%s\n", mysql_error(conn));
//       exit(1);
//    }

char *database = "/data/BUCKLAB/taxonomy/lineages-2018-06-13.reduce.tsv";

FILE *cfp;

 	strcpy(pname,argv[0]);
        optind=1;
          while ((c = getopt(argc, argv, "apndf:s:C:r:o:c:m")) != -1)
               switch (c) {
	       case 'd':
	       		debug=1;
			break;
	       case 'f':  /* input file */
                    strcpy(cname,optarg); errflg--;
                    if (debug) fprintf(stderr,"blasttabfile = %s\n", cname);
                    break;
	       case 'r':  /* read count file 3 columns */
                    strcpy(rname,optarg); errflg--;
                    if (debug) fprintf(stderr,"reference file  = %s\n", rname);
		    break;
	       case 'm':
	       	    max_flag=1;
		    break;
	       case 'o':
	       	    field=atoi(optarg);
		    break;
	       case 'c':   /* number of columns in the file */
	       	    comment_flag=1;;
		    break;

               default :
                    errflg++;
               }

/*            if (errflg)
	  {
		useage();
                exit (2);
          }  */
             /* Create table; no error checking is performed. */

last_id[0]='\0';


        if ((cfp=fopen(cname, "r"))==NULL)
	     {
	     printf ("ERROR: Unable to access cluster (.cdt) info file \"%s\" \n",cname);
	     exit (1);
	     }

	process_list();

	for (i=0;i<key_len;i++)
	{
		printf ("%s\t",key[i]);

	}
	printf ("\n");
count=0;
	while (!(mygetline(line,50000,cfp)))
	{
/* 	printf ("in the mygetline loop\n"); */

		tline=line;

	if (line[0]=='#')
	{
	if (comment_flag) printf ("%s\n",tline);
	}

		/* printf ("%d \n",count); */


/* 		printf ("%s - %s \n", element+3, tline); */

	else
	{
	fld_cnt=0;
	while ((element=(char*)strsep (&tline,"\t"))!=NULL)
		{

		strcpy(in_fields[fld_cnt++],element);
		}

if ((strcmp(last_id,in_fields[1])==0) && (max_flag ==1))
	{ /* do nothing */
	}
else {
		strcpy (last_id,in_fields[1]);

		printf ("%s\t",in_fields[1]);

	found_item=find_critter (in_fields[1],items);
//	fprintf( stderr, "%s\t%d\n", in_fields[1],items );
	if (found_item !=0)
		{
		printf ("%d\t%d\t%f\t",result->length,result->reads,result->coverage);
		}
	else { printf("0\t0\t0.0\t");}

		for (i=2;i<fld_cnt;i++)
		{
		printf ("%s\t",in_fields[i]);
		}

    char buf[BUFSIZE];
    FILE *fp;
    int num_rows=0;
		

        sprintf(query, "awk '$1 == %s {print}' %s", in_fields[0], database);
        fprintf(stderr, "%s\n", query);
        if ((fp = popen(query, "r")) == NULL) {
        	printf("Error opening pipe!\n");
        	exit(1);
    	} 
		while (fgets(buf, BUFSIZE, fp) != NULL) {
			num_rows++;
			fprintf(stderr, "%s", buf);
		}
        
	fprintf(stderr, "%d\n", num_rows);
    if(pclose(fp))  {
        printf("Command not found or exited with error status\n");
        exit(1);
    }

        	 if (num_rows < 1) /* no match */
				{
					printf("%s\t%s\tNO Taxonomy\n", tline,element);
				}
				else     /* got a hit */
				{
					printf("%s", buf);
				}


// 	sprintf (query, "SELECT * FROM `taxid` WHERE `taxon` = '%s';",in_fields[0]);
//         if (mysql_query(conn, query))
//           {
// 	    fprintf(stderr, "ERROR: %s\n", mysql_error(conn));
//             exit(1);
//    	  }
//
//          res = mysql_store_result(conn);
// 	 num_rows=mysql_num_rows(res);

//  	 while ((row = mysql_fetch_row(res)) != NULL)
// 	 	{
//  		while ((spacer=strchr(row[4],';')) != NULL) *spacer='\t';
//    	   	printf("%s\t%s\t%s\n",row[1],row[2],row[3]);
// 	   	}
//	 }

         }
			 }

	  count++;
	}


	     if (cfp!=NULL) fclose(cfp);
    /* close connection */
//    mysql_free_result(res);
//    mysql_close(conn);
	exit(0);

} /* end of main */



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



void process_list()
{
int count;
char *element;
char *tline;
char *spacer;
FILE *rfp;

                   if ((rfp=fopen(rname, "r"))==NULL)
	     {
	     printf ("ERROR: Unable to access name list file \"%s\" \n",rname);
	     exit (1);
	     }
/* read the number of lines in the file */
	while (!(mygetline(line,500,rfp)))
	{
		items++;
	}
rewind(rfp);


  	list_array = (ITEM *) calloc (items, sizeof (ITEM));
count=0;
	while (!(mygetline(line,500,rfp)))
	{
		tline=line;
// printf ("%s\n",line);

		element=(char*)strsep (&tline,"\t");
		strcpy(list_array[count].name,element);
		element=(char*)strsep (&tline,"\t");
		list_array[count].length=atoi(element);
		element=(char*)strsep (&tline,"\t");
		list_array[count].reads=atoi(element);
		element=(char*)strsep (&tline,"\t");
		list_array[count].coverage=atol(element);
//    fprintf (stderr,"%s\t%d\t%d\t%d\t%f\n", list_array[count].name,list_array[count].length,list_array[count].reads,list_array[count].coverage);
 	  count++;

	}

  qsort (list_array, count, sizeof(ITEM), critter_cmp);
/* 	     if (rfp!=NULL) fclose(rfp);  */


/* hash table built now move on */

fprintf (stderr,"readin %d & %d items.\n",items, count);
fflush(NULL);

}

int
critter_cmp (const void *v1, const void *v2)
{
  const ITEM  *c1 = v1;
  const ITEM  *c2 = v2;

/* 	fprintf (stderr,"%s\t%s\n",c1->name,c2->name); */

  return strcmp (c1->name, c2->name);
}


/* Print information about a critter. */

void
print_critter (const ITEM *c)
{
  printf ("%s\n", c->name);
}


/* Do the lookup into the sorted array. */

int find_critter ( const char *name,int count)
{

  strcpy (target.name,name);
  result = bsearch (&target, list_array, count, sizeof(ITEM),critter_cmp);
  if (result)
  {
/*     print_critter (result); */
    return (1);
   }
  else
    {
/*     printf ("Couldn't find %s.\n", name); */
    return(0);
    }
}
