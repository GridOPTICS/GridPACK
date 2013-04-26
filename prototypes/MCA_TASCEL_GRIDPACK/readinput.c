
/*
 * C routine *
 * load entire file into buffer, broadcast buffer
 * then read buffer at receiver processor

 Yousu Chen; Jarek Nieplocha 
 7/6/2007
 */

 
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <mpi.h>
#include "elio.h"



void my_error(int me, char* str, long code)
{
fprintf(stderr,"%d %s %ld\n",me, str,code);
MPI_Abort(MPI_COMM_WORLD,2); 
}


char *nextline(char*buf)
{
char *tmp;

tmp = strchr(buf,'\n');
if(!tmp) my_error(0,"error moving to next line\n",0);

return tmp+1;

}
/*
int binary_search (int array[10000], int n, int key) 
{
  int min = 0;
  int max = n - 1;

  while (max >= min) 
    {
      int i = (min + max) / 2;
      if (key < array[i]) 
        max = i - 1;
      else if (key > array[i]) 
        min = i + 1;
      else 
        return i;
    }

  return -1;
}
*/

void getinput(char *buf, int *slack, int *refgen, double *baseMVA,
              int *nb, int *ng, int *nbrch, int *ns, int *nxfmr_adj,
              int *narea, int *ndc, int *nshunt,
              // bus array
              int *bus_i, int *bus_type, double *pd, double *qd, 
	      double *gs, double *bs, int *bus_area, double *vm, 
	      double *va, char *busname, double *base_kv, 
	      char *zone,
              // generator arrays
	      int *gen_bus, char *gen_id, double *pg, double *qg, 
	      double *qmax, double *qmin, double *vg, double *mbase,
              int *gen_status, double *pmax, double *pmin, 
              //  branch arrays 
	      int *f_bus, int *t_bus, char *br_id, double *br_r,
	      double *br_x, double *br_b, double *rate_a, double *rate_b,
	      double *rate_c, double *tap, double *shift, int *br_status,
              //  shunt arrays 
	      int *shunt_bus, int *modsw, double *vswhi, 
	      double *vswlo, int *swrem, double *binit, int *n1, double *b1,
              int *n2, double *b2, int *n3, double *b3, int *n4, double *b4,
              int *n5, double *b5, int *n6, double *b6, int *n7, double *b7,
              int *n8, double *b8)
{
int s;
int i;
char *ptr=buf, *tmp;
int bu1,bu2,bu7;
double bu3,bu4,bu5,bu6,bu8,bu9,bu11;
char   bu10[20],bu12[20];


int g1,g8,g11;
char g2[20];
char dum[20];
double g3,g4,g5,g6,g7,g9,g10,g12,g13,g14;

int br1,br2,br13;
char br3[20];
double br4,br5,br6,br7,br8,br9,br10,br11,br12;

int s1,s2,s5,s7;
double s3,s4,s6,s8;

double totload, totgen;

char *temps; char *tempt;
char *tmp2;
char *p;
int length;
int isolated[10000];
int isolated_cnt;
int exist_flag1, exist_flag2;

*slack=0;
*refgen=0;


totload =0.0;
totgen = 0.0;

isolated_cnt=0;
bzero(bu10,14);

//printf(" Start reading buffer\n");
sscanf(ptr,"%d %lf",&s,&(*baseMVA));


fflush(stdout);
ptr=nextline(ptr);
ptr=nextline(ptr);
ptr=nextline(ptr);

//printf(" reading bus data...\n");
for (i=0;i<*nb;i++) {
     sscanf(ptr,"%d,%d,%lf,%lf,%lf,%lf,%d,%lf,%lf,%10c,%lf,%s,",
            &bu1,&bu2,&bu3,&bu4,&bu5,&bu6,&bu7,&bu8,&bu9,bu10,&bu11,bu12);
	    
// ignore isolate buses
//     if (b2 <4) {
     bus_i[i]=bu1;
     bus_type[i]=bu2;
     pd[i]=bu3; 
     qd[i]=bu4; 
     gs[i]=bu5/(*baseMVA);
     bs[i]=bu6/(*baseMVA);
//     if ((b5 >=1.0e-16) || (b6 >=1.0e-16)) (*ns)= (*ns) + 1;// s = 1;s ++; // = *ns + 1;    
     if (((double)bu5 !=(double)0.) || ((double)bu6 !=(double)0.0)) (*ns)= (*ns) + 1;// s = 1;s ++; // = *ns + 1;    
     bus_area[i]=bu7;
     vm[i]=bu8;
     va[i]=bu9;
     totload = totload + bu3;

     if (bu2==3) *slack = i;
     base_kv[i]=bu11;
  //   } else {
 //    isolated[isolated_cnt]=b1;
//     isolated_cnt ++;

//     }
     ptr=nextline(ptr);
	    
     fflush(stdout);

}
//printf(" reading bus data done!\n");
     ptr=nextline(ptr);
//printf(" reading gen data...\n");
   for (i=0; i<*ng;i++)
   {
        char temp[5];
        sscanf(ptr, "%d,%4c,%lf,%lf,%lf,%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%d,%lf,%lf,%lf",
	&g1,temp,&g3,&g4,&g5,&g6,&g7,&g8,&g9,&g10,&g10,&g10,&g10,&g10,&g11,&g12,&g13,&g14);
//        if (g11>0) {
//        for (j=0;j<isolated_cnt;j++) 
//        {  
//        exist_flag1=binary_search (isolated,g1,g1);
//        if (exist_flag1 != -1) {
        gen_bus[i]=g1;
        gen_id[0]=temp[1]; /* this strips the opening quote */
        gen_id[1]=temp[2];
	gen_id +=2;        /* this is a fortran array of 2-char strings */
        pg[i]=g3;
        qg[i]=g4;
        qmax[i]=g5;
        qmin[i]=g6;
        vg[i]=g7;
        mbase[i]=g9;
        gen_status[i]=g11;
        pmax[i]=g13;
        pmin[i]=g14;
        ptr=nextline(ptr); /* move to the next line */
	totgen = totgen + g3;
//	}
	if (g1 == *slack) *refgen = i;
  //      printf("%d,%lf,%lf,%lf,%lf,%lf,%d,%lf,%lf,%d,%lf,%lf,%lf,\n",
//	gen_bus[i],pg[i],qg[i],qmax[i],qmin[i],vg[i],g8,mbase[i],g10,gen_status[i],g12,pmax[i],pmin[i]);
 //       printf("OK\n");
   }
   if (totgen <= totload)
   {
        pmax[*refgen] = pg[*refgen] + 1.1 * totload - totgen;
//	printf("*** Insufficient generation, setting Pmax at slack bus to %f\n", pmax[*refgen]);
   }
     ptr=nextline(ptr);

//printf(" reading gen data done\n");
//printf(" reading branch data...\n");
   /* start reading branch arrays */
   for (i=0; i<*nbrch;i++)
   {
        char temp[5];
        char s = '\'';
        sscanf(ptr, "%d,%d,%4c,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
	&br1,&br2,temp,&br4,&br5,&br6,&br7,&br8,&br9,&br10,&br11,&br12,&br12,&br12,&br12,&br13);
//        if (br13 >0) {
//        exist_flag1=binary_search (isolated,br1,br1);
 //       exist_flag2=binary_search (isolated,br2,br2);
  //      if ((exist_flag1 != -1) && (exist_flag2 != -1)) {
        
        if (temp[0]==s)
        {
          br_id[0]=temp[1]; 
          br_id[1]=temp[2];
        } else {
          sscanf(ptr, "%d,%d,%2c,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
	  &br1,&br2,temp,&br4,&br5,&br6,&br7,&br8,&br9,&br10,&br11,&br12,&br12,&br12,&br12,&br13);
          br_id[0]=temp[0]; 
          br_id[1]=temp[1];
        }
	br_id +=2;        /* this is a fortran array of 2-char strings */
        f_bus[i]=abs(br1);
        t_bus[i]=abs(br2);
        br_r[i]=br4;
        br_x[i]=br5;
        br_b[i]=br6;
//	if (br7<=1.0e-3) br7 = 9900.0; 
	if ((double)br7 == (double)0.0) br7 = 9900.0; 
        rate_a[i]=br7;
        rate_b[i]=br8;
        rate_c[i]=br9;
        tap[i]=br10;
        shift[i]=br11;
        br_status[i]=br13;
//        }
        ptr=nextline(ptr);/* move to the next line */
//        printf = ("%s\n", br_id[i]);
  //      printf(" %d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d\n",
//	f_bus[i],t_bus[i],br_r[i],br_x[i],br_b[i],rate_a[i],rate_b[i],rate_c[i],tap[i],shift[i],br_status[i]);
 }
//printf(" reading branch data done!\n");
 ptr=nextline(ptr);
   /* start reading transformer adjustment data, ignore so far */
 for (i=0; i<*nxfmr_adj;i++)
 {   
        ptr=nextline(ptr);
 }
 ptr=nextline(ptr);
   /* start reading area data, ignore so far */
 for (i=0; i<*narea;i++)
 {   
        ptr=nextline(ptr);
 }
 ptr=nextline(ptr);
   /* start reading dc data, ignore so far */
 for (i=0; i<*ndc;i++)
 {   
        ptr=nextline(ptr);
 }
 ptr=nextline(ptr);
   /* start reading shunt data */
 for (i=0; i<*nshunt;i++)
 {   
        sscanf(ptr, "%d,%d,%lf,%lf,%d,%lf,%d,%lf,",
	&s1,&s2,&s3,&s4,&s5,&s6,&s7,&s8);
        ptr=nextline(ptr);
        printf(" %d,%d,%lf,%lf,%d,%lf,%d,%lf\n",
	s1,s2,s3,s4,s5,s6,s7,s8);
 }


}

#define ROOT 0
//void test_file_bcast(char* fname)
void readinput_(char *fname, int *slack, int *refgen, double *baseMVA, 
              int *nb, int *ng, int *nbrch, int *ns, int *nxfmr_adj,
              int *narea, int *ndc, int *nshunt,
              // bus array
              int *bus_i, int *bus_type, double *pd, double *qd, 
	      double *gs, double *bs, int *bus_area, double *vm, 
	      double *va, char  *busname, double *base_kv, 
	      char *zone,
              // generator arrays
	      int *gen_bus, char *gen_id, double *pg, double *qg, 
	      double *qmax, double *qmin, double *vg, double *mbase,
              int *gen_status, double *pmax, double *pmin, 
              //  branch arrays 
	      int *f_bus, int *t_bus, char *br_id, double *br_r,
	      double *br_x, double *br_b, double *rate_a, double *rate_b,
	      double *rate_c, double *tap, double *shift, int *br_status,
              //  shunt arrays 
	      int *shunt_bus, int *modsw, double *vswhi, 
	      double *vswlo, int *swrem, double *binit, int *n1, double *b1,
              int *n2, double *b2, int *n3, double *b3, int *n4, double *b4,
              int *n5, double *b5, int *n6, double *b6, int *n7, double *b7,
              int *n8, double *b8)
{
  Fd_t fd;
  Off_t size;
  long lsize;
  int me;
  char *fbuffer;
  int p;
  double t0,t1;

   MPI_Comm_rank(MPI_COMM_WORLD, &me);
   if(me==ROOT){
     fd = elio_open(fname, ELIO_R, ELIO_PRIVATE);
     if(!fd)  my_error(me, "no such input file", 0);
//     printf ("%s\n",fname);
     if(elio_length(fd, &size) != ELIO_OK) my_error(me,"error length",0); 

     lsize = (long)size;

//     printf("file size is %ld\n",(long)size);
//     printf("broadcasting buffer size\n");
   }
   t0=MPI_Wtime();
   MPI_Bcast(&lsize,1,MPI_LONG, ROOT, MPI_COMM_WORLD); 
   t1=MPI_Wtime();
   if(me==ROOT){
     //printf("broadcasting buffer size done! It elapsed %8.4f seconds\n", t1-t0);
   }  

   if(lsize<1 || lsize > INT_MAX)  my_error(me," bad file size",lsize);

   fbuffer = (char *) malloc((lsize+1)*sizeof(char));
   if(!fbuffer) my_error(me,"could not allocate buffer",lsize);
   bzero(fbuffer, sizeof(fbuffer));
   if(me == ROOT){
      
     // printf("reading file into buffer...\n");
      size = elio_read(fd,0,fbuffer,(Size_t)lsize);
     // printf("reading file into buffer done!\n");
      if((long)size != lsize) my_error(me, "wrong size read", (long)size);

      elio_close(fd);
      //printf("broadcasting buffer...\n");
   }

   MPI_Bcast(fbuffer,(int)lsize, MPI_BYTE, ROOT, MPI_COMM_WORLD);
   if(me==ROOT){
     //printf("broadcasting buffer done!\n");
     //printf("reading buffer...\n");
   }  
 //  printf("%s",fbuffer);
/*   printf("%d bytes=%ld got: %d.*s\n",me,lsize,8,fbuffer);*/
      
   getinput(fbuffer, slack, refgen, baseMVA, 
              nb, ng, nbrch, ns,
              nxfmr_adj, narea, ndc, nshunt,

              // bus array
              bus_i, bus_type, pd, qd, 
	      gs, bs, bus_area, vm, 
	      va, busname, base_kv, 
	      zone,
              // generator arrays
	      gen_bus, gen_id, pg, qg, 
	      qmax, qmin, vg, mbase,
              gen_status, pmax, pmin, 
              //  branch arrays 
	      f_bus, t_bus, br_id, br_r,
	      br_x, br_b, rate_a, rate_b,
	      rate_c, tap, shift, br_status,
              //  shunt arrays 
              shunt_bus, modsw, vswhi, 
	      vswlo, swrem,binit, n1, b1,
              n2, b2, n3, b3, n4, b4,
              n5, b5, n6, b6, n7, b7, n8,b8);

   if(me==ROOT){
     printf("reading buffer done!\n");
     fflush(stdout);
   }  
   free(fbuffer);
   
}

void getca(char *fbuff, char *catype, int *ca1, int *ca2, char *caid, int *nca)
{
char *ptr=fbuff, *tmp;
char *temps; char *tempt;
char *tmp1;
char *tmp2;
char c1[6],c4[6],c[5];
int  i, i2,i3;


bzero(c1,6);
bzero(c4,6);
bzero(c,5);


   ptr=nextline(ptr);
   for (i=0; i<*nca;i++)
   {    
        tmp = ptr;
	sscanf(ptr,"%4c",c);
//	printf("type = %s\n",c);
    	ptr=tmp;
	if (strcmp(c,"brch") ==0)
        { 
           sscanf(ptr,"%4c,%2c,%d,%d",catype,caid,&i2,&i3);
	   ca1[i]=i2;
	   ca2[i]=i3;
	} else if(strcmp(c,"genr") ==0) {
           sscanf(ptr,"%4c,%2c,%d",catype,caid,&i2);
	   ca1[i]=i2;
	}   
	caid +=2;        /* this is a fortran array of 2-char elements */
	catype +=4;        /* this is a fortran array of 4-char elements */
        ptr=nextline(ptr);
   }	    
}





void  readcalist_(char *fname, char *catype, int *ca1, int *ca2, char *caid, int *nca) 
{
  Fd_t fd;
  Off_t size;
  long lsize;
  int me;
  char *fbuffer;
  int p;

   MPI_Comm_rank(MPI_COMM_WORLD, &me);
   if(me==ROOT){
     fd = elio_open(fname, ELIO_R, ELIO_PRIVATE);
      if(!fd)  my_error(me, "no such calist input file", 0);
      printf ("%s\n",fname);
      if(elio_length(fd, &size) != ELIO_OK) my_error(me,"error length",0); 

      lsize = (long)size;

      printf("file size is %ld\n",(long)size);
   }
    MPI_Bcast(&lsize,1,MPI_LONG, ROOT, MPI_COMM_WORLD);
    if(lsize<1 || lsize > INT_MAX)  my_error(me," bad file size",lsize);

   // fbuffer = malloc(lsize+1);
    fbuffer = (char *) malloc((lsize+1)*sizeof(char));
    if(!fbuffer) my_error(me,"could not allocate buffer",lsize);
    bzero(fbuffer, sizeof(fbuffer));

    if(me == ROOT) {
      size = elio_read(fd,0,fbuffer,(Size_t)lsize);
      if((long)size != lsize) my_error(me, "wrong size read", (long)size);
      elio_close(fd);
    }
      MPI_Bcast(fbuffer,(int)lsize, MPI_BYTE, ROOT, MPI_COMM_WORLD);
      getca(fbuffer, catype, ca1, ca2, caid, nca);
   free(fbuffer);
}

void getconfig(char *fbuff, char *fname2,int *M, int *NC, int *caopt, int *outflag, int *ncase, int *icase, char *arg_dir, char *cafile)//)//)//)//)//)//)//)//)//, int *max_it, double *tol)
{
   char *ptr=fbuff;
   
   ptr=nextline(ptr);
   sscanf(ptr,"%s,",fname2);
   printf("fname2=%s\n", fname2);
   ptr=nextline(ptr);
   ptr=nextline(ptr);
   sscanf(ptr,"%d,",M);
   printf("M=%d\n", *M);
   ptr=nextline(ptr);
   ptr=nextline(ptr);
   sscanf(ptr,"%d,",NC);
   printf("NC=%d\n", *NC);
   ptr=nextline(ptr);
   ptr=nextline(ptr);
   sscanf(ptr,"%d,",caopt);
   printf("CAOPT=%d\n", *caopt);
   ptr=nextline(ptr);
   ptr=nextline(ptr);
   sscanf(ptr,"%d,",outflag);
   printf("OUTFLAG=%d\n", *outflag);
   ptr=nextline(ptr);
   ptr=nextline(ptr);
   sscanf(ptr,"%d,",ncase);
   printf("ncase=%d\n", *ncase);
   ptr=nextline(ptr);
   ptr=nextline(ptr);
   sscanf(ptr,"%d,",icase);
   printf("icase=%d\n", *icase);
//   ptr=nextline(ptr);
//   ptr=nextline(ptr);
//   sscanf(ptr,"%lf,",tol);
//   printf("tol=%lf\n", *tol);
//   ptr=nextline(ptr);
//   ptr=nextline(ptr);
//   sscanf(ptr,"%d,",max_it);
//   printf("max_it=%d\n", *max_it);
   ptr=nextline(ptr);
   ptr=nextline(ptr);
   sscanf(ptr,"%s,",cafile);
   printf("cafile=%s\n", cafile);
   ptr=nextline(ptr);
   ptr=nextline(ptr);
   sscanf(ptr,"%s,",arg_dir);
   printf("arg_dir=%s\n", arg_dir);
//   printf("M=%d, NC=%d, caopt=%d, ncase1=%d, ncase2=%d\n", *M,*NC,*caopt,*ncase1,*ncase2);
//   printf("bin=%d, arg_dir=%s, cafile=%s\n", *bin,arg_dir,cafile);
}
      
      
void  readconfig_(char *fname, char *fname2, int *M, int *NC, int *caopt, int *outflag, int *ncase, int *icase, char *arg_dir, char *cafile)//, int *max_it, double *tol) 
{
  Fd_t fd;
  Off_t size;
  long lsize;
  int me;
  char *fbuffer;
  int p;

   MPI_Comm_rank(MPI_COMM_WORLD, &me);
   if(me==ROOT){
     fd = elio_open(fname, ELIO_R, ELIO_PRIVATE);
      if(!fd)  my_error(me, "no such config input file", 0);
      if(elio_length(fd, &size) != ELIO_OK) my_error(me,"error length",0); 
      lsize = (long)size;
   }
    MPI_Bcast(&lsize,1,MPI_LONG, ROOT, MPI_COMM_WORLD);
    if(lsize<1 || lsize > INT_MAX)  my_error(me," bad config file size",lsize);
    fbuffer = (char *) malloc((lsize+1)*sizeof(char));
    if(!fbuffer) my_error(me,"could not allocate config buffer",lsize);
    bzero(fbuffer, sizeof(fbuffer));
    if(me == ROOT) {
      size = elio_read(fd,0,fbuffer,(Size_t)lsize);
      if((long)size != lsize) my_error(me, "wrong size read", (long)size);
      elio_close(fd);
    }
      MPI_Bcast(fbuffer,(int)lsize, MPI_BYTE, ROOT, MPI_COMM_WORLD);
      getconfig(fbuffer, fname2, M, NC, caopt, outflag, ncase, icase, arg_dir, cafile);//,max_it,tol );
   free(fbuffer);
}

void getca2(char *fbuff, int *sel_brid, int *nca)
{
char *ptr=fbuff, *tmp;
int  i, i2;
   sscanf(ptr,"%d",&(*nca));
   ptr=nextline(ptr);
   for (i=0; i<*nca;i++)
   {    
        tmp = ptr;
        sscanf(ptr,"%d",&i2);
	sel_brid[i]=i2;
        ptr=nextline(ptr);
   }	    
}

/* new readcalist for reading branch id only */
void  readcalist2_(char *fname, int *sel_brid, int *nca) 
{
  Fd_t fd;
  Off_t size;
  long lsize;
  int me;
  char *fbuffer;
  int p;

   MPI_Comm_rank(MPI_COMM_WORLD, &me);
   if(me==ROOT){
     fd = elio_open(fname, ELIO_R, ELIO_PRIVATE);
      if(!fd)  my_error(me, "no such calist input file", 0);
      printf ("%s\n",fname);
      if(elio_length(fd, &size) != ELIO_OK) my_error(me,"error length",0); 
      lsize = (long)size;
      printf("file size is %ld\n",(long)size);
   }
    MPI_Bcast(&lsize,1,MPI_LONG, ROOT, MPI_COMM_WORLD);
    if(lsize<1 || lsize > INT_MAX)  my_error(me," bad file size",lsize);
    fbuffer = (char *) malloc((lsize+1)*sizeof(char));
    if(!fbuffer) my_error(me,"could not allocate buffer",lsize);
    bzero(fbuffer, sizeof(fbuffer));
    if(me == ROOT) {
      size = elio_read(fd,0,fbuffer,(Size_t)lsize);
      if((long)size != lsize) my_error(me, "wrong size read", (long)size);
      elio_close(fd);
    }
      MPI_Bcast(fbuffer,(int)lsize, MPI_BYTE, ROOT, MPI_COMM_WORLD);
      getca2(fbuffer, sel_brid, nca);
    free(fbuffer);
}

void getcanew(char *fbuff, int *CA_Index, int *N_of_CA, char ca_name[][5], int *ca_brid, int *ca_from, int *ca_to, int *ca_irow, int *nca)
{
   char *ptr=fbuff, *tmp, *tmpname;
   int  i, j, k,id, i1, i2, i3, icn,s;
   char c[5];

   sscanf(ptr,"%d",&(*nca));
   printf("nca = %d\n", *nca);
   ptr=nextline(ptr);
   icn = 0;
   s=0;
   for (i=0; i<*nca;i++)
   {    
        tmp = ptr;
        sscanf(ptr,"%d,%d,%s",&id,&k,c);
        printf("id = %d, k = %d, c =%s\n", id, k, c);
        s = s + k;
        for (j=0;j<5;j++)
//        for (j=0;j<strlen(c) ;j++)
	{
//         ca_name = c ;
         ca_name[i][j] = c[j] ;
        }
        printf("ca_name = %s\n", ca_name[i]);
        ca_name +=5;
        N_of_CA[i] = k;
	CA_Index[i] = id;
        ca_irow[i]=s;
        ptr=nextline(ptr);
        for (j=0; j<k; j++)
        {
          icn = icn + 1;
          sscanf(ptr,"%d,%d,%d,",&i1,&i2,&i3);
          ca_brid[icn]=i1;
          ca_from[icn]=i2;
          ca_to[icn]=i3;
          printf("caid = %d, ca_brid = %d, ca_from = %d, ca_to = %d\n", id, ca_brid[icn], ca_from[icn],ca_to[icn]);
          ptr=nextline(ptr);
        }
   }	    
}

/* new readcalist 8/29/2012 */
void  readcalistnew_(char *fname, int *CA_Index, int *N_of_CA, char ca_name[][5], int *ca_brid, int *ca_from, int *ca_to, int *ca_irow, int *nca) 
{
  Fd_t fd;
  Off_t size;
  long lsize;
  int me;
  char *fbuffer;
  int p;

   MPI_Comm_rank(MPI_COMM_WORLD, &me);
   if(me==ROOT){
     fd = elio_open(fname, ELIO_R, ELIO_PRIVATE);
      if(!fd)  my_error(me, "no such calist input file", 0);
      printf ("%s\n",fname);
      if(elio_length(fd, &size) != ELIO_OK) my_error(me,"error length",0); 
      lsize = (long)size;
      printf("file size is %ld\n",(long)size);
   }
    MPI_Bcast(&lsize,1,MPI_LONG, ROOT, MPI_COMM_WORLD);
    if(lsize<1 || lsize > INT_MAX)  my_error(me," bad file size",lsize);
    fbuffer = (char *) malloc((lsize+1)*sizeof(char));
    if(!fbuffer) my_error(me,"could not allocate buffer",lsize);
    bzero(fbuffer, sizeof(fbuffer));

    if(me == ROOT) {
      size = elio_read(fd,0,fbuffer,(Size_t)lsize);
      if((long)size != lsize) my_error(me, "wrong size read", (long)size);
      elio_close(fd);
    }
      MPI_Bcast(fbuffer,(int)lsize, MPI_BYTE, ROOT, MPI_COMM_WORLD);
      getcanew(fbuffer, CA_Index, N_of_CA, ca_name, ca_brid, ca_from, ca_to, ca_irow,nca);
   free(fbuffer);
}

