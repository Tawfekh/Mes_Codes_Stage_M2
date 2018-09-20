/***** ms.c     ************************************************
*
*       Generates samples of gametes ( theta given or fixed number
*						of segregating sites.)
*	Usage is shown by typing ms without arguments.
        usage: ms nsam howmany  -t  theta  [options]
		or
	       ms nsam howmany -s segsites  [options]

	   nsam is the number of gametes per sample.
	   howmany is the number of samples to produce.
	   With -t the numbers of segregating sites will randomly vary
		from one sample to the next.
	   with -s segsites,  the number of segregating sites will be
		segsites in each sample.

           Other options: See msdoc.pdf or after downloading and compiling, type ms<CR>.


*	  Arguments of the options are explained here:

           npop:  Number of subpopulations which make up the total population
           ni:  the sample size from the i th subpopulation (all must be
		specified.) The output will have the gametes in order such that
		the first n1 gametes are from the first island, the next n2 are
		from the second island, etc.
           nsites: number of sites between which recombination can occur.
           theta: 4No times the neutral mutation rate
           rho: recombination rate between ends of segment times 4No
	   f: ratio of conversion rate to recombination rate. (Wiuf and Hein model.)
	   track_len:  mean length of conversion track in units of sites.  The
		       total number of sites is nsites, specified with the -r option.
           mig_rate: migration rate: the fraction of each subpop made up of
                 migrants times 4No.
           howmany: howmany samples to generate.

	Note:  In the above definition, No is the total diploid population if
		npop is one, otherwise, No is the diploid population size of each
		subpopulation.
	A seed file called "seedms" will be created  if it doesn't exist. The
		seed(s) in this file will be modified by the program.
		So subsequent runs
		will produce new output.  The initial contents of seedms will be
		printed on the second line of the output.
        Output consists of one line with the command line arguments and one
	 	line with the seed(s).
		The samples appear sequentially following that line.
		Each sample begins with "//", then the number of segregating sites, the positions
		of the segregating sites (on a scale of 0.0 - 1.0). On the following
		lines are the sampled gametes, with mutants alleles represented as
		ones and ancestral alleles as zeros.
	To compile:  cc -o ms  ms.c  streec.c  rand1.c -lm
		or:  cc -o ms ms.c streec.c rand2.c -lm
	 (Of course, gcc would be used instead of cc on some machines.  And -O3 or
		some other optimization switches might be usefully employed with some
		compilers.) ( rand1.c uses drand48(), whereas rand2.c uses rand() ).

*
*   Modifications made to combine ms and mss on 25 Feb 2001
*	Modifications to command line options to use switches  25 Feb 2001
*	Modifications to add // before each sample  25 Feb 2001
	Modifications to add gene conversion 5 Mar 2001
	Added demographic options -d  13 Mar 2001
	Changed ran1() to use rand(). Changed seed i/o to accomodate this change. 20 April.
	Changed cleftr() to check for zero rand() .13 June 2001
	Move seed stuff to subroutine seedit()  11 July 2001
	Modified streec.c to handle zero length demographic intervals 9 Aug 2001
	Corrected problem with negative growth rates (Thanks to D. Posada and C. Wiuf) 13 May 2002
	Changed sample_stats.c to output thetah - pi rather than pi - thetah.  March 8 2003.
	Changed many command line options, allowing arbitrary migration matrix, and subpopulation
	   sizes.  Also allows parameters to come from a file. Option to output trees.  Option to
	   split and join subpopulations.   March 8, 2003. (Old versions saved in msold.tar ).
	!!! Fixed bug in -en option.  Earlier versions may have produced garbage when -en ... used. 9 Dec 2003
	Fixed bug which resulted in incorrect results for the case where
             rho = 0.0 and gene conversion rate > 0.0. This case was not handled
	    correctly in early versions of the program. 5 Apr 2004.  (Thanks to
	    Vincent Plagnol for pointing out this problem.)
	Fixed bug in prtree().  Earlier versions may have produced garbage when the -T option was used.
		 1 Jul 2004.
	Fixed bug in -e. options that caused problems with -f option  13 Aug 2004.
	Fixed bug in -es option, which was a problem when used with -eG. (Thanks again to V. Plagnol.) 6 Nov. 2004
	Added -F option:  -F minfreq  produces output with sites with minor allele freq < minfreq filtered out.  11 Nov. 2004.
	Fixed bug in streec.c (isseg() ).  Bug caused segmentation fault, crash on some machines. (Thanks
	    to Melissa Jane Todd-Hubisz for finding and reporting this bug.)
	Added -seeds option 4 Nov 2006
	Added "tbs" arguments feature 4 Nov 2006
	Added -L option.  10 May 2007
	Changed -ej option to set Mki = 0 pastward of the ej event.  See msdoc.pdf.  May 19 2007.
	fixed bug with memory allocation when using -I option. This caused problems expecially on Windows
          machines.  Thanks to several people, including Vitor Sousa and Stephane De Mita for help on this one.
          Oct. 17, 2007.
     Modified pickb() and pickbmf() to eliminate rare occurrence of fixed mutations Thanks to J. Liechty and K. Thornton. 4 Feb 2010.
	Added -p n  switch to allow position output to be higher precision.  10 Nov 2012.
     Changed spot from int to long in the function re().  29 July 2013.  ( Thanks to Yuri D'Elia for
        this suggestion.)
***************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "ms.h"
#define SITESINC 10

unsigned maxsites = SITESINC ;

struct node{
	int abv;
	int ndes;
	float time;
	int ncs; //number of descendents
	int masked; // have leaf masked or not. 0 means not masked, 1 means masked

	};

struct segl {
	int beg;
	struct node *ptree;
	int next;
	};

double *posit ;
double segfac ;
int count, ntbs, nseeds ;
struct params pars ;

comSS(char *srcSeq[], int *nsam, int *ncol, int *queMarLen, double *fiveSS) {
	char *tempSam;
	int hasEqualed(char* seq, char**seqList, int sizeOfSL);
    int hasEqualedQM(char* seq, char**seqList, int sizeOfSL);
	int observedH = 0, Npi=0, sLSize, i, j, tempDiff=0, ccc=0;
	double Hetero=0.0, sk2=0.0, mean_pi=0.0;
	double average_r2=0.0, average_sk2=0.0, Nmh=0.0, Nvd=0.0;
	void fiveStat(double fv[], char* seg1, char* seg2);
	char **cmatrix(), **segmentList;
	int diff(char* chr1, char* chr2);
	segmentList = cmatrix(nsam[0], ncol[0]+1);
	sLSize = nsam[0];
	memset(segmentList, 0, (unsigned)sLSize*sizeof( char* ) );

	if(*queMarLen > 0) {  // have question marks
		for(i=0;i<sLSize; i++) {
			tempSam = srcSeq[i];
			if(hasEqualedQM(tempSam, segmentList, sLSize) == 0) {
				segmentList[seqListSize(segmentList, sLSize)] = tempSam;
				observedH++;
			}
		}
	} else {
   
		for(i=0;i<sLSize; i++) {
			tempSam = srcSeq[i];
			if(hasEqualed(tempSam, segmentList, sLSize) == 0) {
				segmentList[seqListSize(segmentList, sLSize)] = tempSam;
				observedH++;
			}
		}
	}
	fiveSS[0] = observedH;
  
  
	Npi = 0;
    for(i=0;i<sLSize-1;i++) {
        for(j=i+1;j<sLSize;j++) {
            Npi++;
        }
    }
    int tempDiffArr[Npi];

    Npi = 0;
	Hetero = 0.0;
	mean_pi = 0.0;
	sk2 = 0.0;
	for(i=0;i<sLSize-1;i++) {
		for(j=i+1;j<sLSize;j++) {
			tempDiff = diff(srcSeq[i], srcSeq[j]);
			tempDiffArr[Npi++] = tempDiff;
			mean_pi = mean_pi + tempDiff;
			if (tempDiff > 0) { Hetero ++; }
		}
	}
	Hetero = Hetero/Npi;
	mean_pi = mean_pi/Npi;
	for(i=0;i<Npi;i++) {
		sk2 = sk2 + (tempDiffArr[i]-mean_pi)*(tempDiffArr[i]-mean_pi);
	}
	sk2 = sk2/(Npi-1);
	//fiveSS[1] = sk2; // sk2
	fiveSS[1] = Hetero; // Hetero

	// second initial to zero
	Npi = 0;
	average_r2 = 0.0;
	average_sk2 = 0.0;
	Nmh = 0.0;
	Nvd = 0.0;

	char tempSeg1[sLSize+1];
	for(i=0;i<sLSize+1;i++) {
	  tempSeg1[i] = '\0';
	}
	char tempSeg2[sLSize+1];
	for(i=0;i<sLSize+1;i++) {
	  tempSeg2[i] = '\0';
	}

	double fiveVal[2] = {0.0, 0.0};
	for(i=0;i<ncol[0]-1;i++) {
		for(j=i+1;j<ncol[0];j++) {
			if((srcSeq[0][i] != '-') && (srcSeq[0][j] != '-')) {
  				Npi++;
  				for(ccc=0;ccc<sLSize;ccc++) {
             tempSeg1[ccc] = srcSeq[ccc][i];
             tempSeg2[ccc] = srcSeq[ccc][j];
          }
          fiveStat(fiveVal, tempSeg1, tempSeg2);
          
          //printf("i:%d, j:%d, sk2: %f, r2: %f\n", i,j, fiveVal[1], fiveVal[0]);

  				average_r2 = average_r2 + fiveVal[0];
  				average_sk2 = average_sk2 + fiveVal[1];
          }
      }
  }

	  average_r2 = average_r2/Npi;
	  average_sk2 = average_sk2/Npi;
    fiveSS[2] = average_sk2; // average_sk2
    fiveSS[3] = average_r2; // average_r2
}

//main(argc,argv)
//        int argc;
 //       char *argv[];
mySimulator(char **argv, int *argvLen, int *yiWanHaps, int *quesMarIds, int *queMarLen, double *sixSS, int *oneOrSix)  // for simulator
{

    int argc = *argvLen;   // for simulator
    
    
    //printf("argc:%d\n", argc);
   
   
   
	int i, k, howmany, segsites, ton;
	char **list, **cmatrix(), **segmentList; //, **tbsparamstrs ;
	FILE *pf, *fopen() ;
	double probss, tmrca, ttot ;
	void seedit( const char * ) ;
	void getpars( int argc, char **argv, int *howmany )  ;
	int gensam( char **list, double *probss, double *ptmrca, double *pttot ) ;
	int gensamWithQM( char **list, double *probss, double *ptmrca, double *pttot, int *qmRow, int *qmCol, int queMarLen) ;
    int sLSize, dataCount, yiWanH;
    int j, Npi=0, tempDiff=0, ccc=0;
	double Hetero=0.0, sk2=0.0, mean_pi=0.0;
	double average_r2=0.0, average_sk2=0.0, Nmh=0.0, Nvd=0.0;
    char *tempSam;
    int tSEquOt(char* seg1, char* seg2);
    int hasEqualed(char* seq, char**seqList, int sizeOfSL);
    int hasEqualedQM(char* seq, char**seqList, int sizeOfSL);
    int seqListSize(char**seqList, int sizeOfSL);
    int diff(char* chr1, char* chr2);
	void fiveStat(double fv[], char* seg1, char* seg2);
    ntbs = 0 ;   /* these next few lines are for reading in parameters from a file (for each sample) */
/*	tbsparamstrs = (char **)malloc( argc*sizeof(char *) ) ;

//	for( i=0; i<argc; i++) printf("%s ",argv[i]);  for only data
	for( i =0; i<argc; i++) tbsparamstrs[i] = (char *)malloc(30*sizeof(char) ) ;
	for( i = 1; i<argc ; i++)
			if( strcmp( argv[i],"tbs") == 0 )  argv[i] = tbsparamstrs[ ntbs++] ;
*/
	count=0;


//	if( ntbs > 0 )  for( k=0; k<ntbs; k++)  scanf(" %s", tbsparamstrs[k] );
    getpars( argc, argv, &howmany) ;   /* results are stored in global variable, pars */
// printf("nsam:%d,howmany:%d,singleton:%d,doubleton:%d,otherton:%d,rho:%f\n",pars.cp.nsam,howmany, pars.mfs.singlton,pars.mfs.doubleton,pars.mfs.otherton, pars.cp.r);
 // printf("nsam:%d,howmany:%d,g:%f\n",pars.cp.nsam,howmany, pars.cp.alphag[0]);
 
  sLSize = pars.cp.nsam;
    pars.mfs.strNumber = oneOrSix[0];

   // printf("1 argc:%d\n", argc);
	if( !pars.commandlineseedflag ) seedit("s");
	pf = stdout ;

    if( pars.mp.segsitesin ==  0 ) { // not have -s
        if((pars.mfs.singlton == 0) && (pars.mfs.doubleton == 0) && (pars.mfs.otherton == 0)) { // not have mfs -Y
            list = cmatrix(pars.cp.nsam,maxsites+1);
            posit = (double *)malloc( (unsigned)( maxsites*sizeof( double)) ) ;
        } else {
            ton = pars.mfs.doubleton + pars.mfs.singlton + pars.mfs.otherton;
            if(pars.mfs.moreLength > 0) {
                list = cmatrix(pars.cp.nsam*2, ton+1); //increase the size two times
            } else {
                list = cmatrix(pars.cp.nsam, ton+1);
            }
            segmentList = cmatrix(pars.cp.nsam, ton+1);
            posit = (double *)malloc( (unsigned)( ton*sizeof( double)) ) ;
            if( pars.mp.theta > 0.0 ){
                segfac = 1.0 ;
                for(  i= ton; i > 1; i--) segfac *= i ;
            }
        }
	} else { //have -s
	    list = cmatrix(pars.cp.nsam, pars.mp.segsitesin+1 ) ;
        posit = (double *)malloc( (unsigned)( pars.mp.segsitesin*sizeof( double)) ) ;
        if( pars.mp.theta > 0.0 ){
             segfac = 1.0 ;
             for(  i= pars.mp.segsitesin; i > 1; i--) segfac *= i ;
        }
	}
  
  
    //printf("2 argc:%d\n", argc);
	//fprintf(pf,"befff[0]: %d\n", 7);
//fprintf(pf,"queMarLen[0]: %d\n", *queMarLen);
int testLen = *queMarLen;
int qmRow[testLen]; // question marks row index
int qmCol[testLen]; // question marks column index

if(*queMarLen > 0) {  // have question marks
        //int testLen = 6;
        //int qMInd[] = {0, 1, 2, 3, 4 , 5};

        for(dataCount=0;dataCount<testLen;dataCount++) {
           // int qmTmp = atoi(queMarInd[dataCount]);  // atoi(queMarInd[dataCount]); is String to int
            int qmTmp = quesMarIds[dataCount]; //int qmTmp = qMInd[dataCount];
            //fprintf(pf,"markIndex: %d\n", qmTmp);
            qmRow[dataCount] = qmTmp%sLSize;
            qmCol[dataCount] = qmTmp/sLSize;
        }
}


    dataCount = 0;
    while( howmany-count++ ) {
  //          printf("dataCount:%d",  dataCount);
	/*   if( (ntbs > 0) && (count >1 ) ){
	         for( k=0; k<ntbs; k++){
			    if( scanf(" %s", tbsparamstrs[k]) == EOF ){
			       if( !pars.commandlineseedflag ) seedit( "end" );
				   exit(0);
				}
			 }
			 getpars( argc, argv, &howmany) ;
	   } */
 //      fprintf(pf,"\n//");     for only data
	/*   if( ntbs >0 ){
			for(k=0; k< ntbs; k++) printf("\t%s", tbsparamstrs[k] ) ;
		} */
//		printf("\n");  for only data



   // printf("3 argc:%d\n", argc);
//printf("nsam:%d,howmany:%d,singleton:%d,doubleton:%d,otherton:%d,rho:%f\n",pars.cp.nsam,howmany, pars.mfs.singlton,pars.mfs.doubleton,pars.mfs.otherton, pars.cp.r);



if(*queMarLen > 0) {  // have question marks
    segsites = gensamWithQM( list, &probss, &tmrca, &ttot, qmRow, qmCol, testLen) ;
} else {
    segsites = gensam( list, &probss, &tmrca, &ttot ) ;
}

      /*  if( pars.mp.timeflag ) fprintf(pf,"time:\t%lf\t%lf\n",tmrca, ttot ) ;
        if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
   	       if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 ))
		       fprintf(pf,"prob: %g\n", probss ) ;
           fprintf(pf,"segsites: %d\n",segsites);  //for only data
		   if( segsites > 0 )	fprintf(pf,"positions: ");   // for only data
		   for( i=0; i<segsites; i++)
              fprintf(pf,"%6.*lf ", pars.output_precision,posit[i] ); // for only data
           fprintf(pf,"\n");// for only data 
	       if( segsites > 0 ) { */
	         // for(i=0;i<pars.cp.nsam; i++) { fprintf(pf,"%s\n", list[i] ); }
	   //    }
	  //  } 



  //  printf("4 argc:%d\n", argc);

//printf("gensam before:%d", 7);
        yiWanH = 0;
        //clear
        memset(segmentList, 0, (unsigned)sLSize*sizeof( char* ) );



       // if(1) {
       if(*queMarLen > 0) {  // have question marks
            for(i=0;i<sLSize; i++) {
                tempSam = list[i];
                if(hasEqualedQM(tempSam, segmentList, sLSize) == 0) {
                    segmentList[seqListSize(segmentList, sLSize)] = tempSam;
                    yiWanH++;
                }
            }
        } else {
            for(i=0;i<sLSize; i++) {
                tempSam = list[i];
                if(hasEqualed(tempSam, segmentList, sLSize) == 0) {
                    segmentList[seqListSize(segmentList, sLSize)] = tempSam;
                    yiWanH++;
                }
            }
        }
        
        
   // printf("5 argc:%d\n", argc);

        //fprintf(pf,"H: %d\n", yiWanH);
        //if(*queMarLen > 0) {  // have question marks only H
        //    yiWanHaps[dataCount++] = yiWanH;
        //} else {
            if(oneOrSix[0] == 1) {
                yiWanHaps[dataCount++] = yiWanH;   // for simulator
            } else {
                sixSS[dataCount++] = yiWanH; // H

                Npi = 0;
                for(i=0;i<sLSize-1;i++) {
                    for(j=i+1;j<sLSize;j++) {
                        Npi++;
                    }
                }
                int tempDiffArr[Npi];

                Npi = 0;
                Hetero = 0.0;
                mean_pi = 0.0;
                sk2 = 0.0;
                for(i=0;i<sLSize-1;i++) {
                    for(j=i+1;j<sLSize;j++) {
                        tempDiff = diff(list[i], list[j]);
                        tempDiffArr[Npi++] = tempDiff;
                        mean_pi = mean_pi + tempDiff;
                        if (tempDiff > 0) { Hetero ++; }
                    }
                }
                Hetero = Hetero/Npi;
                mean_pi = mean_pi/Npi;
                
				/*for(i=0;i<Npi;i++) {
                    sk2 = sk2 + (tempDiffArr[i]-mean_pi)*(tempDiffArr[i]-mean_pi);
                }
                sk2 = sk2/(Npi-1);
                sixSS[dataCount++] = sk2; // sk2
				*/
				
				
                sixSS[dataCount++] = Hetero; // Hetero

                // second initial to zero
                Npi = 0;
                average_r2 = 0.0;
                average_sk2 = 0.0;
                Nmh = 0.0;
                Nvd = 0.0;

                char tempSeg1[sLSize+1];
                for(i=0;i<sLSize+1;i++) {
                  tempSeg1[i] = '\0';
                }
                char tempSeg2[sLSize+1];
                for(i=0;i<sLSize+1;i++) {
                  tempSeg2[i] = '\0';
                }
                //double fiveVal[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
				double fiveVal[2] = {0.0, 0.0};
                for(i=0;i<segsites-1;i++) {
                    for(j=i+1;j<segsites;j++) {
                        if((list[0][i] != '-') && (list[0][j] != '-')) {
                            Npi++;
                            for(ccc=0;ccc<sLSize;ccc++) {
                                //if(list[ccc][i] == '2') {
                                //        tempSeg1[ccc] = '1';
                                //    } else {
                                        tempSeg1[ccc] = list[ccc][i];
                                //   }
                                    //fprintf(pf, "tempSeg1: %s\n", tempSeg1);
                                //    if(list[ccc][j] == '2') {
                                //        tempSeg2[ccc] = '1';
                                //    } else {
                                        tempSeg2[ccc] = list[ccc][j];
                                //    }
                                    //fprintf(pf, "tempSeg2: %s\n", tempSeg2);
                            }
                                //printf("strlent strlen(tempSeg1): %d \n", strlen(tempSeg1));

                            fiveStat(fiveVal, tempSeg1, tempSeg2);

                                //	fprintf(pf,"r2: %f, var_pi: %f, LD1: %f,very_diverse_flag: %f, missing_haplotype_flag: %f\n",fiveVal[0], fiveVal[1], fiveVal[2], fiveVal[3], fiveVal[4]);


                            average_r2 = average_r2 + fiveVal[0];
                            average_sk2 = average_sk2 + fiveVal[1];
                            /*if(fiveVal[3] == 1.0) {

                                Nvd = Nvd + 1.0;
                                if(fiveVal[4] == 1.0) {
                                    Nmh = Nmh + 1.0;
                                }
                            }*/
                        }
                    }
                }

                average_r2 = average_r2/Npi;
                average_sk2 = average_sk2/Npi;

                //if(Nvd == 0.0) {
                //    average_LD1 = 1.0;
               // } else {
               //     average_LD1 = average_LD1/Nvd;
               // }

                sixSS[dataCount++] = average_sk2; // average_sk2
                sixSS[dataCount++] = average_r2; // average_r2
                //sixSS[dataCount++] = average_LD1; // average_LD1
            }
        //}


    }
	if( !pars.commandlineseedflag ) seedit( "end" );


}

int hasEqualed(char* seq, char**seqList, int sizeOfSL) {
    int flag = 0;
    int i = 0;
    for(i=0;i<sizeOfSL; i++) {
        if(seqList[i] != NULL) {
            if(strcmp(seq, seqList[i]) == 0) {
                flag = 1;
                break;
            }
        } else {
            break;
        }
    }
    return flag;
}

int hasEqualedQM(char* seq, char**seqList, int sizeOfSL) {
    int flag = 0;
    int i = 0;
    for(i=0;i<sizeOfSL; i++) {
        if(seqList[i] != NULL) {
            /*if(strcmp(seq, seqList[i]) == 0) {  // // if equal, return 1
                flag = 1;
                break;
            }*/
			if(tSEquOt(seq, seqList[i]) == 0) {  // // if equal, return 1
                flag = 1;
            }
        } else {
            break;
        }

    }
    return flag;
}

int tSEquOt(char* seg1, char* seg2) { // similiar function to strcmp
	int flag = 0;
	int i = 0;
	int len = strlen(seg1);
	for(i=0;i<len;i++) {
		if((seg1[i] == '?') || (seg2[i] == '?')) {
			continue;
		}

		if(seg1[i] != seg2[i]) {
			flag = 1;
			break;
		}
	}
	return flag;
}

int seqListSize(char**seqList, int sizeOfSL) {
    int noNullSize = 0;
    int i = 0;
    for(i=0;i<sizeOfSL; i++) {
        if(seqList[i] != NULL) {
            noNullSize++;
        } else {
            break;
        }
    }
    return noNullSize;
}


void fiveStat(double fv[], char* seg1, char* seg2) {
	int len = strlen(seg1);
	//printf("seg1: %s", seg1);
	//printf("seg2: %s", seg2);
	//printf("strlent len: %d \n", len);


	int i=0;
	double n00=0.0, n01=0.0, n10=0.0, n11=0.0;
	double p00=0.0, p01=0.0, p10=0.0, p11=0.0, p0_=0.0, p_0=0.0, p1_=0.0, p_1=0.0;
	double D=0.0, r2=0.0, mean_pi=0.0, var_pi=0.0, LD1=0.0, very_diverse_flag=0.0, missing_haplotype_flag=0.0;
	for(i=0;i<len;i++) {
		if((seg1[i] == '0') && (seg2[i] == '0')) {
			n00 = n00 + 1.0;
		} else if ((seg1[i] == '0') && (seg2[i] == '1')) {
			n01 = n01 + 1.0;
		} else if ((seg1[i] == '1') && (seg2[i] == '0')) {
			n10 = n10 + 1.0;
		} else if ((seg1[i] == '1') && (seg2[i] == '1')) {
			n11 = n11 + 1.0;
		} else if ((seg1[i] == '?') && (seg2[i] == '0')) {
			n00 <- n00 + 1;
		} else if ((seg1[i] == '?') && (seg2[i] == '1')) {
			n11 <- n11 + 1;
		} else if ((seg1[i] == '0') && (seg2[i] == '?')) {
			n00 <- n00 + 1;
		} else if ((seg1[i] == '1') && (seg2[i] == '?')) {
			n11 <- n11 + 1;
		}
	}

	//printf("n00: %f, n01: %f, n10: %f, n11: %f\n", n00, n01, n10, n11);

	p00 = n00/len;

	//printf("n00: %f, len: %f, p00: %f, n00/len: %f \n", n00, len, p00, n00/len);

	p01 = n01/len;
	p10 = n10/len;
	p11 = n11/len;

	//printf("p00: %f, p01: %f, p10: %f, p11: %f\n", p00, p01, p10, p11);

	p0_ = (n00+n01)/len;
	p_0 = (n00+n10)/len;
	p1_ = 1.0 - p0_;
	p_1 = 1.0 - p_0;
	D = p00 - p0_*p_0;



	if(D == 0.0) {
		r2 = 0.0;
	} else {
		r2 = D*D/(p0_*p1_*p_0*p_1);
	}

	//printf("r2: %f\n", r2);

	mean_pi = (n00*n01 + n00*n10 + 2*n00*n11 + 2*n01*n10 + n01*n11 + n10*n11) / (len*(len-1)/2);

	//printf("mean_pi: %f\n", mean_pi);

	var_pi = ((n00*(n00-1)/2 + n01*(n01-1)/2 + n10*(n10-1)/2 + n11*(n11-1)/2) * (0-mean_pi)*(0-mean_pi)+
               n00*n01*(1-mean_pi)*(1-mean_pi)+
               n00*n10*(1-mean_pi)*(1-mean_pi)+
               n00*n11*(2-mean_pi)*(2-mean_pi)+
               n01*n10*(2-mean_pi)*(2-mean_pi)+
               n01*n11*(1-mean_pi)*(1-mean_pi)+
               n10*n11*(1-mean_pi)*(1-mean_pi)) / (len*(len-1)/2-1);


	//printf("var_pi: %f\n", var_pi);


/*	if((p00<=p01) && (p00<=p10) && (p00<=p11)) {
		if((p00==0.0) || (p0_==0.0) || (p_0==0.0)) {
			LD1 = 1.0;
		} else {
			LD1 = 1.0 - p00/p0_/p_0;
		}
	} else if((p01<=p10) && (p01<=p11)) {
		if((p01==0.0) || (p0_==0.0) || (p_1==0.0)) {
			LD1 = 1.0;
		} else {
			LD1 = 1.0 - p01/p0_/p_1;
		}
	} else if(p10<=p11) {
		if((p10==0.0) || (p1_==0.0) || (p_0==0.0)) {
			LD1 = 1.0;
		} else {
			LD1 = 1.0 - p10/p1_/p_0;
		}
	} else {
		if((p11==0.0) || (p1_==0.0) || (p_1==0.0)) {
		  LD1 = 1.0;
		} else {
		  LD1 = 1.0 - p11/p1_/p_1;
		}
	}
*/
	//printf("LD1: %f\n", LD1);

	//if((p0_>0.25) && (p1_>0.25) && (p_0>0.25) && (p_1>0.25)) {
	//	very_diverse_flag = 1.0;
	//}

	//if((p00<0.04) || (p01<0.04) || (p10<0.04) || (p11<0.04)) {
	//	missing_haplotype_flag = 1.0;
	//}

	fv[0] = r2;
	fv[1] = var_pi;
	//fv[2] = LD1;
	//fv[3] = very_diverse_flag;
	//fv[4] = missing_haplotype_flag;
}

int diff(char* chr1, char* chr2) {
	int diffNum = 0;
	int len = strlen(chr1);
	int i=0;
	for(i=0;i<len;i++) {
		if((chr1[i] == '?') || (chr2[i] == '?')) {

		} else if(chr1[i] != chr2[i]) {
			diffNum++;
		}
	}
	return diffNum;
}


	int
gensam( char **list, double *pprobss, double *ptmrca, double *pttot )
{
	int nsegs, h, i, k, j, seg, ns, start, end, len, segsit ;
	struct segl *seglst, *segtre_mig(struct c_params *p, int *nsegs ) ; /* used to be: [MAXSEG];  */
	double nsinv,  tseg, tt, ttime(struct node *, int nsam), ttimemf(struct node *, int nsam, int mfreq) ;

	//stt is singleton total time of all trees, dtt is doubleton total time of all trees
	//ott is otherton total time of all trees
	double stt, dtt, ott;

   	double *pk;
    // spk[k] stores total time of total tips in each segment tree, 0=<k<nsegs
	// dpk[k] stores total time of total nodes (which has two desendents) in each segment tree, 0=<k<nsegs
	// opk[k] stores total time of total nodes (which has >= three desendents) in each segment tree, 0=<k<nsegs
	double *spk, *dpk, *opk;

	int *ss;
	int *sdouton, *ssinton, *sothton;
    int douton, sinton, othton;
	int segsitesin,nsites, mhlen;
	double theta, es ;
	int nsam, mfreq ;
	void prtree( struct node *ptree, int nsam);
	int make_gametes(int nsam, int mfreq, double *posit, long int nsites, struct node *ptree, double tt, int newsites, int ns, char **list );
    int make_MFSgametes(int nsam, int mfreq, double *posit, long int nsites, struct node *ptree, double sintontt, double doutontt, double othtontt, int stonnewsites, int dtonnewsites, int otonnewsites, int ns, char **list, int loop);

 	void ndes_setup( struct node *, int nsam );

    int totalSeg;
    void sdottime(struct node *, int nsam, double* sdot, int unfo);
    void sdottimemf(struct node *, int nsam, double* sdot, int mfreq, int unfo);
    double sdotime[3];

	nsites = pars.cp.nsites ;
	nsinv = 1./nsites;

	//clock_t starts, finishs;
	//starts = clock();
    seglst = segtre_mig(&(pars.cp),  &nsegs ) ;
  //  finishs = clock();
  //  double duration = (double)(finishs - starts) / CLOCKS_PER_SEC;
 //   printf( "%f seconds\n", duration );



	nsam = pars.cp.nsam;
	segsitesin = pars.mp.segsitesin ;
	theta = pars.mp.theta ;
	mfreq = pars.mp.mfreq ;
    mhlen = pars.mfs.mhlength;

	if( pars.mp.treeflag ) {
	  	ns = 0 ;
	    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
	      if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) ){
		     end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		     start = seglst[seg].beg ;
		     len = end - start + 1 ;
		     fprintf(stdout,"[%d]", len);
	      }
	      prtree( seglst[seg].ptree, nsam ) ;
	      if( (segsitesin == 0) && ( theta == 0.0 ) && ( pars.mp.timeflag == 0 ) )
	  	      free(seglst[seg].ptree) ;
	    }
	}

	if( pars.mp.timeflag ) {
      tt = 0.0 ;
	  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
		if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		if( (nsegs==1) || ( ( start <= nsites/2) && ( end >= nsites/2 ) ) )
		  *ptmrca = (seglst[seg].ptree + 2*nsam-2) -> time ;
		len = end - start + 1 ;
		tseg = len/(double)nsites ;
		if( mfreq == 1 ) tt += ttime(seglst[seg].ptree,nsam)*tseg ;
		else tt += ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
		if( (segsitesin == 0) && ( theta == 0.0 )  )
	  	      free(seglst[seg].ptree) ;
	    }
		*pttot = tt ;
	 }

    if(segsitesin == 0) {
        if((theta > 0.0) && ((pars.mfs.singlton == 0) && (pars.mfs.doubleton == 0) && (pars.mfs.otherton == 0))) {
            ns = 0 ;
            for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
                if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
                end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
                start = seglst[seg].beg ;
                len = end - start + 1 ;
                tseg = len*(theta/nsites) ;
                if( mfreq == 1) tt = ttime(seglst[seg].ptree, nsam);
                        else tt = ttimemf(seglst[seg].ptree, nsam, mfreq );
                segsit = poisso( tseg*tt );
                if( (segsit + ns) >= maxsites ) {
                    maxsites = segsit + ns + SITESINC ;
                    posit = (double *)realloc(posit, maxsites*sizeof(double) ) ;
                      biggerlist(nsam, list) ;
                }
                make_gametes(nsam,mfreq, posit, mhlen, seglst[seg].ptree,tt, segsit, ns, list );
                free(seglst[seg].ptree) ;
                locate(segsit,start*nsinv, len*nsinv,posit+ns);
                ns += segsit;
            }
        } else if((pars.mfs.singlton > 0) || (pars.mfs.doubleton > 0) || (pars.mfs.otherton > 0)) {

             douton = pars.mfs.doubleton;
             sinton = pars.mfs.singlton;
             othton = pars.mfs.otherton;

             // for mfs
             spk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
             dpk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
             opk = (double *)malloc((unsigned)(nsegs*sizeof(double)));

             sdouton = (int *)malloc((unsigned)(nsegs*sizeof(int)));
             ssinton = (int *)malloc((unsigned)(nsegs*sizeof(int)));
             sothton = (int *)malloc((unsigned)(nsegs*sizeof(int)));

             if( (spk==NULL) || (dpk==NULL) || (opk==NULL) || (ssinton==NULL) || (sdouton==NULL) || (sothton==NULL)) perror("malloc error. gensam.2");

             stt = 0.0;
             dtt = 0.0;
             ott = 0.0;

             for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {

                if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
                end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
                start = seglst[seg].beg ;
                len = end - start + 1 ;
                tseg = len/(double)nsites ;

                if(mfreq == 1) {
                    sdottime(seglst[seg].ptree,nsam, sdotime, pars.mfs.unfolded);
                } else {
                    sdottimemf(seglst[seg].ptree,nsam, sdotime, mfreq, pars.mfs.unfolded);
                }
               // printf("stime:%f, dtime:%f, xtime:%f\n", sdotime[0], sdotime[1], sdotime[2] );
                spk[k] = sdotime[0]*tseg;
                dpk[k] = sdotime[1]*tseg;
                opk[k] = sdotime[2]*tseg;
                stt += spk[k];
                dtt += dpk[k];
                ott += opk[k];


              }

              if( theta > 0.0 ) {
                es = theta * (stt+dtt+ott);
                *pprobss = exp( -es )*pow( es, (double)(douton+sinton+othton)) / segfac ;
              }

              //for singlton
	          if(stt > 0.0) {
                for (k=0;k<nsegs;k++) spk[k] /= stt;
                //init ssinton, which stores each singleton number in each segment tree
                mnmial(sinton,nsegs,spk,ssinton);
              } else {
                 for( k=0; k<nsegs; k++) ssinton[k] = 0 ;
              }

              //for doubleton
              if(dtt > 0.0) {
                for (k=0;k<nsegs;k++) dpk[k] /= dtt;
                //init douton, which stores each doubleton number in each segment tree
                mnmial(douton,nsegs,dpk,sdouton);
              } else {
                for( k=0; k<nsegs; k++) sdouton[k] = 0 ;
              }

              //for otherton
              if(ott > 0.0) {
                for (k=0;k<nsegs;k++) opk[k] /= ott;
                mnmial(othton,nsegs,opk,sothton);
              } else {
                for( k=0; k<nsegs; k++) sothton[k] = 0 ;
              }

              ns = 0 ;
              for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
                 end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
                 start = seglst[seg].beg ;
                 len = end - start + 1 ;
                 tseg = len/(double)nsites;
                 totalSeg = ssinton[k] + sdouton[k] + sothton[k];

                 locate(totalSeg, start*nsinv, len*nsinv, posit+ns);

                 make_MFSgametes(nsam, mfreq, posit, mhlen, seglst[seg].ptree, stt*spk[k]/tseg, dtt*dpk[k]/tseg, ott*opk[k]/tseg, ssinton[k], sdouton[k], sothton[k], ns, list, 0);

                 ns = ns + totalSeg;
                 free(seglst[seg].ptree) ;
              }

              free(spk);
              free(dpk);
              free(opk);
              free(ssinton);
              free(sdouton);
              free(sothton);
        }
    } else if( segsitesin > 0 ) {

        pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
        ss = (int *)malloc((unsigned)(nsegs*sizeof(int)));
        if( (pk==NULL) || (ss==NULL) ) perror("malloc error. gensam.2");


        tt = 0.0 ;
      for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
        if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
        end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
        start = seglst[seg].beg ;
        len = end - start + 1 ;
        tseg = len/(double)nsites ;
               if( mfreq == 1 ) pk[k] = ttime(seglst[seg].ptree,nsam)*tseg ;
               else pk[k] = ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
                 tt += pk[k] ;
      }
	  if( theta > 0.0 ) {
	    es = theta * tt ;
	    *pprobss = exp( -es )*pow( es, (double) segsitesin) / segfac ;
	  }
	  if( tt > 0.0 ) {
          for (k=0;k<nsegs;k++) pk[k] /= tt ;
          mnmial(segsitesin,nsegs,pk,ss);
	  }
	  else
            for( k=0; k<nsegs; k++) ss[k] = 0 ;
	  ns = 0 ;
	  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
		 end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		 start = seglst[seg].beg ;
		 len = end - start + 1 ;
		 tseg = len/(double)nsites;

		 make_gametes(nsam,mfreq, posit, mhlen, seglst[seg].ptree,tt*pk[k]/tseg, ss[k], ns, list);
		 free(seglst[seg].ptree) ;
		 locate(ss[k],start*nsinv, len*nsinv,posit+ns);

		 ns += ss[k] ;
	  }
	  free(pk);
	  free(ss);

    }
	for(i=0;i<nsam;i++) list[i][ns] = '\0' ;
	if(pars.mfs.moreLength > 0) {
        for(i=nsam;i<nsam*2;i++) list[i][ns] = '\0' ;
	}
	return( ns ) ;
}



int
gensamWithQM( char **list, double *pprobss, double *ptmrca, double *pttot, int *qmRow, int *qmCol, int queMarLen )
{
	int nsegs, h, i, k, j, seg, ns, start, end, len, segsit ;
	struct segl *seglst, *segtre_mig(struct c_params *p, int *nsegs ) ; /* used to be: [MAXSEG];  */
	double nsinv,  tseg, tt, ttime(struct node *, int nsam), ttimemf(struct node *, int nsam, int mfreq) ;

	double stt, dtt, ott;

   	double *pk;

	double *spk, *dpk, *opk;

	int *ss;
	int *sdouton, *ssinton, *sothton;
    int douton, sinton, othton;
	int segsitesin,nsites, mhlen;
	double theta, es ;
	int nsam, mfreq ;
	void prtree( struct node *ptree, int nsam);
	int make_MFSgametesQM(int nsam, double *posit, long int mhlen, struct node *ptree, double* sdoTemTime, int* sdoTemSite, int ns, char **list, int* currTreRow, int* currTreCol, int currTreLen);

 	void ndes_setup( struct node *, int nsam );

    int totalSeg;
    void sdottime(struct node *, int nsam, double* sdot, int unfo);
    void sdottimemf(struct node *, int nsam, double* sdot, int mfreq, int unfo);
    double sdotime[3];
    int currTreQM, qmCount, currCount;

	nsites = pars.cp.nsites ;
	nsinv = 1./nsites;


    seglst = segtre_mig(&(pars.cp),  &nsegs ) ;

	nsam = pars.cp.nsam;
	segsitesin = pars.mp.segsitesin ;
	theta = pars.mp.theta ;
	mfreq = pars.mp.mfreq ;
    mhlen = pars.mfs.mhlength;

    douton = pars.mfs.doubleton;
    sinton = pars.mfs.singlton;
    othton = pars.mfs.otherton;

    // for mfs
    spk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
    dpk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
    opk = (double *)malloc((unsigned)(nsegs*sizeof(double)));

    sdouton = (int *)malloc((unsigned)(nsegs*sizeof(int)));
    ssinton = (int *)malloc((unsigned)(nsegs*sizeof(int)));
    sothton = (int *)malloc((unsigned)(nsegs*sizeof(int)));

    if( (spk==NULL) || (dpk==NULL) || (opk==NULL) || (ssinton==NULL) || (sdouton==NULL) || (sothton==NULL)) perror("malloc error. gensam.2");

     stt = 0.0;
     dtt = 0.0;
     ott = 0.0;

     for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {

        if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
        end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
        start = seglst[seg].beg ;
        len = end - start + 1 ;
        tseg = len/(double)nsites ;

        if(mfreq == 1) {
            sdottime(seglst[seg].ptree,nsam, sdotime, pars.mfs.unfolded);
        } else {
            sdottimemf(seglst[seg].ptree,nsam, sdotime, mfreq, pars.mfs.unfolded);
        }
       // printf("stime:%f, dtime:%f, xtime:%f\n", sdotime[0], sdotime[1], sdotime[2] );
        spk[k] = sdotime[0]*tseg;
        dpk[k] = sdotime[1]*tseg;
        opk[k] = sdotime[2]*tseg;
        stt += spk[k];
        dtt += dpk[k];
        ott += opk[k];
      }

      //for singlton
      if(stt > 0.0) {
        for (k=0;k<nsegs;k++) spk[k] /= stt;
        //init ssinton, which stores each singleton number in each segment tree
        mnmial(sinton,nsegs,spk,ssinton);
      } else {
         for( k=0; k<nsegs; k++) ssinton[k] = 0 ;
      }

      //for doubleton
      if(dtt > 0.0) {
        for (k=0;k<nsegs;k++) dpk[k] /= dtt;
        //init douton, which stores each doubleton number in each segment tree
        mnmial(douton,nsegs,dpk,sdouton);
      } else {
        for( k=0; k<nsegs; k++) sdouton[k] = 0 ;
      }

      //for otherton
      if(ott > 0.0) {
        for (k=0;k<nsegs;k++) opk[k] /= ott;
        mnmial(othton,nsegs,opk,sothton);
      } else {
        for( k=0; k<nsegs; k++) sothton[k] = 0 ;
      }

      ns = 0 ;
      double sdoTemTime[3] = {0.0, 0.0, 0.0};
      int sdoTemSite[3] = {0, 0, 0};

      for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
         end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
         start = seglst[seg].beg ;
         len = end - start + 1 ;
         tseg = len/(double)nsites;
         totalSeg = ssinton[k] + sdouton[k] + sothton[k];


         currTreQM = 0;
         for(qmCount=0;qmCount<queMarLen;qmCount++) {
            if((qmCol[qmCount] >= ns) && (qmCol[qmCount] < (ns+totalSeg))) {
                currTreQM++;
            }
         }
         int currTreRow[currTreQM];
         int currTreCol[currTreQM];
         currCount = 0;
         for(qmCount=0;qmCount<queMarLen;qmCount++) {
            if((qmCol[qmCount] >= ns) && (qmCol[qmCount] < (ns+totalSeg))) {
                currTreRow[currCount] = qmRow[qmCount];
                currTreCol[currCount] = qmCol[qmCount];
                currCount++;
            }
         }


         locate(totalSeg, start*nsinv, len*nsinv, posit+ns);

         sdoTemTime[0] = stt*spk[k]/tseg;
         sdoTemTime[1] = dtt*dpk[k]/tseg;
         sdoTemTime[2] = ott*opk[k]/tseg;

         sdoTemSite[0] = ssinton[k];
         sdoTemSite[1] = sdouton[k];
         sdoTemSite[2] = sothton[k];

         make_MFSgametesQM(nsam, posit, mhlen, seglst[seg].ptree, sdoTemTime, sdoTemSite, ns, list, currTreRow, currTreCol, currTreQM);

         ns = ns + totalSeg;
         free(seglst[seg].ptree) ;
      }

      free(spk);
      free(dpk);
      free(opk);
      free(ssinton);
      free(sdouton);
      free(sothton);


	for(i=0;i<nsam;i++) list[i][ns] = '\0' ;

	return( ns ) ;
}

/************ make_gametes.c  *******************************************
*
*
*****************************************************************************/

#define STATE0 '0'
#define STATE1 '1'
#define STATE2 '2'
#define STAR   '-'
#define QUMA   '?'

int
make_MFSgametesQM(int nsam, double *posit, long int mhlen, struct node *ptree, double* sdoTemTime, int* sdoTemSite, int ns, char **list, int* currTreRow, int* currTreCol, int currTreLen)
{
	int  tip, j, i, node, ttsites, randomN, lastNode;
	double y, ran1(), randN;
    int pickb_StonQM(int nsam, struct node *ptree, double tt);
    int pickb_DtonQM(int nsam, struct node *ptree, double tt);
    int pickb_OtonQM(int nsam, struct node *ptree, double tt);
    double sintontt, doutontt, othtontt;
    int stonnewsites, dtonnewsites, otonnewsites;
    int currSiteNsam = 0;
    int k=0, cc=0, ccnn=0;
    double st = 0.0, dt = 0.0, ot = 0.0;
    int nc = 0, ancest = 0;

    sintontt = sdoTemTime[0];
    doutontt = sdoTemTime[1];
    othtontt = sdoTemTime[2];

    stonnewsites = sdoTemSite[0];
    dtonnewsites = sdoTemSite[1];
    otonnewsites = sdoTemSite[2];

	// for( tip=0; tip < 2*nsam-1 ; tip++) {
 //   printf("node: %d, parent: %d, time: %f, nc: %d\n", tip, (ptree+tip)->abv, (ptree + (ptree+tip)->abv )->time - (ptree+tip)->time, (ptree+tip)->ncs);
 //  }

    lastNode = -1;

    j = ns;
    ttsites = stonnewsites+dtonnewsites+otonnewsites;
    while(ttsites>0) {
        st = 0.0;
        dt = 0.0;
        ot = 0.0;

        currSiteNsam = 0;
        for(i=0;i<currTreLen;i++) {
            if(currTreCol[i] == j) {
                currSiteNsam++;
            }
        }

        if(currSiteNsam > 0) { // there are some question marks for this site
            int currSams[currSiteNsam];

            currSiteNsam = 0;
            for(i=0;i<currTreLen;i++) {
                if(currTreCol[i] == j) {
                    currSams[currSiteNsam++] = currTreRow[i];
                }
            }

            //update ncs with question marks
            for( i=0; i<2*nsam-2; i++) (ptree+i)->ncs = 0;
            for( i=0; i<2*nsam-2; i++) (ptree+i)->masked = 0; // very important


            for( tip=0; tip < nsam; tip++) {
                for(i=0;i<currSiteNsam;i++) {
                    if(currSams[i] == tip) {
                        (ptree+tip)->masked = 1; // the leaf is masked
                      //  (ptree+(ptree+tip)->abv)->masked = 1; // the direct parent of masked leaf is masked
                        break;
                    }
                }

                if(!((ptree+tip)->masked)) { // update ncs
                    k = (ptree+tip)->abv;
                    while(k < 2*nsam-2) {
                        cc = (ptree+k) -> ncs;
                        (ptree+k) -> ncs = cc+1;
                        k = (ptree+k)->abv;
                    }
                }
            }

	// for( tip=0; tip < 2*nsam-1 ; tip++) {
   // printf("node: %d, parent: %d, time: %f, nc: %d, masked: %d\n", tip, (ptree+tip)->abv, (ptree + (ptree+tip)->abv )->time - (ptree+tip)->time, (ptree+tip)->ncs, (ptree+tip)->masked);
  // }
            for( i=0; i<2*nsam-2; i++) {
                if((ptree+i)->masked == 1) {
                    continue;
                }

                nc = (ptree+i)->ncs;

                if(pars.mfs.unfolded != 1) { // folded
                    if(nc == 0) {
                        st += (ptree + (ptree+i)->abv )->time - (ptree+i)->time; // singleton
                        (ptree+i)->masked = 1;
                        ancest = (ptree+i)->abv;
                        while(ancest < 2*nsam-2) {
                            ccnn = (ptree+ancest) -> ncs;
                            if(ccnn > 1) {
                                break;
                            }
                            if((ccnn == 1) && ((ptree+ancest)->masked == 0)) { //if(ccnn == 1) {
                                st += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time; // singleton
                                (ptree+ancest)->masked = 1;
                            }
                            ancest = (ptree+ancest)->abv;
                        }
                    } else if(nc == (nsam-1)) {
                        st += (ptree + (ptree+i)->abv )->time - (ptree+i)->time; // singleton
                        (ptree+i)->masked = 1;
                        ancest = (ptree+i)->abv;
                        while(ancest < 2*nsam-2) {
                            ccnn = (ptree+ancest) -> ncs;
                            if(ccnn > (nsam-1)) {
                                break;
                            }
                            if((ccnn == (nsam-1)) && ((ptree+ancest)->masked == 0)) { // if(ccnn == (nsam-1)) {
                                st += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time; // singleton
                                (ptree+ancest)->masked = 1;
                            }
                            ancest = (ptree+ancest)->abv;
                        }
                    } else if(nc == 2) {
                        dt += (ptree + (ptree+i)->abv )->time - (ptree+i)->time; // doubleton
                        (ptree+i)->masked = 1;
                        ancest = (ptree+i)->abv;
                        while(ancest < 2*nsam-2) {
                            ccnn = (ptree+ancest) -> ncs;
                            if(ccnn > 2) {
                                break;
                            }
                            if((ccnn == 2) && ((ptree+ancest)->masked == 0)) {   // if(ccnn == 2) {
                                dt += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time; // doubleton
                                (ptree+ancest)->masked = 1;
                            }
                            ancest = (ptree+ancest)->abv;
                        }
                    } else if(nc == (nsam-2)) {
                        dt += (ptree + (ptree+i)->abv )->time - (ptree+i)->time; // doubleton
                        (ptree+i)->masked = 1;
                        ancest = (ptree+i)->abv;
                        while(ancest < 2*nsam-2) {
                            ccnn = (ptree+ancest) -> ncs;
                            if(ccnn > (nsam-2)) {
                                break;
                            }
                            if((ccnn == (nsam-2)) && ((ptree+ancest)->masked == 0)) {  //if(ccnn == (nsam-2)) {
                                dt += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time; // doubleton
                                (ptree+ancest)->masked = 1;
                            }
                            ancest = (ptree+ancest)->abv;
                        }
                    } else {
                        ot += (ptree + (ptree+i)->abv )->time - (ptree+i)->time; // xton
                        (ptree+i)->masked = 1;
                        ancest = (ptree+i)->abv;
                        while(ancest < 2*nsam-2) {
                            ccnn = (ptree+ancest) -> ncs;
                            if((ccnn == (nsam-1)) || (ccnn == (nsam-2))) {
                                break;
                            }
                            if((ccnn >= 3) && (ccnn != (nsam-1)) && (ccnn != (nsam-2)) && ((ptree+ancest)->masked == 0)) {  // if((ccnn >= 3) && (ccnn != (nsam-1)) && (ccnn != (nsam-2))) {
                                ot += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time; // xton
                                (ptree+ancest)->masked = 1;
                            }
                            ancest = (ptree+ancest)->abv;
                        }
                    }
                }
            }

            for( i=0; i<2*nsam-2; i++) (ptree+i)->masked = 0; // very important
            for( tip=0; tip < nsam; tip++) {
                for(i=0;i<currSiteNsam;i++) {
                    if(currSams[i] == tip) {
                        (ptree+tip)->masked = 1; // the leaf is masked
                      //  (ptree+(ptree+tip)->abv)->masked = 1; // the direct parent of masked leaf is masked
                        break;
                    }
                }
            }

        } else {
            //update ncs with question marks
            for( i=0; i<2*nsam-2; i++) (ptree+i)->ncs = 0;
            for( i=0; i<2*nsam-2; i++) (ptree+i)->masked = 0; // very important

            for( tip=0; tip < nsam; tip++) {
                k = (ptree+tip)->abv;
                while(k < 2*nsam-2) {
                    cc = (ptree+k) -> ncs;
                    (ptree+k) -> ncs = cc+1;
                    k = (ptree+k)->abv;
                }
            }

            st = sintontt;
            dt = doutontt;
            ot = othtontt;
        }

//	 for( tip=0; tip < 2*nsam-1 ; tip++) {
 //   printf("node: %d, parent: %d, time: %f, nc: %d, masked: %d\n", tip, (ptree+tip)->abv, (ptree + (ptree+tip)->abv )->time - (ptree+tip)->time, (ptree+tip)->ncs, (ptree+tip)->masked);
 //  }
        randomN = (int)(1+ran1()*ttsites);
        if(randomN <= stonnewsites) {
            stonnewsites--;
            node = pickb_StonQM(nsam, ptree, st);
        } else if(randomN <= stonnewsites+dtonnewsites) {
            dtonnewsites--;
            node = pickb_DtonQM(nsam, ptree, dt);
        } else {
            otonnewsites--;
            node = pickb_OtonQM(nsam, ptree, ot);
        }

//printf("node: %d", 1111111);
        // when have question marks, there are not multiple hit
        for(tip=0; tip < nsam ; tip++) {
            if((ptree+tip)->masked == 1) {
                list[tip][j] = QUMA;
            } else if( tdesn(ptree, tip, node) ) {
                list[tip][j] = STATE1 ;
            } else {
                list[tip][j] = STATE0 ;
            }
        }

        j++;
        ttsites = stonnewsites+dtonnewsites+otonnewsites;
    }
}
	void
ndes_setup(struct node *ptree, int nsam )
{
	int i ;

	for( i=0; i<nsam; i++) (ptree+i)->ndes = 1 ;
	for( i = nsam; i< 2*nsam -1; i++) (ptree+i)->ndes = 0 ;
	for( i= 0; i< 2*nsam -2 ; i++)  (ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes ;

}

	int
biggerlist(nsam,  list )
	int nsam ;
	char ** list ;
{
	int i;

/*  fprintf(stderr,"maxsites: %d\n",maxsites);  */
	for( i=0; i<nsam; i++){
	   list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
	   if( list[i] == NULL ) perror( "realloc error. bigger");
	   }
}



/* allocates space for gametes (character strings) */
	char **
cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned) nsam*sizeof( char* ) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) len*sizeof( char ) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}



	int
locate(n,beg,len,ptr)
	int n;
	double beg, len, *ptr;
{
	int i;

	ordran(n,ptr);
	for(i=0; i<n; i++)
		ptr[i] = beg + ptr[i]*len ;

}

int NSEEDS = 3 ;

  void
getpars(int argc, char **argv, int *phowmany )
{
	int arg, i, j, sum , pop , argstart, npop , npop2, pop2 ;
	double migr, mij, psize, palpha ;
	void addtoelist( struct devent *pt, struct devent *elist );
	void argcheck( int arg, int argc, char ** ) ;
	int commandlineseed( char ** ) ;
	void free_eventlist( struct devent *pt, int npop );
	struct devent *ptemp , *pt ;
	FILE *pf ;
	char ch3 ;

//printf("in getPars argc:%d\n", argc);
  if( count == 0 ) {
	if( argc < 4 ){ fprintf(stderr,"Too few command line arguments\n"); usage();}
	pars.cp.nsam = atoi( argv[1] );
//	printf("in getPars  nsam:%d\n", pars.cp.nsam);
	if( pars.cp.nsam <= 0 ) { fprintf(stderr,"First argument error. nsam <= 0. \n"); usage();}
	*phowmany = atoi( argv[2] );
//	printf("in getPars  rep:%d\n", atoi( argv[2] ));
	if( *phowmany  <= 0 ) { fprintf(stderr,"Second argument error. howmany <= 0. \n"); usage();}
	pars.commandlineseedflag = 0 ;
	  pars.output_precision = 4 ;
	pars.cp.r = pars.mp.theta =  pars.cp.f = 0.0 ;
	pars.cp.track_len = 0. ;
	pars.cp.npop = npop = 1 ;
	pars.cp.mig_mat = (double **)malloc( (unsigned) sizeof( double *) );
	pars.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double ));
	pars.cp.mig_mat[0][0] =  0.0 ;
	pars.mp.segsitesin = 0 ;

	//init added arguments
	pars.mfs.unfolded = 1;
    pars.mfs.multipleHit = 0;
    pars.mfs.mhlength = 0;
    pars.mfs.moreLength = 0;
    pars.mfs.singlton = 0;
    pars.mfs.doubleton = 0;
    pars.mfs.otherton = 0;

	pars.mp.treeflag = 0 ;
 	pars.mp.timeflag = 0 ;
       pars.mp.mfreq = 1 ;
	pars.cp.config = (int *) malloc( (unsigned)(( pars.cp.npop +1 ) *sizeof( int)) );
	(pars.cp.config)[0] = pars.cp.nsam ;
	pars.cp.size= (double *) malloc( (unsigned)( pars.cp.npop *sizeof( double )) );
	(pars.cp.size)[0] = 1.0  ;
	pars.cp.alphag = (double *) malloc( (unsigned)(( pars.cp.npop ) *sizeof( double )) );
	(pars.cp.alphag)[0] = 0.0  ;
	pars.cp.nsites = 2 ;
  }
  else{
	npop = pars.cp.npop ;
	free_eventlist( pars.cp.deventlist, npop );
  }
  	pars.cp.deventlist = NULL ;

	arg = 3 ;

	while( arg < argc ){
		if( argv[arg][0] != '-' ) { fprintf(stderr," argument should be -%s ?\n", argv[arg]); usage();}
		switch ( argv[arg][1] ){
			case 'f' :
				if( ntbs > 0 ) { fprintf(stderr," can't use tbs args and -f option.\n"); exit(1); }
				arg++;
				argcheck( arg, argc, argv);
				pf = fopen( argv[arg], "r" ) ;
				if( pf == NULL ) {fprintf(stderr," no parameter file %s\n", argv[arg] ); exit(0);}
				arg++;
				argc++ ;
				argv = (char **)malloc(  (unsigned)(argc+1)*sizeof( char *) ) ;
				argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
				argstart = arg ;
				while( fscanf(pf," %s", argv[arg]) != EOF ) {
					arg++;
					argc++;
					argv = (char **)realloc( argv, (unsigned)argc*sizeof( char*) ) ;
				        argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
					}
				fclose(pf);
				argc--;
				arg = argstart ;
				break;
			case 'r' :
				arg++;
				argcheck( arg, argc, argv);
				pars.cp.r = atof(  argv[arg++] );
				argcheck( arg, argc, argv);
				pars.cp.nsites = atoi( argv[arg++]);
				if( pars.cp.nsites <2 ){
					fprintf(stderr,"with -r option must specify both rec_rate and nsites>1\n");
					usage();
					}
				break;
			case 'p' :
				arg++;
				argcheck(arg,argc,argv);
				pars.output_precision = atoi( argv[arg++] ) ;
				break;
			case 'c' :
				arg++;
				argcheck( arg, argc, argv);
				pars.cp.f = atof(  argv[arg++] );
				argcheck( arg, argc, argv);
				pars.cp.track_len = atof( argv[arg++]);
				if( pars.cp.track_len <1. ){
					fprintf(stderr,"with -c option must specify both f and track_len>0\n");
					usage();
					}
				break;
            case 'Y' :
                arg++;
                pars.mfs.singlton = atoi( argv[arg++]);
                pars.mfs.doubleton = atoi( argv[arg++]);
                pars.mfs.otherton = atoi( argv[arg++]);
                break;
            case 'U' :
                arg++;
                pars.mfs.unfolded = atoi( argv[arg++]);
                break;
            case 'H' :
                arg++;
                pars.mfs.multipleHit = atoi( argv[arg++]);
                pars.mfs.mhlength = atoi( argv[arg++]);
                break;
            case 'A' :
                arg++;
                pars.mfs.moreLength = atoi( argv[arg++]);
                break;
			case 't' :
				arg++;
				argcheck( arg, argc, argv);
				pars.mp.theta = atof(  argv[arg++] );
				break;
			case 's' :
				arg++;
				argcheck( arg, argc, argv);
				if( argv[arg-1][2] == 'e' ){  /* command line seeds */
					pars.commandlineseedflag = 1 ;
					if( count == 0 ) nseeds = commandlineseed(argv+arg );
					arg += nseeds ;
				}
				else {
				    pars.mp.segsitesin = atoi(  argv[arg++] );
				}
				break;
			case 'F' :
				arg++;
				argcheck( arg, argc, argv);
				pars.mp.mfreq = atoi(  argv[arg++] );
                                if( (pars.mp.mfreq < 2 ) || (pars.mp.mfreq > pars.cp.nsam/2 ) ){
                                    fprintf(stderr," mfreq must be >= 2 and <= nsam/2.\n");
                                    usage();
                                    }
				break;
			case 'T' :
				pars.mp.treeflag = 1 ;
				arg++;
				break;
			case 'L' :
				pars.mp.timeflag = 1 ;
				arg++;
				break;
			case 'I' :
			    arg++;
			    if( count == 0 ) {
				argcheck( arg, argc, argv);
			       	pars.cp.npop = atoi( argv[arg]);
			        pars.cp.config = (int *) realloc( pars.cp.config, (unsigned)( pars.cp.npop*sizeof( int)));
				npop = pars.cp.npop ;
				}
			    arg++;
			    for( i=0; i< pars.cp.npop; i++) {
				argcheck( arg, argc, argv);
				pars.cp.config[i] = atoi( argv[arg++]);
				}
			    if( count == 0 ){
				pars.cp.mig_mat =
                                        (double **)realloc(pars.cp.mig_mat, (unsigned)(pars.cp.npop*sizeof(double *) )) ;
				pars.cp.mig_mat[0] =
                                         (double *)realloc(pars.cp.mig_mat[0], (unsigned)( pars.cp.npop*sizeof(double)));
				for(i=1; i<pars.cp.npop; i++) pars.cp.mig_mat[i] =
                                         (double *)malloc( (unsigned)( pars.cp.npop*sizeof(double)));
				pars.cp.size = (double *)realloc( pars.cp.size, (unsigned)( pars.cp.npop*sizeof( double )));
				pars.cp.alphag =
                                          (double *) realloc( pars.cp.alphag, (unsigned)( pars.cp.npop*sizeof( double )));
			        for( i=1; i< pars.cp.npop ; i++) {
				   (pars.cp.size)[i] = (pars.cp.size)[0]  ;
				   (pars.cp.alphag)[i] = (pars.cp.alphag)[0] ;
				   }
			        }
			     if( (arg <argc) && ( argv[arg][0] != '-' ) ) {
				argcheck( arg, argc, argv);
				migr = atof(  argv[arg++] );
				}
			     else migr = 0.0 ;
			     for( i=0; i<pars.cp.npop; i++)
				    for( j=0; j<pars.cp.npop; j++) pars.cp.mig_mat[i][j] = migr/(pars.cp.npop-1) ;
			     for( i=0; i< pars.cp.npop; i++) pars.cp.mig_mat[i][i] = migr ;
			     break;
			case 'm' :
			     if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
			     if( argv[arg][2] == 'a' ) {
				    arg++;
				    for( pop = 0; pop <npop; pop++)
				      for( pop2 = 0; pop2 <npop; pop2++){
					     argcheck( arg, argc, argv);
					     pars.cp.mig_mat[pop][pop2]= atof( argv[arg++] ) ;
					  }
				    for( pop = 0; pop < npop; pop++) {
					  pars.cp.mig_mat[pop][pop] = 0.0 ;
					  for( pop2 = 0; pop2 < npop; pop2++){
					    if( pop2 != pop ) pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2] ;
					  }
				    }
				}
			    else {
		             arg++;
			         argcheck( arg, argc, argv);
		             i = atoi( argv[arg++] ) -1;
			         argcheck( arg, argc, argv);
		             j = atoi( argv[arg++] ) -1;
			         argcheck( arg, argc, argv);
		             mij = atof( argv[arg++] );
		             pars.cp.mig_mat[i][i] += mij -  pars.cp.mig_mat[i][j]  ;
		             pars.cp.mig_mat[i][j] = mij;
			    }
				break;
			case 'n' :
			     if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
			    arg++;
			    argcheck( arg, argc, argv);
			    pop = atoi( argv[arg++] ) -1;
			    argcheck( arg, argc, argv);
			    psize = atof( argv[arg++] );
			    pars.cp.size[pop] = psize ;
			   break;
			case 'g' :
			     if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
			    arg++;
			    argcheck( arg, argc, argv);
			    pop = atoi( argv[arg++] ) -1;
			    if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); usage(); }
			    palpha = atof( argv[arg++] );
			    pars.cp.alphag[pop] = palpha ;
			   break;
			case 'G' :
			    arg++;
			    if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); usage(); }
			    palpha = atof( argv[arg++] );
			    for( i=0; i<pars.cp.npop; i++)
			       pars.cp.alphag[i] = palpha ;
			   break;
			case 'e' :
			    pt = (struct devent *)malloc( sizeof( struct devent) ) ;
			    pt->detype = argv[arg][2] ;
			    ch3 = argv[arg][3] ;
			    arg++;
			    argcheck( arg, argc, argv);
			    pt->time = atof( argv[arg++] ) ;
			    pt->nextde = NULL ;
			    if( pars.cp.deventlist == NULL )
				    pars.cp.deventlist = pt ;
			    else if ( pt->time < pars.cp.deventlist->time ) {
				    ptemp = pars.cp.deventlist ;
				    pars.cp.deventlist = pt ;
				    pt->nextde = ptemp ;
				}
			    else
				   addtoelist( pt, pars.cp.deventlist ) ;
			    switch( pt->detype ) {
				case 'N' :
			          argcheck( arg, argc, argv);
				      pt->paramv = atof( argv[arg++] ) ;
				      break;
				case 'G' :
				  if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eG.\n"); usage(); }
				  pt->paramv = atof( argv[arg++] ) ;
				  break;
				case 'M' :
				    argcheck( arg, argc, argv);
				    pt->paramv = atof( argv[arg++] ) ;
				    break;
				case 'n' :
			          argcheck( arg, argc, argv);
				  pt->popi = atoi( argv[arg++] ) -1 ;
			          argcheck( arg, argc, argv);
				  pt->paramv = atof( argv[arg++] ) ;
				  break;
				case 'g' :
			          argcheck( arg, argc, argv);
				  pt->popi = atoi( argv[arg++] ) -1 ;
				  if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eg.\n"); usage(); }
				  pt->paramv = atof( argv[arg++] ) ;
				  break;
				case 's' :
			          argcheck( arg, argc, argv);
				  pt->popi = atoi( argv[arg++] ) -1 ;
			          argcheck( arg, argc, argv);
				  pt->paramv = atof( argv[arg++] ) ;
				  break;
				case 'm' :
				  if( ch3 == 'a' ) {
				     pt->detype = 'a' ;
				     argcheck( arg, argc, argv);
				     npop2 = atoi( argv[arg++] ) ;
				     pt->mat = (double **)malloc( (unsigned)npop2*sizeof( double *) ) ;
				     for( pop =0; pop <npop2; pop++){
					   (pt->mat)[pop] = (double *)malloc( (unsigned)npop2*sizeof( double) );
					   for( i=0; i<npop2; i++){
					     if( i == pop ) arg++;
					     else {
				               argcheck( arg, argc, argv);
					       (pt->mat)[pop][i] = atof( argv[arg++] ) ;
					     }
					   }
				     }
				     for( pop = 0; pop < npop2; pop++) {
					    (pt->mat)[pop][pop] = 0.0 ;
					    for( pop2 = 0; pop2 < npop2; pop2++){
					       if( pop2 != pop ) (pt->mat)[pop][pop] += (pt->mat)[pop][pop2] ;
					    }
				     }
				  }
				  else {
			            argcheck( arg, argc, argv);
				        pt->popi = atoi( argv[arg++] ) -1 ;
			            argcheck( arg, argc, argv);
				        pt->popj = atoi( argv[arg++] ) -1 ;
			            argcheck( arg, argc, argv);
				        pt->paramv = atof( argv[arg++] ) ;
				  }
				  break;
				case 'j' :
			          argcheck( arg, argc, argv);
				  pt->popi = atoi( argv[arg++] ) -1 ;
			          argcheck( arg, argc, argv);
				  pt->popj = atoi( argv[arg++] ) -1 ;
				  break;
				default: fprintf(stderr,"e event\n");  usage();
			    }
			 break;
			default: fprintf(stderr," option default\n");  usage() ;
			}
		}
		if( (pars.mp.theta == 0.0) && ( pars.mp.segsitesin == 0 ) && ( pars.mp.treeflag == 0 ) && (pars.mp.timeflag == 0) && (pars.mfs.singlton == 0) && (pars.mfs.doubleton == 0) && (pars.mfs.otherton == 0) ) {
			fprintf(stderr," either -s or -t or -Y or -T option must be used. \n");
			usage();
			exit(1);
        }
		sum = 0 ;
		for( i=0; i< pars.cp.npop; i++) sum += (pars.cp.config)[i] ;
		if( sum != pars.cp.nsam ) {
			fprintf(stderr," sum sample sizes != nsam\n");
			usage();
			exit(1);
        }

        if((pars.mp.segsitesin > 0) && (pars.mfs.unfolded != 1)){
            fprintf(stderr,"\nCannot set -s and -U simultaneous.\n");
            usage();
			exit(1);
        }

        if((pars.mfs.unfolded != 1) && (pars.cp.nsam < 6)) {
            fprintf(stderr,"\nWith folded MFS, sample size must be equal or greater than 6. \n");
            usage();
			exit(1);
        }

        if((pars.mfs.multipleHit == 0) && (pars.mfs.moreLength > 0)) {
            fprintf(stderr,"\nWithout multiple hit , you must not specify -A moreLength. \n");
            usage();
			exit(1);
        }

       if((pars.mfs.moreLength > 0) && ((pars.mfs.singlton == 0) && (pars.mfs.doubleton == 0) && (pars.mfs.otherton == 0))){
            fprintf(stderr,"\nWithout -Y, you must not specify -A moreLength. \n");
            usage();
			exit(1);
        }

        if((pars.mp.segsitesin > 0) && ((pars.mfs.singlton > 0) || (pars.mfs.doubleton > 0) || (pars.mfs.otherton > 0))){
            fprintf(stderr,"\nCannot set -s and -Y simultaneous.\n");
            usage();
			exit(1);
        }
}


	void
argcheck( int arg, int argc, char *argv[] )
{
	if( (arg >= argc ) || ( argv[arg][0] == '-') ) {
	   fprintf(stderr,"not enough arguments after %s\n", argv[arg-1] ) ;
	   fprintf(stderr,"For usage type: ms<return>\n");
	   exit(0);
	  }
}

	int
usage()
{
fprintf(stderr,"usage: ms nsam howmany \n");
fprintf(stderr,"  Options: \n");
fprintf(stderr,"\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
fprintf(stderr,"\t -s segsites   (fixed number of segregating sites)\n");
fprintf(stderr,"\t -Y singleton doubleton tripleton  (fixed number of MFS)\n");
fprintf(stderr,"\t -U 1/0 (1 represents unfolded MFS; 0 represents folded MFS. Default is 1.)\n");
fprintf(stderr,"\t -H 0/1 mhlength (0 represents no multiple hit; 1 represents multiple hit. Default is 0. If you specify 1, you must specify 'length'.)\n");
fprintf(stderr,"\t -A morelength (If you specify multiple hit and you want to specify more length.)\n");
fprintf(stderr,"\t -T          (Output gene tree.)\n");
fprintf(stderr,"\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
fprintf(stderr,"\t -r rho nsites     (rho here is 4Nc)\n");
fprintf(stderr,"\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
fprintf(stderr,"\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n");
fprintf(stderr,"\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");
fprintf(stderr,"\t -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");
fprintf(stderr,"\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n");
fprintf(stderr,"\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");
fprintf(stderr,"\t\t -n i size_i   (popi has size set to size_i*N0 \n");
fprintf(stderr,"\t\t -g i alpha_i  (If used must appear after -M option.)\n");
fprintf(stderr,"\t   The following options modify parameters at the time 't' specified as the first argument:\n");
fprintf(stderr,"\t -eG t alpha  (Modify growth rate of all pop's.)\n");
fprintf(stderr,"\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");
fprintf(stderr,"\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n");
fprintf(stderr,"\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n");
fprintf(stderr,"\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");
fprintf(stderr,"\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");
fprintf(stderr,"\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
fprintf(stderr,"\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");
fprintf(stderr,"\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
fprintf(stderr,"\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
fprintf(stderr,"\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
fprintf(stderr,"\t\t  size, alpha and M are unchanged.\n");
fprintf(stderr,"\t  -f filename     ( Read command line arguments from file filename.)\n");
fprintf(stderr,"\t  -p n ( Specifies the precision of the position output.  n is the number of digits after the decimal.)\n");
fprintf(stderr," See msdoc.pdf for explanation of these parameters.\n");

exit(1);
}
	void
addtoelist( struct devent *pt, struct devent *elist )
{
	struct devent *plast, *pevent, *ptemp  ;

	pevent = elist ;
	while(  (pevent != NULL ) && ( pevent->time <= pt->time ) )  {
		plast = pevent ;
		pevent = pevent->nextde ;
		}
	ptemp = plast->nextde ;
	plast->nextde = pt ;
	pt->nextde = ptemp ;
}

	void
free_eventlist( struct devent *pt, int npop )
{
   struct devent *next ;
   int pop ;

   while( pt != NULL){
	  next = pt->nextde ;
	  if( pt->detype == 'a' ) {
	     for( pop = 0; pop < npop; pop++) free( (pt->mat)[pop] );
		 free( pt->mat );
	  }
	  free(pt);
	  pt = next ;
   }
}




	int
make_gametes(int nsam, int mfreq, double *posit, long int nsites, struct node *ptree, double tt, int newsites, int ns, char **list )
{
	int  tip, j,  node, lastNode;
        int pickb(int nsam, struct node *ptree, double tt),
            pickbmf(int nsam, int mfreq, struct node *ptree, double tt) ;
    double ran1(), randN;
    lastNode = -1;

	//for( tip=0; tip < 2*nsam-1 ; tip++) {
   //   printf("node: %d, parent: %d, time: %f\n", tip, (ptree+tip)->abv, (ptree + (ptree+tip)->abv )->time - (ptree+tip)->time);
  // }

	for(  j=ns; j< ns+newsites ;  j++ ) {
		if( mfreq == 1 ) node = pickb(  nsam, ptree, tt);
		else node = pickbmf(  nsam, mfreq, ptree, tt);

		if(pars.mfs.multipleHit == 1) { //have multiple hit
            if( j>ns ) {
                if(lastNode == -1) {
                    for( tip=0; tip < nsam ; tip++) {
                        if( tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
                        else list[tip][j] = STATE0 ;
                    }
                    lastNode = node;
                } else {
                    if(floor(posit[j]*nsites) == floor(posit[j-1]*nsites)) {
                        if(lastNode == node) { // ??һ????????��??node??ͬ
                        } else if(tdesn(ptree, lastNode, node) || tdesn(ptree, node, lastNode)) { // ?ڶ?????????��??node?д?????ϵ
                            for( tip=0; tip < nsam ; tip++) {
                                if( tdesn(ptree, tip, lastNode) && tdesn(ptree, tip, node)) {
                                    list[tip][j-1] = STATE2;
                                } else if(tdesn(ptree, tip, lastNode) || tdesn(ptree, tip, node)) {
                                    list[tip][j-1] = STATE1;
                                } else {
                                    list[tip][j-1] = STATE0;
                                }
                            }
                        } else { // ????????????��??nodeû?д?????ϵ
                            randN = ran1();
                            for( tip=0; tip < nsam ; tip++) {
                                if( tdesn(ptree, tip, node) ) {
                                     if(randN<1.0/3.0) {
                                        list[tip][j-1] = STATE1;
                                     } else {
                                        list[tip][j-1] = STATE2;
                                     }
                                }
                            }
                        }
                        for( tip=0; tip < nsam ; tip++) {
                            list[tip][j] = STAR;
                        }
                        lastNode = -1;
                    } else {
                        for( tip=0; tip < nsam ; tip++) {
                            if( tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
                            else list[tip][j] = STATE0 ;
                        }
                        lastNode = node;
                    }
                }
            } else {
                for( tip=0; tip < nsam ; tip++) {
                    if( tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
                    else list[tip][j] = STATE0 ;
                }
                lastNode = node;
            }
		} else {
            for( tip=0; tip < nsam ; tip++) {
               if( tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
               else list[tip][j] = STATE0 ;
            }
        }
    }
}

int
make_MFSgametes(int nsam, int mfreq, double *posit, long int mhlen, struct node *ptree, double sintontt, double doutontt, double othtontt, int stonnewsites, int dtonnewsites, int otonnewsites, int ns, char **list, int loop)
{
	int  tip, j, i, node, ttsites, randomN, lastNode;
	double y, ran1(), randN;
    int pickb_Ston(int nsam, struct node *ptree, double tt),  pickbmf_Ston(int nsam, int mfreq, struct node *ptree, double tt) ;
    int pickb_Dton(int nsam, struct node *ptree, double tt),  pickbmf_Dton(int nsam, int mfreq, struct node *ptree, double tt) ;
    int pickb_Oton(int nsam, struct node *ptree, double tt),  pickbmf_Oton(int nsam, int mfreq, struct node *ptree, double tt) ;

/*	 for( tip=0; tip < 2*nsam-1 ; tip++) {
    printf("node: %d, parent: %d, time: %f, nc: %d\n", tip, (ptree+tip)->abv, (ptree + (ptree+tip)->abv )->time - (ptree+tip)->time, (ptree+tip)->ncs);
   }
*/
    lastNode = -1;

    j = ns;
    ttsites = stonnewsites+dtonnewsites+otonnewsites;
    while(ttsites>0) {
        randomN = (int)(1+ran1()*ttsites);
        if(randomN <= stonnewsites) {
			stonnewsites--;
			if(mfreq == 1) {
                node = pickb_Ston(nsam, ptree, sintontt);
			} else node = pickbmf_Ston(nsam, mfreq, ptree, sintontt);
        } else if(randomN <= stonnewsites+dtonnewsites) {
			dtonnewsites--;
			if(mfreq == 1) {
                node = pickb_Dton(nsam, ptree, doutontt);
			} else node = pickbmf_Dton(  nsam, mfreq, ptree, doutontt);
		} else {
			otonnewsites--;
			if(mfreq == 1) {
                node = pickb_Oton(nsam, ptree, othtontt);
			} else node = pickbmf_Oton(  nsam, mfreq, ptree, othtontt);
		}

        /*if((pars.mfs.multipleHit == 1) && (pars.mfs.strNumber == 1)) { //have multiple hit
            if( j>ns ) {
                if(lastNode == -1) {
                    for( tip=0; tip < nsam ; tip++) {
                        if( tdesn(ptree, tip, node) ) list[tip+loop*nsam][j] = STATE1 ;
                        else list[tip+loop*nsam][j] = STATE0 ;
                    }
                    lastNode = node;
                } else {
                    if(floor(posit[j]*mhlen) == floor(posit[j-1]*mhlen)) {
                        if(lastNode == node) { // ??һ????????��??node??ͬ
                        } else if(tdesn(ptree, lastNode, node) || tdesn(ptree, node, lastNode)) { // ?ڶ?????????��??node?д?????ϵ
                            for( tip=0; tip < nsam ; tip++) {
                                if( tdesn(ptree, tip, lastNode) && tdesn(ptree, tip, node)) {
                                    list[tip+loop*nsam][j-1] = STATE2;
                                } else if(tdesn(ptree, tip, lastNode) || tdesn(ptree, tip, node)) {
                                    list[tip+loop*nsam][j-1] = STATE1;
                                } else {
                                    list[tip+loop*nsam][j-1] = STATE0;
                                }
                            }
                        } else { // ????????????��??nodeû?д?????ϵ
                            randN = ran1();
                            for( tip=0; tip < nsam ; tip++) {
                                if( tdesn(ptree, tip, node) ) {
                                     if(randN<1.0/3.0) {
                                        list[tip+loop*nsam][j-1] = STATE1;
                                     } else {
                                        list[tip+loop*nsam][j-1] = STATE2;
                                     }
                                }
                            }
                        }
                        for( tip=0; tip < nsam ; tip++) {
                            list[tip+loop*nsam][j] = STAR;
                        }
                        lastNode = -1;
                    } else {
                        for( tip=0; tip < nsam ; tip++) {
                            if( tdesn(ptree, tip, node) ) list[tip+loop*nsam][j] = STATE1 ;
                            else list[tip+loop*nsam][j] = STATE0 ;
                        }
                        lastNode = node;
                    }
                }

            } else {
                for( tip=0; tip < nsam ; tip++) {
                    if( tdesn(ptree, tip, node) ) list[tip+loop*nsam][j] = STATE1 ;
                    else list[tip+loop*nsam][j] = STATE0 ;
                }

                lastNode = node;
            }
        } else { *///not have multiple hit
            for(tip=0; tip < nsam ; tip++) {
                if(tdesn(ptree, tip, node) ) list[tip+loop*nsam][j] = STATE1;
                else list[tip+loop*nsam][j] = STATE0;
            }
        //}

        j++;
        ttsites = stonnewsites+dtonnewsites+otonnewsites;
    }
}



/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/

	double
ttime(ptree, nsam)
	struct node *ptree;
	int nsam;
{
	double t;
	int i;

	t = (ptree + 2*nsam-2) -> time ;
	for( i=nsam; i< 2*nsam-1 ; i++)
		t += (ptree + i)-> time ;
	return(t);
}

/***  sdottime: Returns the total time of singleton, doubleton, xton in the tree **/
void sdottime(struct node *ptree, int nsam, double* dotime, int unfo)
{
    double st, dt, ot;
	int i, tip, cc, k, nc;
	st = 0.0;
	dt = 0.0;
	ot = 0.0;


	//for( tip=0; tip < 2*nsam-1 ; tip++) {
   //   printf("node: %d, parent: %d, time: %f\n", tip, (ptree+tip)->abv, (ptree + (ptree+tip)->abv )->time - (ptree+tip)->time);
   // }
    // init number of children
    for( i=0; i<2*nsam-2; i++) (ptree+i)->ncs = 0;
    for( tip=0; tip < nsam; tip++) {
        k = (ptree+tip)->abv;
        while(k < 2*nsam-2) {
            cc = (ptree+k) -> ncs;
            (ptree+k) -> ncs = cc+1;
            k = (ptree+k)->abv;
        }
    }

//	for( tip=0; tip < 2*nsam-1 ; tip++) {
 //       printf("node: %d, ncs: %d\n", tip, (ptree+tip)->ncs);
 //   }

	for( i=0; i<2*nsam-2; i++) {
        nc = (ptree+i)->ncs;

        if(unfo != 1) { // folded
            if((nc == 0) || (nc == (nsam-1)))  {
                st += (ptree + (ptree+i)->abv )->time - (ptree+i)->time; // singleton
            } else if ((nc == 2) || (nc == (nsam-2))) {
                dt += (ptree + (ptree+i)->abv )->time - (ptree+i)->time; // doubleton
            } else {
                ot += (ptree + (ptree+i)->abv )->time - (ptree+i)->time; // xton
            }

        } else { // unfolded
            if(nc == 0) {
                st += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
            } else if (nc == 2) {
                dt += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
            } else {
                ot += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
            }
        }
    }

    *(dotime) = st;
    *(dotime+1) = dt;
    *(dotime+2) = ot;
}



	double
ttimemf( ptree, nsam, mfreq)
	struct node *ptree;
	int nsam, mfreq;
{
	double t;
	int i;

	t = 0. ;
	for( i=0;  i< 2*nsam-2  ; i++)
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	return(t);
}

/***  sdottimemf: With mfreq, Returns the total time of singleton, doubleton, xton in the tree **/
void sdottimemf(struct node *ptree, int nsam, double* dotime, int mfreq, int unfo)
{
    double st, dt, ot;
	int i, tip, cc, k, nc;
	st = 0.0;
	dt = 0.0;
	ot = 0.0;


    // init number of children
    for( i=0; i<2*nsam-2; i++) (ptree+i)->ncs = 0;
    for( tip=0; tip < nsam; tip++) {
        k = (ptree+tip)->abv;
        while(k < 2*nsam-2) {
            cc = (ptree+k) -> ncs;
            (ptree+k) -> ncs = cc+1;
            k = (ptree+k)->abv;
        }
    }


    for( i=0; i<2*nsam-2; i++) {
        nc = (ptree+i)->ncs;

        if(unfo != 1) { // folded
            if((nc == 0) || (nc == (nsam-1)))  {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ) {
                    st += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
                }
            } else if ((nc == 2) || (nc == (nsam-2))) {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ) {
                    dt += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
                }
            } else {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ) {
                    ot += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
                }
            }

        } else { // unfolded
            if(nc == 0) {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ) {
                    st += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
                }
            } else if (nc == 2) {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ) {
                    dt += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
                }
            } else {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ) {
                    ot += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
                }
            }
        }
    }

    *(dotime) = st;
    *(dotime+1) = dt;
    *(dotime+2) = ot;


}


	void
prtree( ptree, nsam)
	struct node *ptree;
	int nsam;
{
	double t;
	int i, *descl, *descr ;
	void parens( struct node *ptree, int *descl, int *descr, int noden );

	descl = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	descr = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	for( i=0; i<2*nsam-1; i++) descl[i] = descr[i] = -1 ;
	for( i = 0; i< 2*nsam-2; i++){
	  if( descl[ (ptree+i)->abv ] == -1 ) descl[(ptree+i)->abv] = i ;
	  else descr[ (ptree+i)->abv] = i ;
	 }
	parens( ptree, descl, descr, 2*nsam-2);
	free( descl ) ;
	free( descr ) ;
}

	void
parens( struct node *ptree, int *descl, int *descr,  int noden)
{
	double time ;

   if( descl[noden] == -1 ) {
	printf("%d:%5.3lf", noden+1, (ptree+ ((ptree+noden)->abv))->time );
	}
   else{
	printf("(");
	parens( ptree, descl,descr, descl[noden] ) ;
	printf(",");
	parens(ptree, descl, descr, descr[noden] ) ;
	if( (ptree+noden)->abv == 0 ) printf(");\n");
	else {
	  time = (ptree + (ptree+noden)->abv )->time - (ptree+noden)->time ;
	  printf("):%5.3lf", time );
	  }
        }
}




int
pickb_StonQM(int nsam, struct node *ptree, double tt)
{
	double x, y=0.0, ran1();
	int i, nc, ancest, ccnn;
    x = ran1()*tt;
	//if(pars.mfs.unfolded != 1) {
        for( i=0, y=0; i < nsam; i++) {
            if((ptree+i)->masked == 1) {
                continue;
            }

            nc = (ptree+i)->ncs;

            if(nc == 0) {

                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                (ptree+i)->masked = 1;
                if( y >= x ) return( i ) ;

                ancest = (ptree+i)->abv;
                while(ancest < 2*nsam-2) {
                    ccnn = (ptree+ancest) -> ncs;
                    if(ccnn > 1) {
                        break;
                    }
                    if((ccnn == 1) && ((ptree+ancest)->masked == 0)) {
                        y += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time;
                        (ptree+ancest)->masked = 1;

                        if( y >= x ) return( ancest ) ;
                    }
                    ancest = (ptree+ancest)->abv;
                }

            } else if(nc == (nsam-1)) {

                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                (ptree+i)->masked = 1;
                if( y >= x ) return( i ) ;

                ancest = (ptree+i)->abv;
                while(ancest < 2*nsam-2) {
                    ccnn = (ptree+ancest) -> ncs;
                    if(ccnn > (nsam-1)) {
                        break;
                    }
                    if((ccnn == (nsam-1)) && ((ptree+ancest)->masked == 0)) {
                        y += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time;
                        (ptree+ancest)->masked = 1;

                        if( y >= x ) return( ancest ) ;
                    }
                    ancest = (ptree+ancest)->abv;
                }

            }

        }
    /*} else {
        for( i=0, y=0; i < nsam; i++) {
            y += (ptree + (ptree+i)->abv )->time;
            if( y >= x ) return( i ) ;
        }
    }*/
	return( nsam-1 );
}


int
pickb_DtonQM(nsam, ptree, tt)
	int nsam;
	struct node *ptree;
	double tt;
{
    double x, y=0.0, ran1();
	int i, j, nc, ccnn;
	int tip, ancest;
    x = ran1()*tt;

  //  if(pars.mfs.unfolded != 1) {
        for( i=nsam, y=0; i < 2*nsam-2 ; i++) {
            if((ptree+i)->masked == 1) {
                continue;
            }

            nc = (ptree+i)->ncs;

            if(nc == 2) {

                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                (ptree+i)->masked = 1;

                if( y >= x ) return( i ) ;
                // only start from original not masked node
                ancest = (ptree+i)->abv;
                while(ancest < 2*nsam-2) {
                    ccnn = (ptree+ancest) -> ncs;
                    if(ccnn > 2) {
                        break;
                    }
                    if((ccnn == 2) && ((ptree+ancest)->masked == 0)) {   // if(ccnn == 2) {
                        y += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time;
                        (ptree+ancest)->masked = 1;

                        if( y >= x ) return( ancest ) ;
                    }
                    ancest = (ptree+ancest)->abv;
                }



            } else if(nc == (nsam-2)) {

                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                (ptree+i)->masked = 1;

                if( y >= x ) return( i ) ;
                // only start from original not masked node
                ancest = (ptree+i)->abv;
                while(ancest < 2*nsam-2) {
                    ccnn = (ptree+ancest) -> ncs;
                    if(ccnn > (nsam-2)) {
                        break;
                    }
                     if((ccnn == (nsam-2)) && ((ptree+ancest)->masked == 0)) {  //if(ccnn == (nsam-2)) {
                        y += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time;
                        (ptree+ancest)->masked = 1;

                        if( y >= x ) return( ancest ) ;
                    }
                    ancest = (ptree+ancest)->abv;
                }


            }
        }
   // } else {
    	/*for( i=nsam, y=0; i < 2*nsam-2 ; i++) {
            if((ptree+i)->ncs == 2) {
                if((ptree+i)->masked == 0) {
                    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                    if( y >= x ) return( i ) ;
                    // only start from original not masked node
                    ancest = (ptree+i)->abv;
                    while(ancest < 2*nsam-2) {
                        ccnn = (ptree+ancest) -> ncs;
                        if(ccnn > 2) {
                            break;
                        }
                        if((ccnn == 2) && ((ptree+ancest)->masked == 0)) {
                            y += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time;
                            if( y >= x ) return( ancest ) ;
                        }
                        ancest = (ptree+ancest)->abv;
                    }
                }

            }
        }*/
   // }
	//printf("Warning: Have selected the LastNode: %d\n", lastNode);
	//return( lastNode );
}

int
pickb_OtonQM(nsam, ptree, tt)
	int nsam;
	struct node *ptree;
	double tt;
{
    double x, y=0.0, ran1();
	int i, j, nc, ccnn;
	int tip, ancest;
    x = ran1()*tt;

   // if(pars.mfs.unfolded != 1) {
         for( i=nsam, y=0; i < 2*nsam-2; i++) {
            if((ptree+i)->masked == 1) {
                continue;
            }
            nc = (ptree+i)->ncs;

            if((nc >= 3) && (nc != (nsam-1)) && (nc != (nsam-2))) {
                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                (ptree+i)->masked = 1;

                if( y >= x ) return( i ) ;
                // only start from original not masked node
                ancest = (ptree+i)->abv;
                while(ancest < 2*nsam-2) {
                    ccnn = (ptree+ancest) -> ncs;
                    if((ccnn == (nsam-1)) || (ccnn == (nsam-2))) {
                        break;
                    }
                    if((ccnn >= 3) && (ccnn != (nsam-1)) && (ccnn != (nsam-2)) && ((ptree+ancest)->masked == 0)) {  // if((ccnn >= 3) && (ccnn != (nsam-1)) && (ccnn != (nsam-2))) {
                        y += (ptree + (ptree+ancest)->abv )->time - (ptree+ancest)->time;
                        (ptree+ancest)->masked = 1;

                        if( y >= x ) return( ancest ) ;
                    }
                    ancest = (ptree+ancest)->abv;
                }
            }

        }
   /* } else {
        for( i=nsam, y=0; i < 2*nsam-2; i++) {
            if((ptree+i)->ncs >=3 ) {
                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                if( y >= x ) return( i ) ;
            }
        }
    } */


	//	printf("Warning: Have selected the node : %d\n", 2*nsam-3);
	return( 2*nsam-3 );
}


int
pickb_Ston(nsam, ptree, tt)
	int nsam;
	struct node *ptree;
	double tt;
{
	double x, y=0.0, ran1();
	int i, nc;
    x = ran1()*tt;
	if(pars.mfs.unfolded != 1) {
        for( i=0, y=0; i < 2*nsam-2 ; i++) {
            nc = (ptree+i)->ncs;
            if((nc == 0) || (nc == (nsam-1))) {
                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                if( y >= x ) return( i ) ;
            }
        }
    } else {
        for( i=0, y=0; i < nsam; i++) {
            y += (ptree + (ptree+i)->abv )->time;
            if( y >= x ) return( i ) ;
        }
    }
	return( nsam-1 );
}


	int
pickb_Dton(nsam, ptree, tt)
	int nsam;
	struct node *ptree;
	double tt;
{
    double x, y=0.0, ran1();
	int i, j, nc;
	int tip;
    x = ran1()*tt;

    if(pars.mfs.unfolded != 1) {
        for( i=nsam, y=0; i < 2*nsam-2 ; i++) {
            nc = (ptree+i)->ncs;
            if((nc == 2) || (nc == (nsam-2))) {
                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                if( y >= x ) return( i ) ;
            }
        }
    } else {
    	for( i=nsam, y=0; i < 2*nsam-2 ; i++) {
            if((ptree+i)->ncs == 2) {
                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                if( y >= x ) return( i ) ;
            }
        }
    }
	//printf("Warning: Have selected the LastNode: %d\n", lastNode);
	//return( lastNode );
}

int
pickb_Oton(nsam, ptree, tt)
	int nsam;
	struct node *ptree;
	double tt;
{
    double x, y=0.0, ran1();
	int i, j, nc;
	int tip;
    x = ran1()*tt;

    if(pars.mfs.unfolded != 1) {
         for( i=nsam, y=0; i < 2*nsam-2; i++) {
            nc = (ptree+i)->ncs;
            if((nc >= 3) && (nc != (nsam-1)) && (nc != (nsam-2))) {
                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                if( y >= x ) return( i ) ;
            }
        }
    } else {
        for( i=nsam, y=0; i < 2*nsam-2; i++) {
            if((ptree+i)->ncs >=3 ) {
                y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                if( y >= x ) return( i ) ;
            }
        }
    }


	//	printf("Warning: Have selected the node : %d\n", 2*nsam-3);
	return( 2*nsam-3 );
}


/***  pickb : returns a random branch from the tree. The probability of picking
              a particular branch is proportional to its duration. tt is total
	      time in tree.   ****/

	int
pickb(nsam, ptree, tt)
	int nsam;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
		if( y >= x ) return( i ) ;
		}
	return( 2*nsam - 3  );  /* changed 4 Feb 2010 */
}


	int
pickbmf(nsam, mfreq, ptree, tt )
	int nsam, mfreq;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i, lastbranch = 0 ;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ){
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
		lastbranch = i ;    /* changed 4 Feb 2010 */
	  }
	  if( y >= x ) return( i ) ;
	}
	return( lastbranch );   /*  changed 4 Feb 2010 */
}

int
pickbmf_Oton(nsam, mfreq, ptree, tt )
	int nsam, mfreq;
	struct node *ptree;
	double tt;
{
    double x, y, ran1();
	int i, lastbranch = 0 ;
	int tip;
    x = ran1()*tt;

    if(pars.mfs.unfolded != 1) {
         for( i=nsam, y=0; i < 2*nsam-2; i++) {
            if(((ptree+i)->ncs >=3) && ((ptree+i)->ncs != (nsam-1)) && ((ptree+i)->ncs != (nsam-2))) {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ){
                    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                    lastbranch = i ;
                }
                if( y >= x ) return( i ) ;
            }
        }
    } else {
        for( i=nsam, y=0; i < 2*nsam-2; i++) {
            if((ptree+i)->ncs >=3 ) {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ){
                    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                    lastbranch = i ;
                }
                if( y >= x ) return( i ) ;
            }
        }
    }
	return( lastbranch );
}

	int
pickbmf_Dton(nsam, mfreq, ptree, tt )
	int nsam, mfreq;
	struct node *ptree;
	double tt;
{
    double x, y, ran1();
	int i, lastbranch = 0 ;

	int tip;
    x = ran1()*tt;

    if(pars.mfs.unfolded != 1) {
        for( i=nsam, y=0; i < 2*nsam-2 ; i++) {
            if(((ptree+i)->ncs == 2) || ((ptree+i)->ncs == (nsam-2))) {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ){
                    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                    lastbranch = i ;
                }
                if( y >= x ) return( i ) ;
            }
        }
    } else {
    	for( i=nsam, y=0; i < 2*nsam-2 ; i++) {
            if((ptree+i)->ncs == 2) {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ){
                    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                    lastbranch = i ;
                }
                if( y >= x ) return( i ) ;
            }
        }
    }
	//printf("Warning: Have selected the LastNode: %d\n", lastNode);
	return( lastbranch );
}


int
pickbmf_Ston(nsam, mfreq, ptree, tt )
	int nsam, mfreq;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i, lastbranch = 0 ;
    x = ran1()*tt;
	if(pars.mfs.unfolded != 1) {
        for( i=0, y=0; i < 2*nsam-2 ; i++) {
            if(((ptree+i)->ncs == 0) || ((ptree+i)->ncs == (nsam-1))) {
                if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ) {
                    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time;
                    lastbranch = i ;
                }
                if( y >= x ) return( i ) ;
            }
        }
    } else {
        for( i=0, y=0; i < nsam; i++) {
            if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ) {
                y += (ptree + (ptree+i)->abv )->time;
                lastbranch = i ;
            }
            if( y >= x ) return( i ) ;
        }
    }
	return(lastbranch);
}


/****  tdesn : returns 1 if tip is a descendant of node in *ptree, otherwise 0. **/

	int
tdesn(ptree, tip, node )
	struct node *ptree;
	int tip, node;
{
	int k;

	for( k= tip ; k < node ; k = (ptree+k)->abv ) ;
	if( k==node ) return(1);
	else return(0);
}


/* pick2()  */

	int
pick2(n,i,j)
	int n, *i, *j;
{
	double ran1();

	*i = n * ran1() ;
	while( ( *j = n * ran1() ) == *i )
		;
	return(0) ;
}

/**** ordran.c  ***/

	int
ordran(n,pbuf)
	int n;
	double pbuf[];
{
	ranvec(n,pbuf);
	order(n,pbuf);
	return;
}


	int
mnmial(n,nclass,p,rv)
	int n, nclass, rv[];
	double p[];
{
	double ran1();
	double x, s;
	int i, j;

	for(i=0; i<nclass; i++) rv[i]=0;
	for(i=0; i<n ; i++) {
	   x = ran1();
	   j=0;
	   s = p[0];
	   while( (x > s) && ( j<(nclass-1) ) )  s += p[++j];
	   rv[j]++;
	   }
	return(j);
}

        int
order(n,pbuf)
        int n;
        double pbuf[];
{
        int gap, i, j;
        double temp;

        for( gap= n/2; gap>0; gap /= 2)
           for( i=gap; i<n; i++)
                for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
                   temp = pbuf[j];
                   pbuf[j] = pbuf[j+gap];
                   pbuf[j+gap] = temp;
                   }
}


	int
ranvec(n,pbuf)
	int n;
	double pbuf[];
{
	int i;
	double ran1();

	for(i=0; i<n; i++)
		pbuf[i] = ran1();

	return;
}



	int
poisso(u)
	double u;
{
	double  cump, ru, ran1(), p, gasdev(double, double) ;
	int i=1;

	if( u > 30. ){
	    i =  (int)(0.5 + gasdev(u,u)) ;
	    if( i < 0 ) return( 0 ) ;
	    else return( i ) ;
	  }

	ru = ran1();
	p = exp(-u);
	if( ru < p) return(0);
	cump = p;

	while( ru > ( cump += (p *= u/i ) ) )
		i++;
	return(i);
}


/* a slight modification of crecipes version */

double gasdev(m,v)
	double m, v;
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	double ran1();

	if  (iset == 0) {
		do {
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset= v1*fac;
		iset=1;
		return( m + sqrt(v)*v2*fac);
	} else {
		iset=0;
		return( m + sqrt(v)*gset ) ;
	}
}

