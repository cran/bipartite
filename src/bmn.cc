#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
//#include <string.h>

//#include "Tstring.h"
//#include "IOfiles.h"

/*	Input andput files */

/* caracter LF:fin de linea, CE:cero y UN:uno */
#define LF 10
#define CE 48
#define UN 49

//	Para Numerical Recipes
#define NR_END 1
#define FREE_ARG char*

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 3.e-7
#define RNMX (1.0-EPS)

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M_IND 7
#define NSTACK 50
   
#define ITMAX 100
#define EPSB 3.0e-8
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//	El valor más alto de U econtrado es:
const double UMAX = 0.04145;

extern "C" {

	

//	para el algoritmo genetico:
unsigned long NBRS;
int POPSIZE;	// Must be at least 15
int TOURSIZE;
int NBGENER;
int OVER;
int *toc;
int *tor;


/*
 ***************************************************************** 
 */

void matrixSize(char inputFile, int& nrows, int& ncols, int& skip);

void endnote(FILE *f);

void readMatrix(char inputFile, int nrows, int ncols, int skip, int **m);

double matrixTemperature(bool &success,int sp, int **mat, int *c_ord, int *r_ord, int nr, int nc, long& idum);

void orderMatrix(int **mat, int *c_ord, int *r_ord, int nrows, int ncols, int& effInsect, int& effPlant);

void removeBlacks(int **mat, int *c_ord, int *r_ord, int effInsect, int effPlant, int& nbRows, 
				  int& nbCols, double& phi);

void calcZ(double phi, double& z);

void calcDistance(double z, double border[],double **mat,int nr,int nc);

double packMatrix(int sp, int **dataMat, int **intMat, double **d, int *c_ord, int *r_ord, int nr, int nc, 
				  int nrows, int ncols, long& idum);

double calcTemp(double **d, int **mat, int indr[], int indc[], int nr, int nc);

void calcIdiosyncTemp(double **d, int **mat, int indr[], int indc[], int nr, int nc);

void prePackMatrix(int **mat, int indr[], int indc[], int nr, int nc, double x);

void prePackrows(int **mat, int indr[], int indc[], int nr, int nc, double x);

void prePackcols(int **mat, int indr[], int indc[], int nr, int nc, double x);

void prePackNTC(int **mat, int indr[], int indc[], int nr, int nc);

void prePackNTCrows(int **mat, int indr[], int indc[], int nr, int nc);

void prePackNTCcols(int **mat, int indr[], int indc[], int nr, int nc);

void permute(long& idum, int n, int index[]);

void mutate(long& idum, int n, int index[]);

void choosePlayers(long& idum, int n,int m,int arr[]);

void crossOver(long& idum, int nr, int nc, int pr[], int pc[], int or1[], int oc[] );

void indexx(int n, int arr[], int indx[]);

void indexxD(int n, double arr[], int indx[]);

double zbrent(int nr, int nc, double u1, double v1, double z);

double func(double yy, int nr, int nc, double x1, double y1, double z);

double ran1(long& idum);

void avevar(double data[], unsigned long n, double& ave, double& var);

int *ivector(long nl, long nh);

void free_ivector(int *v, long nl, long nh);

double *vector(long nl, long nh);

void free_vector(double *v, long nl, long nh);

int **imatrix(long nrl, long nrh, long ncl, long nch);

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);

double **matrix(long nrl, long nrh, long ncl, long nch);

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);

void nrerror();

/*
 ***************************************************************** 
 */

int bmn(int *matrix, int *n_rows, int *n_cols,double *temperature,
        int *n_nullmodels, int *p_size, int *n_individuals,
        int *n_generations, int *nullmodels,
        double *p_null1, double *avt_null1, double 	*avv_null1,
        double *p_null2, double *avt_null2, double 	*avv_null2,
        double *p_null3, double *avt_null3, double 	*avv_null3,
        int *poc, int *por
        )
{


	

	bool success;
	int ncols,nrows;//,skip;
        
	int i,j;
	int nPresI,nPresJ=0;
	int *colOrder,*rowOrder;
	int **interactMat;
	int **tempMat;
	unsigned long nr;
	long idum=-48367;
	double fill;
	double p,p1,p2,p3,r;
	double tInt,tTemp;
	double ave,var;
	double *pCol,*pRow;
	double *nullTemp;
	
  toc=poc;     // move pointer of column and row order
  tor=por;     // to R objects in bmn

//int *mat = matrix;



/*
	int stoch;
	fprintf(stdout, "This program can run in two modes: stochastic and deterministic.\n");
	fprintf(stdout, "In the deterministic mode, if you analyse the same matrix several times\n");
	fprintf(stdout, "it will always give the same results. In stochastic mode, \n");
	fprintf(stdout, "it will find slightly different results each time you analyse the data.\n");
	fprintf(stdout, "Type 1 if you want to use the stochastic mode, \n");
	fprintf(stdout, "and type any other number if you want to use the deterministic version.\n");
	cin >> stoch;
	fprintf(stdout, "\n\n");

	if (stoch==1)
	{
		idum = time(0);			// auto seeding using system time
		idum *= -1;
	}
*/

//char  INFILE='i';

//char OUTFILE= 'o';


nrows=*n_rows;
ncols=*n_cols;

NBRS =*n_nullmodels;

POPSIZE =*p_size;



TOURSIZE =*n_individuals;

	
NBGENER = *n_generations;

//	fprintf(stdout, "\n\n\n");
//	fprintf(stdout, "BINMATNEST is calculating the temperature of your matrix\n");
//	fprintf(stdout, "and the corresponding p-values. Please be patient.\n");
//	fprintf(stdout, "\n\n\n");

//	matrixSize(INFILE,nrows,ncols,skip);
	interactMat = imatrix(1,nrows,1,ncols);
//	readMatrix(INFILE,nrows,ncols,skip,interactMat);

// hand over R matrix....


for (int j=1; j<=ncols;j++)
	{
for (int i=1; i<=nrows;i++)
{
	interactMat[i][j] = matrix[(j-1)*nrows+(i-1)];
	}
}



	//	Some variables are initialised
	fill = 0.0;
	for (i=1;i<=nrows;i++)
	{
		for (j=1;j<=ncols;j++) 
		{
			fill += interactMat[i][j];
		}
	}
	fill /= (1.0*nrows*ncols);

	pRow = vector(1,nrows);
	pCol = vector(1,ncols);
	for (i=1;i<=nrows;i++)
	{
		pRow[i] = 0;
		for (j=1;j<=ncols;j++) pRow[i] += interactMat[i][j];
		pRow[i] /= 1.0*ncols;
	}
	for (j=1;j<=ncols;j++)
	{
		pCol[j] = 0;
		for (i=1;i<=nrows;i++) pCol[j] += interactMat[i][j];
		pCol[j] /= 1.0*nrows;
	}

	colOrder = ivector(1,ncols);
	rowOrder = ivector(1,nrows);

	

	tInt = matrixTemperature(success,1,interactMat,colOrder,rowOrder,nrows,ncols,idum);

	//hand over matrix temperature to R
	
	*temperature = tInt;

  // return packing order by pointers (tor and toc via poc and por)
  

  

	free_imatrix(interactMat,1,nrows,1,ncols);


  if (*nullmodels==1) {  //should nullmodels be calculated??
  
	tempMat = imatrix(1,nrows,1,ncols);
	nullTemp = vector(1,NBRS);
	//	First null model:
	p1 = 0.0;
	for (nr=1;nr<=NBRS;nr++)
	{
		do {
			for (i=1;i<=nrows;i++)
			{
				nPresI = 0;
				while (nPresI==0)
				{
					for (j=1;j<=ncols;j++) 
					{
						p = fill;
						r = ran1(idum);
						tempMat[i][j] = (r<p ? 1 : 0);
						nPresI += tempMat[i][j];
					}
				}
			}
			for (j=1;j<=ncols;j++)
			{
				nPresJ = 0;
				for (i=1;i<=nrows;i++) nPresJ += tempMat[i][j];
				if (nPresJ==0) break;
			}
		} while (nPresJ==0);
		tTemp = matrixTemperature(success,0,tempMat,colOrder,rowOrder,nrows,ncols,idum);
		if (success)
		{
			nullTemp[nr] = tTemp;
			if (tTemp < tInt) p1 += 1.0;
		}
		else nr--;
	}
	p1 /= (1.0*NBRS);
	avevar(nullTemp,NBRS,ave,var);
	
	*p_null1 = p1;
	*avt_null1 = ave;
	*avv_null1 = var;
	
	//out = fopen(OUTFILE, "a");
	//out = fopen("OUTFILE", "a");
	//fprintf(out,"  Null model     p-value   Average T    Variance\n");
	//fprintf(out,"       First %11.5f %11.5f %11.5f\n",p1,ave,var);
	//fprintf(stdout,"First null model:  p1 = %11.5f\n",p1);

	//close_file(out);
	//fclose(out);

	//	Second null model:
	p2 = 0.0;
	for (nr=1;nr<=NBRS;nr++)
	{
		do {
			for (i=1;i<=nrows;i++)
			{
				nPresI = 0;
				while (nPresI==0)
				{
					for (j=1;j<=ncols;j++) 
					{
						p = pCol[j];
						r = ran1(idum);
						tempMat[i][j] = (r<p ? 1 : 0);
						nPresI += tempMat[i][j];
					}
				}
			}
			for (j=1;j<=ncols;j++)
			{
				nPresJ = 0;
				for (i=1;i<=nrows;i++) nPresJ += tempMat[i][j];
				if (nPresJ==0) break;
			}
		} while (nPresJ==0);
		tTemp = matrixTemperature(success,0,tempMat,colOrder,rowOrder,nrows,ncols,idum);
		if (success)
		{
			nullTemp[nr] = tTemp;
			if (tTemp < tInt) p2 += 1.0;
		}
		else nr--;
	}
	p2 /= (1.0*NBRS);
	avevar(nullTemp,NBRS,ave,var);
	
		*p_null2 = p2;
	*avt_null2 = ave;
	*avv_null2 = var;
	//out = open_append(OUTFILE);
	//out = fopen("OUTFILE", "a");
	//fprintf(out,"      Second %11.5f %11.5f %11.5f\n",p2,ave,var);
	//fprintf(stdout,"Second null model: p2 = %11.5f\n",p2);
	//close_file(out);
	//fclose(out);
	//	Third null model:
	p3 = 0.0;
	for (nr=1;nr<=NBRS;nr++)
	{
		do {
			for (i=1;i<=nrows;i++)
			{
				nPresI = 0;
				while (nPresI==0)
				{
					for (j=1;j<=ncols;j++) 
					{
						p = 0.5*(pRow[i]+pCol[j]);
						r = ran1(idum);
						tempMat[i][j] = (r<p ? 1 : 0);
						nPresI += tempMat[i][j];
					}
				}
			}
			for (j=1;j<=ncols;j++)
			{
				nPresJ = 0;
				for (i=1;i<=nrows;i++) nPresJ += tempMat[i][j];
				if (nPresJ==0) break;
			}
		} while (nPresJ==0);
		tTemp = matrixTemperature(success,0,tempMat,colOrder,rowOrder,nrows,ncols,idum);
		if (success)
		{
			nullTemp[nr] = tTemp;
			if (tTemp < tInt) p3 += 1.0;
		}
		else nr--;
	}
	p3 /= (1.0*NBRS);
	avevar(nullTemp,NBRS,ave,var);
	
		*p_null3 = p3;
	*avt_null3 = ave;
	*avv_null3 = var;
	//out = open_append(OUTFILE);
	//out = fopen("OUTFILE", "a");
	//fprintf(out,"       Third %11.5f %11.5f %11.5f\n",p3,ave,var);
	//fprintf(stdout,"Third null model:  p3 = %11.5f\n",p3);
	//close_file(out);
	//fclose(out);

	free_ivector(rowOrder,1,nrows);
	free_ivector(colOrder,1,ncols);
	free_imatrix(tempMat,1,nrows,1,ncols);
	free_vector(pRow,1,nrows);
	free_vector(pCol,1,ncols);
	free_vector(nullTemp,1,NBRS);
 }
 // end of nullmodels.....
	//fprintf(stdout, "\n\n\n");
	//fprintf(stdout, "BINMATNEST has saved the results.\n");
	//fprintf(stdout, "Type a number and press ENTER to continue\n");
//        cin >> TOURSIZE;

        //res[1] = 1.11;
        //matrix[8] = 2.234;
	return EXIT_SUCCESS;

}
/***
 *****************************************
 ***/
void matrixSize(char inputFile, int& nrows, int& ncols, int& skip)
{


skip=1;
}
/***
 *****************************************
 ***/
void readMatrix(char inputFile, int nrows, int ncols, int skip, int **m)
{



	
}
/***
 *****************************************
 ***/
void endnote(FILE *f)
{

}
/***
 *****************************************
 ***/
double matrixTemperature(bool &success,int sp, int **dataMat, int *c_ord, int *r_ord, int nr, int nc, long& idum)
{
	int i,j;
	int effPlant,effInsect;
	int nbRows,nbCols;
	static int count = 0;
	int **mat;
	int **packed;
	double phi,z;
	double tMat=0.0;
	double *border;
	double **distance;

	success = true;
	mat = imatrix(1,nr,1,nc);
	for (i=1;i<=nr;i++)
	{
		for (j=1;j<=nc;j++) 
		{
			mat[i][j] = dataMat[i][j];
		}
	}

	//	Permuting rows and columns for pre-packing
	orderMatrix(mat,c_ord,r_ord,nr,nc,effInsect,effPlant);
	removeBlacks(mat,c_ord,r_ord,effInsect,effPlant,nbRows,nbCols,phi);
	packed = imatrix(1,nbRows,1,nbCols);
	for (i=1;i<=nbRows;i++)
	{
		for (j=1;j<=nbCols;j++) 
		{
			packed[i][j] = mat[i][j];
		}
	}

	if (nbCols>2 && nbRows>2)
	{
		border = vector(1,nbCols);
		calcZ(phi,z);

		//	Calculate distance matrix
		distance = matrix(1,nbRows,1,nbCols);
		calcDistance(z,border,distance,nbRows,nbCols);
		//	Final packing using GA. Returns the temperature of the matrix
		tMat = packMatrix(sp,dataMat,packed,distance,c_ord,r_ord,nbRows,nbCols,nr,nc,idum);

		free_matrix(distance,1,nbRows,1,nbCols);
		free_vector(border,1,nbCols);
	}
	else 
	{
		if (sp) nrerror();
		else
		{
			tMat = 0;
			success = false;
			count++;
			if (count>1000) nrerror();
		}
	}
	
	free_imatrix(mat,1,nr,1,nc);
	free_imatrix(packed,1,nbRows,1,nbCols);

	return tMat;
}
/***
 *****************************************
 ***/
void orderMatrix(int **mat, int *c_ord, int *r_ord, int nrows, int ncols, int& effInsect, int& effPlant)
{
	//	We permute rows and columns in such a way that rows with more interactions are closer to the top, and
	//	columns with more interactions are closer to the left.
	int i,j;
	int *lP,*lI;
	int *indP,*indI;
	int **temp;

	lP = ivector(1,ncols);
	indP = ivector(1,ncols);
	lI = ivector(1,nrows);
	indI = ivector(1,nrows);
	temp = imatrix(1,nrows,1,ncols);

	effInsect = effPlant = 0;

	//	We count the number of insects that interact with at least one plant
	for (i=1;i<=nrows;i++)
	{
		indI[i] = i;
		lI[i] = 0;
		for (j=1;j<=ncols;j++) lI[i] -= mat[i][j];
		if (lI[i]<0) effInsect++;
	}
	indexx(nrows,lI,indI);
	for (i=1;i<=nrows;i++) r_ord[i] = indI[i];

	//	We count the number of plants that interact with at least one insect
	for (j=1;j<=ncols;j++)
	{
		indP[j] = j;
		lP[j] = 0;
		for (i=1;i<=nrows;i++) lP[j] -= mat[i][j];
		if (lP[j]<0) effPlant++;
	}
	indexx(ncols,lP,indP);
	for (j=1;j<=ncols;j++) c_ord[j] = indP[j];

	for (i=1;i<=nrows;i++)
	{
		for (j=1;j<=ncols;j++) temp[i][j] = mat[i][j];
	}
	for (i=1;i<=nrows;i++)
	{
		for (j=1;j<=ncols;j++) mat[i][j]=temp[indI[i]][indP[j]];
	}

	free_ivector(lP,1,ncols);
	free_ivector(indP,1,ncols);
	free_ivector(lI,1,nrows);
	free_ivector(indI,1,nrows);
	free_imatrix(temp,1,nrows,1,ncols);
}
/***
 *****************************************
 ***/
void removeBlacks(int **mat, int *c_ord, int *r_ord, int effInsect, int effPlant, int& nbRows, 
				  int& nbCols, double& phi)
{
	int i,j,k;
	int blackInsect=0,blackPlant=0;
	int nbinter=0;
	int **temp;

	temp = imatrix(1,effInsect,1,effPlant);

	//	We count the number of black rows
	for (i=1;i<=effInsect;i++)
	{
		k = 0;
		for (j=1;j<=effPlant;j++) k += mat[i][j];
		if (k==effPlant) blackInsect++;
		else break;
	}
	nbRows = (blackInsect>1 ? effInsect - blackInsect + 1 : effInsect);

	//	We count the number of black columns
	for (j=1;j<=effPlant;j++)
	{
		k = 0;
		for (i=1;i<=effInsect;i++) k += mat[i][j];
		if (k==effInsect) blackPlant++;
		else break;
	}
	nbCols = (blackPlant>1 ? effPlant - blackPlant + 1 : effPlant);

	if (blackInsect>1)
	{
		for (i=1;i<=effInsect;i++)
		{
			for (j=1;j<=effPlant;j++) temp[i][j] = mat[i][j];
		}
		for (i=1;i<=effInsect-blackInsect+1;i++)
		{
			r_ord[i] = r_ord[i+blackInsect-1];
			for (j=1;j<=effPlant;j++) mat[i][j] = temp[i+blackInsect-1][j];
		}
		for (i=effInsect-blackInsect+2;i<=effInsect;i++)
		{
			r_ord[i] = 0;
			for (j=1;j<=effPlant;j++) mat[i][j] = 0;
		}
	}

	if (blackPlant>1)
	{
		for (j=1;j<=effPlant;j++)
		{
			for (i=1;i<=effInsect;i++) temp[i][j] = mat[i][j];
		}
		for (j=1;j<=effPlant-blackPlant+1;j++)
		{
			c_ord[j] = c_ord[j+blackPlant-1];
			for (i=1;i<=effInsect;i++) mat[i][j] = temp[i][j+blackPlant-1];
		}
		for (j=effPlant-blackPlant+2;j<=effPlant;j++)
		{
			c_ord[j] = 0;
			for (i=1;i<=effInsect;i++) mat[i][j] = 0;
		}
	}

	//	We count the number of interactions
	for (i=1;i<=nbRows;i++)
	{
		for (j=1;j<=nbCols;j++) nbinter += mat[i][j];
	}
	if (nbRows>1 && nbCols>1) phi = (2.0*nbinter-nbRows-nbCols+1)/(2.0*(nbRows-1)*(nbCols-1));
	else phi = -1;

	free_imatrix(temp,1,effInsect,1,effPlant);
}
/***
 *****************************************
 ***/
void calcZ(double phi, double& z)
	//	Calculate the parameter for the "line of perfect order"
{
	int i;
	static double znVal[41] = {0.2, 0.224, 0.2509, 0.281, 0.3147, 0.3525, 0.3948, 0.4421, 0.4952, 
		0.5546, 0.6212, 0.6957, 0.7792, 0.8727, 0.9774, 1.0947, 1.2261, 1.3732, 1.538, 1.7226, 
		1.9293, 2.1608, 2.4201, 2.7105, 3.0357, 3.4, 3.808, 4.265, 4.7768, 5.35, 5.992, 6.711, 
		7.5163, 8.4183, 9.4285, 10.5599, 11.8271, 13.2464, 14.8359, 16.6162, 18.6102};
	static double propOc[41] = {0.996, 0.9921, 0.9855, 0.9751, 0.9599, 0.9389, 0.9116, 0.8776, 0.8371, 
		0.7907, 0.7394, 0.6844, 0.6272, 0.5691, 0.5114, 0.4555, 0.402, 0.3519, 0.3057, 0.2636, 0.2258, 
		0.1922, 0.1626, 0.1368, 0.1146, 0.0956, 0.0793, 0.0656, 0.0541, 0.0445, 0.0365, 0.0297, 
		0.0242, 0.0197, 0.0161, 0.013, 0.0106, 0.0086, 0.007, 0.0057, 0.0046};

	if (phi>=1) z = 1000;
	else if (phi>0)
	{
		if (phi>=propOc[0])
		{
			z = znVal[0]*(1-phi)/(1-propOc[0]);
		}
		else if (phi<=propOc[40])
		{
			z = znVal[40];
		}
		else
		{
			for (i=1;i<=40;i++)
			{
				if (propOc[i]<=phi) break;
			}
			z = znVal[i-1]+(znVal[i]-znVal[i-1])*(propOc[i-1]-phi)/(propOc[i-1]-propOc[i]);
		}
	}
	else z = -1;
}
/***
 *****************************************
 ***/
void calcDistance(double z,double border[],double **mat,int nr,int nc)
{
	int i,j;
	double d2,dmax2;
	double s,x,x0,x1,xx,y,y0,y1,yy;

	if (z>0 && z<100)
	{
		for (j=1;j<=nc;j++)
		{
			x = (j-0.5)/nc;
			xx = (x-0.5/nc)*nc/(nc-1);
			s = pow(1-xx,z);
			yy = pow(1-s,1/z);
			y = (0.5+(nr-1)*yy)/nr;
			border[j] = y*nr;
		}

		for (i=1;i<=nr;i++)
		{
			for (j=1;j<=nc;j++)
			{
				if ((i==1 && j==nc) || (i==nr && j==1)) mat[i][j] = 0;
				else
				{
					x1 = (j-0.5)/nc;
					y1 = (nr-i+0.5)/nr;
//					dmax2 = (x1+y1<1 ? 2*(x1+y1)*(x1+y1) : 2*(2-x1-y1)*(2-x1-y1));
//					but we calculate only half of the distance:
					dmax2 = (x1+y1<1 ? (x1+y1)*(x1+y1) : (2-x1-y1)*(2-x1-y1));
					yy = zbrent(nr,nc,x1,y1,z);
					y = (0.5+(nr-1)*yy)/nr;
/*	
					d2 = (x1-x)*(x1-x) + (y1-y)*(y1-y);
					with:
					x = x1 + y1 - y;
					so that x1-x = y-y1
					and d2 = (y-y1)*(y-y1) + (y1-y)*(y1-y) = 2*(y1-y)*(y1-y);
					We calculate only half of the distance:
*/	
					d2 = (y1-y)*(y1-y);
					mat[i][j] = d2/dmax2;
					if ((nr-i+0.5)<border[j]) mat[i][j] *= -1;
				}
			}
		}
	}
	else if (z>100)
	{
		for (j=1;j<=nc;j++) border[j] = 0.5;

		for (i=1;i<=nr;i++)
		{
			for (j=1;j<=nc;j++)
			{
				if ((i==1 && j==nc) || (i==nr && j==1)) mat[i][j] = 0;
				else
				{
					x1 = j-0.5;
					y1 = nr-i+0.5;
					if (y1<=nr*(1-x1/nc))
					{
						x0 = x1 + (nc*y1)/nr;
						y0 = y1 + (nr*x1)/nc;
						dmax2 = x0*x0 + y0*y0;
						d2 = (1+(1.0*nc*nc)/(nr*nr))*(nr-i)*(nr-i);
						mat[i][j] = d2/dmax2;
						if ((nr-i+0.5)<border[j]) mat[i][j] *= -1;
					}
					else
					{
						x0 = 2*nr - x1 - (nc*y1)/nr;
						y0 = 2*nc - y1 - (nr*x1)/nc;
						dmax2 = x0*x0 + y0*y0;
						d2 = (1+(1.0*nr*nr)/(nc*nc))*(nc-j)*(nc-j);
						mat[i][j] = d2/dmax2;
						if ((nr-i+0.5)<border[j]) mat[i][j] *= -1;
					}
				}
			}
		}
	}
	else
	{
		border[1] = 0.5;
		for (j=2;j<=nc;j++) border[j] = nr - 0.5;

		for (i=1;i<=nr;i++)
		{
			for (j=1;j<=nc;j++)
			{
				if ((i==1 && j==nc) || (i==nr && j==1)) mat[i][j] = 0;
				else
				{
					x1 = j-0.5;
					y1 = nr-i+0.5;
					if (y1<=nr*(1-x1/nc))
					{
						x0 = x1 + (nc*y1)/nr;
						y0 = y1 + (nr*x1)/nc;
						dmax2 = x0*x0 + y0*y0;
						d2 = (1+(1.0*nr*nr)/(nc*nc))*(j-1)*(j-1);
						mat[i][j] = d2/dmax2;
						if ((nr-i+0.5)<border[j]) mat[i][j] *= -1;
					}
					else
					{
						x0 = 2*nr - x1 - (nc*y1)/nr;
						y0 = 2*nc - y1 - (nr*x1)/nc;
						dmax2 = x0*x0 + y0*y0;
						d2 = (1+(1.0*nc*nc)/(nr*nr))*(i-1)*(i-1);
						mat[i][j] = d2/dmax2;
						if ((nr-i+0.5)<border[j]) mat[i][j] *= -1;
					}
				}
			}
		}
	}
}
/***
 *****************************************
 ***/
double packMatrix(int sp, int **dataMat, int **intMat, double **d, int *c_ord, int *r_ord, int nr, int nc, 
				  int nrows, int ncols, long& idum)
	//	indr[i] indica la fila de la matriz que deben ocupar los datos del insecto i
	//	indc[j] indica la columna de la matriz que deben ocupar los datos de la planta j
{
	int i,j,l,k1,k2,ng;
	int *indr,*indc;
	int *pr,*or1,*pc,*oc;
	int *players,*sorted;
	int *bestr,*bestc;
	int *unused;
	int **pRows,**pCols;
	int **mat;
	double tMin;
	double *temp,*tPlayers;


	indr = ivector(1,nr);
	indc = ivector(1,nc);
	pr = ivector(1,nr);
	pc = ivector(1,nc);
	or1 = ivector(1,nr);
	oc = ivector(1,nc);
	bestr = ivector(1,nr);
	bestc = ivector(1,nc);
	players = ivector(1,TOURSIZE);
	sorted = ivector(1,TOURSIZE);
	tPlayers = vector(1,TOURSIZE);
	temp = vector(1,POPSIZE);
	pRows = imatrix(1,POPSIZE,1,nr);
	pCols = imatrix(1,POPSIZE,1,nc);

	for (i=1;i<=nr;i++) pRows[1][i] = indr[i] = bestr[i] = i;
	for (j=1;j<=nc;j++) pCols[1][j] = indc[j] = bestc[j] = j;
	tMin = temp[1] = calcTemp(d,intMat,indr,indc,nr,nc);
	for (l=2;l<=12;l++)
	{
		for (i=1;i<=nr;i++) pRows[l][i] = i;
		for (j=1;j<=nc;j++) pCols[l][j] = j;
		prePackMatrix(intMat,pRows[l],pCols[l],nr,nc,0.1*(l-2));
	}
	//	this one is created with the algorithm from NTC
	for (i=1;i<=nr;i++) pRows[13][i] = i;
	for (j=1;j<=nc;j++) pCols[13][j] = j;
	prePackNTC(intMat,pRows[13],pCols[13],nr,nc);
	for (l=14;l<=POPSIZE;l++)
	{
		k1 = 1 + l%13;
		for (i=1;i<=nr;i++) pRows[l][i] = pRows[k1][i];
		for (j=1;j<=nc;j++) pCols[l][j] = pCols[k1][j];
		if (ran1(idum)<0.7) mutate(idum,nr,pRows[l]);
		if (ran1(idum)<0.7) mutate(idum,nc,pCols[l]);
	}

	for (l=2;l<=POPSIZE;l++)
	{
		temp[l] = calcTemp(d,intMat,pRows[l],pCols[l],nr,nc);
		if (temp[l]<tMin)
		{
			tMin = temp[l];
			for (i=1;i<=nr;i++) bestr[i] = pRows[l][i];
			for (j=1;j<=nc;j++) bestc[j] = pCols[l][j];
		}
	}

	for (ng=1;ng<=NBGENER;ng++)
	{
		choosePlayers(idum,TOURSIZE,POPSIZE,players);
		for (l=1;l<=TOURSIZE;l++) tPlayers[l] = temp[players[l]];
		indexxD(TOURSIZE,tPlayers,sorted);

		k1 = players[sorted[1]];
		for (i=1;i<=nr;i++) or1[i] = pRows[k1][i];
		for (j=1;j<=nc;j++) oc[j] = pCols[k1][j];
		do
		{
			k2 = 1 + int(POPSIZE*ran1(idum));
		} while (k1==k2 || k2>POPSIZE);
		for (i=1;i<=nr;i++) pr[i] = pRows[k2][i];
		for (j=1;j<=nc;j++) pc[j] = pCols[k2][j];
		crossOver(idum,nr,nc,pr,pc,or1,oc);
		l = players[sorted[TOURSIZE]];
		for (i=1;i<=nr;i++) pRows[l][i] = or1[i];
		for (j=1;j<=nc;j++) pCols[l][j] = oc[j];
		temp[l] = calcTemp(d,intMat,or1,oc,nr,nc);
		if (temp[l]<tMin)
		{
			tMin = temp[l];
			for (i=1;i<=nr;i++) bestr[i] = or1[i];
			for (j=1;j<=nc;j++) bestc[j] = oc[j];
		}

		k1 = players[sorted[2]];
		for (i=1;i<=nr;i++) or1[i] = pRows[k1][i];
		for (j=1;j<=nc;j++) oc[j] = pCols[k1][j];
		do
		{
			k2 = 1 + int(POPSIZE*ran1(idum));
		} while (k1==k2 || k2>POPSIZE);
		for (i=1;i<=nr;i++) pr[i] = pRows[k2][i];
		for (j=1;j<=nc;j++) pc[j] = pCols[k2][j];
		crossOver(idum,nr,nc,pr,pc,or1,oc);
		l = players[sorted[TOURSIZE-1]];
		for (i=1;i<=nr;i++) pRows[l][i] = or1[i];
		for (j=1;j<=nc;j++) pCols[l][j] = oc[j];
		temp[l] = calcTemp(d,intMat,or1,oc,nr,nc);
		if (temp[l]<tMin)
		{
			tMin = temp[l];
			for (i=1;i<=nr;i++) bestr[i] = or1[i];
			for (j=1;j<=nc;j++) bestc[j] = oc[j];
		}
	}

	if (sp)
	{
		//out = open_append(OUTFILE);
		//out = fopen("OUTFILE", "a");
		//fprintf(out,"Packed matrix:\n");
		for (i=1;i<=nr;i++)
		{
			for (j=1;j<=nc;j++) 
			{
				;//fprintf(out,"%1i",intMat[bestr[i]][bestc[j]]);
			}
			//fprintf(out,"\n");
		}
		//fprintf(out,"\n");
        
		//close_file(out);
		//fclose(out);
		calcIdiosyncTemp(d,intMat,or1,oc,nr,nc);
	}

	free_ivector(indr,1,nr);
	free_ivector(indc,1,nc);
	free_ivector(pr,1,nr);
	free_ivector(pc,1,nc);
	free_ivector(or1,1,nr);
	free_ivector(oc,1,nc);
	free_ivector(players,1,TOURSIZE);
	free_ivector(sorted,1,TOURSIZE);
	free_vector(tPlayers,1,TOURSIZE);
	free_vector(temp,1,POPSIZE);
	free_imatrix(pRows,1,POPSIZE,1,nr);
	free_imatrix(pCols,1,POPSIZE,1,nc);

	if (sp)
	{
		//out = open_append(OUTFILE);
		//out = fopen("OUTFILE", "a");

		sorted = ivector(1,nr);
		for (i=1;i<=nr;i++) sorted[i] = r_ord[i];
		for (i=1;i<=nr;i++) r_ord[i] = sorted[bestr[i]];
		free_ivector(sorted,1,nr);

		unused = ivector(1,nrows);
		for (i=1;i<=nrows;i++) unused[i] = 1;
		for (i=1;i<=nr;i++) unused[r_ord[i]] = 0;
		for (i=1;i<=nrows;i++) 
		{
			if (unused[i])
			{
				k1 = 0;
				for (j=1;j<=ncols;j++) k1 += dataMat[i][j];
				unused[i] = (k1 > 0? 1 : nr+1);
			}
		}

		//fprintf(out,"Permutation of rows:\n");
		//fprintf(out,"Row number            ");
		for (i=1;i<=nrows;i++)
		{
			if (unused[i]==1) ;//fprintf(out,"%5i",unused[i]);
		}
		for (i=1;i<=nr;i++) ;//fprintf(out,"%5i",i);
		for (i=1;i<=nrows;i++)
		{
			if (unused[i]>1);// fprintf(out,"%5i",unused[i]);
		}
		//fprintf(out,"\n");
		//fprintf(out,"comes from position   ");
		for (i=1;i<=nrows;i++)
		{
			if (unused[i]==1) tor[i-1]=i;//fprintf(out,"%5i",i);
		}
		for (i=1;i<=nr;i++) tor[i-1] = r_ord[i];// fprintf(out,"%5i",r_ord[i]);
		for (i=1;i<=nrows;i++)
		{
			if (unused[i]>1) tor[i-1]=i;//fprintf(out,"%5i",i);
		}
		//fprintf(out,"\n\n");
		free_ivector(unused,1,nrows);

		sorted = ivector(1,nc);
		for (j=1;j<=nc;j++) sorted[j] = c_ord[j];
		for (j=1;j<=nc;j++) c_ord[j] = sorted[bestc[j]];
		free_ivector(sorted,1,nc);

		unused = ivector(1,ncols);
		for (j=1;j<=ncols;j++) unused[j] = 1;
		for (j=1;j<=nc;j++) unused[c_ord[j]] = 0;
		for (j=1;j<=ncols;j++) 
		{
			if (unused[j])
			{
				k1 = 0;
				for (i=1;i<=nrows;i++) k1 += dataMat[i][j];
				unused[j] = (k1 > 0? 1 : nc+1);
			}
		}

		//fprintf(out,"Permutation of columns:\n");
		//fprintf(out,"Column number         ");
		for (j=1;j<=ncols;j++)
		{
			if (unused[j]==1);// fprintf(out,"%5i",unused[j]);
		}
		for (j=1;j<=nc;j++) ;//fprintf(out,"%5i",j);
		for (j=1;j<=ncols;j++)
		{
			if (unused[j]>1) ;//fprintf(out,"%5i",unused[j]);
		}
		//fprintf(out,"\n");
		//fprintf(out,"comes from position   ");

		for (j=1;j<=ncols;j++)
		{
			if (unused[j]==1)  toc[j-1] = j;// fprintf(out,"%5i",j);
		}
		for (j=1;j<=nc;j++) toc[j-1] = c_ord[j];//c_ord[j];// fprintf(out,"%5i",c_ord[j]);
		for (j=1;j<=ncols;j++)
		{
			if (unused[j]>1)  toc[j-1] = j;// fprintf(out,"%5i",j);
		}
		//fprintf(out,"\n\n");
		free_ivector(unused,1,ncols);

		//close_file(out);
		//fclose(out);
	}

	mat = imatrix(1,nr,1,nc);
	for (i=1;i<=nr;i++)
	{
		for (j=1;j<=nc;j++) mat[i][j] = intMat[i][j];
	}
	for (i=1;i<=nr;i++)
	{
		for (j=1;j<=nc;j++) intMat[i][j] = mat[bestr[i]][bestc[j]];
	}
	free_imatrix(mat,1,nr,1,nc);
	free_ivector(bestr,1,nr);
	free_ivector(bestc,1,nc);
	return tMin;
}
/***
 *****************************************
 ***/
void prePackMatrix(int **mat, int indr[], int indc[], int nr, int nc, double x)
{
	int fold;

	if (nc > nr)
	{
		for (fold=1;fold<=4;fold++)
		{
			prePackcols(mat,indr,indc,nr,nc,x);
			prePackrows(mat,indr,indc,nr,nc,x);
		}
	}
	else
	{
		for (fold=1;fold<=4;fold++)
		{
			prePackrows(mat,indr,indc,nr,nc,x);
			prePackcols(mat,indr,indc,nr,nc,x);
		}
	}
}
/***
 *****************************************
 ***/
void prePackcols(int **mat, int indr[], int indc[], int nr, int nc, double x)
{
	int i,j;
	double *w;

	w = vector(1,nc);

	//	We calculate an index for each column (plant)
	for (j=1;j<=nc;j++)
	{
		w[j] = 0;
		for (i=1;i<=nr;i++) 
		{
			if (mat[indr[i]][j]==1) w[j] -= x*i*i;
			else w[j] += (1-x)*(nr-i+1)*(nr-i+1);
		}
	}

	indexxD(nc,w,indc);

	free_vector(w,1,nc);
}
/***
 *****************************************
 ***/
void prePackrows(int **mat, int indr[], int indc[], int nr, int nc, double x)
{
	int i,j;
	double *w;

	w = vector(1,nr);

	//	We calculate an index for each column (plant)
	for (i=1;i<=nr;i++)
	{
		w[i] = 0;
		for (j=1;j<=nc;j++)
		{
			if (mat[i][indc[j]]==1) w[i] -= x*j*j;
			else w[i] += (1-x)*(nc-j+1)*(nc-j+1);
		}
	}

	indexxD(nr,w,indr);

	free_vector(w,1,nr);
}
/***
 *****************************************
 ***/
void prePackNTC(int **mat, int indr[], int indc[], int nr, int nc)
{
	int fold;

	if (nc > nr)
	{
		for (fold=1;fold<=8;fold++)
		{
			prePackNTCcols(mat,indr,indc,nr,nc);
			prePackNTCrows(mat,indr,indc,nr,nc);
		}
	}
	else
	{
		for (fold=1;fold<=8;fold++)
		{
			prePackNTCrows(mat,indr,indc,nr,nc);
			prePackNTCcols(mat,indr,indc,nr,nc);
		}
	}
}
/***
 *****************************************
 ***/
void prePackNTCcols(int **mat, int indr[], int indc[], int nr, int nc)
	//	indr[i] indica la fila de la matriz que deben ocupar los datos del insecto i
	//	indc[j] indica la columna de la matriz que deben ocupar los datos de la planta j
{
	int i,j;
	int *ndum,*ns,*nt;
	double *s,*t;

	ndum = ivector(1,nc);
	ns = ivector(1,nc);
	nt = ivector(1,nc);
	s = vector(1,nc);
	t = vector(1,nc);

	//	We calculate an index for each column (plant)
	for (j=1;j<=nc;j++)
	{
		s[j] = t[j] = 0;
		for (i=1;i<=nr;i++) 
		{
			if (mat[indr[i]][j]==1) s[j] -= i*i;
			else t[j] += (nr-i+1)*(nr-i+1);
		}
	}
	indexxD(nc,s,ndum);
	indexx(nc,ndum,ns);
	indexxD(nc,t,ndum);
	indexx(nc,ndum,nt);
	for (j=1;j<=nc;j++) s[j] = ns[j] + nt[j];
	indexxD(nc,s,indc);

	free_ivector(ndum,1,nc);
	free_ivector(ns,1,nc);
	free_ivector(nt,1,nc);
	free_vector(s,1,nc);
	free_vector(t,1,nc);
}
/***
 *****************************************
 ***/
void prePackNTCrows(int **mat, int indr[], int indc[], int nr, int nc)
	//	indr[i] indica la fila de la matriz que deben ocupar los datos del insecto i
	//	indc[j] indica la columna de la matriz que deben ocupar los datos de la planta j
{
	int i,j;
	int *ndum,*ns,*nt;
	double *s,*t;

	ndum = ivector(1,nr);
	ns = ivector(1,nr);
	nt = ivector(1,nr);
	s = vector(1,nr);
	t = vector(1,nr);

	//	We calculate an index for each column (plant)
	for (i=1;i<=nr;i++)
	{
		s[i] = t[i] = 0;
		for (j=1;j<=nc;j++) 
		{
			if (mat[i][indc[j]]==1) s[i] -= j*j;
			else t[i] += (nc-j+1)*(nc-j+1);
		}
	}
	indexxD(nr,s,ndum);
	indexx(nr,ndum,ns);
	indexxD(nr,t,ndum);
	indexx(nr,ndum,nt);
	for (i=1;i<=nr;i++) s[i] = ns[i] + nt[i];
	indexxD(nr,s,indr);
}
/***
 *****************************************
 ***/
double calcTemp(double **d, int **mat, int indr[], int indc[], int nr, int nc)
{
	//	indr[i] indica la fila de la matriz que deben ocupar los datos del insecto i
	//	indc[j] indica la columna de la matriz que deben ocupar los datos de la planta j
	int i,j;
	int i1,j1;
	double unex;
	double tMat;
	unex = 0.0;
	for (i=1;i<=nr;i++)
	{
		i1 = indr[i];
		for (j=1;j<=nc;j++)
		{
			j1 = indc[j];
			if ((mat[i1][j1]==1 && d[i][j]<0) || (mat[i1][j1]==0 && d[i][j]>0))
			{
				unex += fabs(d[i][j]);
			}
		}
	}
	unex /= (nr*nc);
	tMat = 100*unex/UMAX;

	return tMat;
}
/***
 *****************************************
 ***/
void calcIdiosyncTemp(double **d, int **mat, int indr[], int indc[], int nr, int nc)
{
	//	indr[i] indica la fila de la matriz que deben ocupar los datos del insecto i
	//	indc[j] indica la columna de la matriz que deben ocupar los datos de la planta j
	int i,j;
	int i1,j1;
	double unex;
	double tIdiosync;


	//out = open_append(OUTFILE);
	//out = fopen("OUTFILE", "a");
	//	Idiosyncratic temperature for rows:
	//fprintf(out,"Idiosyncratic temperature for rows:\n");
	//fprintf(out,"Row:                       ");
	for (i=1;i<=nr;i++) ;//fprintf(out,"%10i",i);
	//fprintf(out,"\n");
	//fprintf(out,"Idiosyncratic temperature: ");
	for (i=1;i<=nr;i++)
	{
		unex = 0.0;
		i1 = indr[i];
		for (j=1;j<=nc;j++)
		{
			j1 = indc[j];
			if ((mat[i1][j1]==1 && d[i][j]<0) || (mat[i1][j1]==0 && d[i][j]>0))
			{
				unex += fabs(d[i][j]);
			}
		}
		unex /= nc;
		tIdiosync = 100*unex/UMAX;
		//fprintf(out,"%10.5f",tIdiosync);
	}
	//fprintf(out,"\n\n");

	//	Idiosyncratic temperature for columns:
	//fprintf(out,"Idiosyncratic temperature for columns:\n");
	//fprintf(out,"Column:                    ");
	for (j=1;j<=nc;j++);// fprintf(out,"%10i",j);
	//fprintf(out,"\n");
	//fprintf(out,"Idiosyncratic temperature: ");
	for (j=1;j<=nc;j++)
	{
		unex = 0.0;
		j1 = indc[j];
		for (i=1;i<=nr;i++)
		{
			i1 = indr[i];
			if ((mat[i1][j1]==1 && d[i][j]<0) || (mat[i1][j1]==0 && d[i][j]>0))
			{
				unex += fabs(d[i][j]);
			}
		}
		unex /= nr;
		tIdiosync = 100*unex/UMAX;
		//fprintf(out,"%10.5f",tIdiosync);
	}
	//fprintf(out,"\n\n");

	//close_file(out);
	//fclose(out);
}
/***
 *****************************************
 ***/
void permute(long& idum, int n, int index[])
//	Takes an array index[i] as imput, with 1=1,..., n and permutes it at random
{
	int i,j,k,m=n;
	int *perm;
	perm = ivector(1,n);
	for (i=1;i<=n;i++)
	{
		j = 1 + int(m*ran1(idum));
		if (j<=m)
		{
			perm[i] = index[j];
			for (k=j;k<m;k++) index[k] = index[k+1];
			m--;
		}
		else i--;
	}
	for (i=1;i<=n;i++) index[i] = perm[i];
	free_ivector(perm,1,n);
}
/***
 *****************************************
 ***/
void mutate(long& idum, int n, int index[])
//	Takes an array index[i] as imput, with 1=1,..., n and returns a random mutation
//	Mutations consist of a cyclic permutation of a part of the sequence
{
	int i,i1,i2,itoto;

	i1 = 1 + int(n*ran1(idum));
	i2 = 1 + int(n*ran1(idum));
	if (i1 < i2)
	{
		itoto = index[i2];
		for (i=i2;i>i1;i--) index[i] = index[i-1];
		index[i1] = itoto;
	}
	else if (i1 > i2)
	{
		itoto = index[i2];
		for (i=i2;i<i1;i++) index[i] = index[i+1];
		index[i1] = itoto;
	}
}
/***
 *****************************************
 ***/
void choosePlayers(long& idum, int n,int m,int arr[])
//	Chooses n numbers from integers 1, 2... m, and returns them as the elements of arr[]
{
	int i,j,k;
	int *index;

	index = ivector(1,m);
	if (n>m) nrerror();
	else if (n==m)
	{
		for (i=1;i<=n;i++) arr[i] = i;
	}
	else
	{
		for (i=1;i<=m;i++) index[i] = i;
		for (i=1;i<=n;i++)
		{
			j = 1 + int(m*ran1(idum));
			if (j<=m)
			{
				arr[i] = index[j];
				for (k=j;k<m;k++) index[k] = index[k+1];
				m--;
			}
			else i--;
		}
	}
	free_ivector(index,1,n);
}
/***
 *****************************************
 ***/
void crossOver(long& idum, int nr, int nc, int pr[], int pc[], int or1[], int oc[])
//	At input, the arrays pr[] and or1[] contain a permutation of row indexes from two individuals of the population,
//	while pc[] and oc[] contain permutations of the column indexes. At output, pr and pc are not changed, but or and oc
//	contain new permutations, obtained from "crossing" over the corresponding input strategies, and possibly by mutating 
//	the resulting strategy.
{
	int i,j,k,m;
	int flag=0;
	int *ind;

	while (flag==0)
	{
		if (ran1(idum)<0.5)
		{
			//	crossing over of rows
			ind = ivector(1,nr);
			for (i=1;i<=nr;i++) ind[i] = 0;
			k = 2 + int((nr-2)*ran1(idum));
			for (i=1;i<=k;i++) ind[or1[i]] = 1;
			for (i=k+1;i<=nr;i++)
			{
				if (ind[pr[i]]==0)
				{
					ind[pr[i]] = 1;
					or1[i] = pr[i];
				}
				else or1[i] = 0;
			}
			m = 0;
			for (i=1;i<=nr;i++)
			{
				if (ind[i]==0)
				{
					m++;
					ind[m] = i;
				}
			}
			if (m>1) permute(idum,m,ind);
			if (m>0)
			{
				for (i=1;i<=nr;i++)
				{
					if (or1[i]==0)
					{
						if (m<=0) nrerror();
						or1[i] = ind[m];
						m--;
					}
				}
			}
			free_ivector(ind,1,nr);
			flag = 1;
		}
		if (ran1(idum)<0.5)
		{
			//	crossing over of columns
			ind = ivector(1,nc);
			for (j=1;j<=nc;j++) ind[j] = 0;
			k = 2 + int((nc-2)*ran1(idum));
			for (j=1;j<=k;j++) ind[oc[j]] = 1;
			for (j=k+1;j<=nc;j++)
			{
				if (ind[pc[j]]==0)
				{
					ind[pc[j]] = 1;
					oc[j] = pc[j];
				}
				else oc[j] = 0;
			}
			m = 0;
			for (j=1;j<=nc;j++)
			{
				if (ind[j]==0)
				{
					m++;
					ind[m] = j;
				}
			}
			if (m>1) permute(idum,m,ind);
			if (m>0)
			{
				for (j=1;j<=nc;j++)
				{
					if (oc[j]==0)
					{
						if (m<=0) nrerror();
						oc[j] = ind[m];
						m--;
					}
				}
			}
			free_ivector(ind,1,nc);
			flag = 1;
		}
	}
	//	now mutate:
	if (ran1(idum)<0.1) mutate(idum,nr,or1);
	if (ran1(idum)<0.1) mutate(idum,nc,oc);
}
/***
 *****************************************
 ***/
void indexx(int n, int arr[], int indx[])
{
	int i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	int a;
   
	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M_IND) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1])
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) nrerror();
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
/***
 *****************************************
 ***/
void indexxD(int n, double arr[], int indx[])
{
	int i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	double a;
   
	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M_IND) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1])
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) nrerror();
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
/***
 *****************************************
 ***/
double zbrent(int nr, int nc, double u1, double v1, double z)
{
	int iter;
	double x1,x2=1,tol=1.0e-5;
	double a,b=x2,c=x2,d=0.0,e=0.0,min1,min2;
	double fa,fb,fc,p,q,r,s,tol1,xm;

	x1 = (u1+v1<1 ? 0 : u1+v1-1);
	a = x1;
	fa=func(a,nr,nc,u1,v1,z);
	fb=func(b,nr,nc,u1,v1,z);
   
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
	{
		nrerror();
	}
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPSB*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=func(b,nr,nc,u1,v1,z);
	}
	nrerror();
	return 0.0;
}
/***
 *****************************************
 ***/
double func(double yy, int nr, int nc, double x1, double y1, double z)
{
	double f,s,x;
	x = nc*(x1+y1-(nr-1)*yy/nr-0.5/nc-0.5/nr)/(nc-1.0);
	if (fabs(x)<EPS) s = 1;
	else if (x < 1) s = pow(1-x,z);
	else s = 0;
	if (fabs(1-yy)<EPS) f = s;
	else if (yy>0) f = pow(yy,z) + s - 1;
	else f = s - 1;
	return f;
}
/***
 *****************************************
 ***/
double ran1(long& idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
   
	if (idum <= 0 || !iy) {
		if (-idum < 1) idum=1;
		else idum = -idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
/***
 *****************************************
 ***/
void avevar(double data[], unsigned long n, double& ave, double& var)
{
	unsigned long j;
	double s,ep;
   
	for (ave=0.0,j=1;j<=n;j++) ave += data[j];
	ave /= n;
	var=ep=0.0;
	for (j=1;j<=n;j++) {
		s=data[j]-ave;
		ep += s;
		var += s*s;
	}
	var=(var-ep*ep/n)/(n-1);
}
/***
 *****************************************
 ***/
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;
   
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror();
	return v-nl+NR_END;
}
/***
 *****************************************
 ***/
double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
   
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror();
	return v-nl+NR_END;
}
/***
 *****************************************
 ***/
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
/***
 *****************************************
 ***/
void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
/***
 *****************************************
 ***/
int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
   
	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror();
	m += NR_END;
	m -= nrl;
   
	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror();
	m[nrl] += NR_END;
	m[nrl] -= ncl;
   
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
   
	/* return pointer to array of pointers to rows */
	return m;
}
/***
 *****************************************
 ***/
double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
   
	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror();
	m += NR_END;
	m -= nrl;
   
	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror();
	m[nrl] += NR_END;
	m[nrl] -= ncl;
   
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
   
	/* return pointer to array of pointers to rows */
	return m;
}
/***
 *****************************************
 ***/
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an integer matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
/***
 *****************************************
 ***/
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
/***
 *****************************************
 ***/
void nrerror()
/* Numerical Recipes standard error handler */
{
	//fprintf(stderr,"\n");
	//fprintf(stderr,"\n");
	//fprintf(stderr,"Run-time error...\n");
	//fprintf(stderr,"\n");
	//fprintf(stderr,"fehler\n");
	//fprintf(stderr,"\n\n");
	//fprintf(stdout, "Type a number and press ENTER to continue\n");
	//fprintf(stderr,"\n");
//        cin >> TOURSIZE;

	exit(1);
}
/***
 *****************************************
 ***/

}

   
