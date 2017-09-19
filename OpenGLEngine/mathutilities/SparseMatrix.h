// SparseMatrix.h: interface for the CSparseMatrix class.
//
//////////////////////////////////////////////////////////////////////
//new

#ifndef _SPARSE_MATRIX_
#define _SPARSE_MATRIX_

#include <stdio.h>
#include <assert.h>
//#include "Utility\App.h"
#include "BasicDataStructure.h"


//#include "stdafx.h"

//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif


// User-defined tolerancy
#define TOL 0.00005

class CMatrixElement {
public:
	int i;
	int j;
	double value;
	CMatrixElement *rowNext;
	CMatrixElement *colNext;

	CMatrixElement(int newi=0, int newj=0, double newValue=0.0){
		i = newi;
		j = newj;
		value = newValue;
		rowNext = colNext = 0;
	}

	//virtual ~CMatrixElement() { 
	~CMatrixElement() { 
		//if (rowNext != 0)  
		//	delete rowNext; 
	}
};

class CSparseMatrix{
//protected:
public:
	int numRows;
	int numCols;
	CMatrixElement* *rowList;
	CMatrixElement* *colList;
	double* diagonal;
public:
	CSparseMatrix(int nRows){
		numRows = numCols = 0;
		rowList = colList = NULL;
		diagonal = NULL;
		setDimensions(nRows);
	}

	CSparseMatrix(int nRows, int nCols){
		numRows = numCols = 0;
		rowList = colList = NULL;
		diagonal = NULL;
		setDimensions(nRows,nCols);
	}

	CSparseMatrix(int nRows, int numEl, int i[], int j[],  double vals[]){
		numRows = numCols = 0;
		rowList = colList = NULL;
		diagonal = NULL;
		setDimensions(nRows);
		setValues(numEl,i,j,vals);
	}

	~CSparseMatrix(){
		//CleanUp();
	}

	void CleanUp(){
		//Old Method
		//Only delete rows since the matrix element deletes rowNext
		//for(int i = 0; i < numRows; i++){
		//   if(rowList[i] != NULL)
		//   delete rowList[i];
		//}

		//New Method
		//need validation
		for(int rowIndex = 0;  rowIndex < numRows; rowIndex++){ 
			for(CMatrixElement *theElem = rowList[rowIndex]; theElem != NULL; ){
				CMatrixElement *tmp = theElem;
				theElem = theElem->rowNext;
				delete tmp;
			}
		}

	   if(rowList != NULL)  delete [] rowList;  rowList  = NULL;
	   if(colList != NULL)  delete [] colList;  colList  = NULL;
	   if(diagonal != NULL) delete [] diagonal; diagonal = NULL;
	}

	void CSparseMatrix::setValues(int numEl, int i[], int j[], double vals[]){
	   for(int idx = 0; idx < numEl; idx++)
		  set1Value(i[idx],j[idx],vals[idx]);
	}

void
CSparseMatrix::set1Mat3(int i, int j, XMat33 m){
	/*set1Value(i + 0, j + 0, m.v[0].n[0]);
	set1Value(i + 0, j + 1, m.v[0].n[1]);
	set1Value(i + 0, j + 2, m.v[0].n[2]);

	set1Value(i + 1, j + 0, m.v[1].n[0]);
	set1Value(i + 1, j + 1, m.v[1].n[1]);
	set1Value(i + 1, j + 2, m.v[1].n[2]);

	set1Value(i + 2, j + 0, m.v[2].n[0]);
	set1Value(i + 2, j + 1, m.v[2].n[1]);
	set1Value(i + 2, j + 2, m.v[2].n[2]);*/


	modify1Value(i + 0, j + 0, m(0,0));
	modify1Value(i + 0, j + 1, m(0,1));
	modify1Value(i + 0, j + 2, m(0,2));

	modify1Value(i + 1, j + 0, m(1,0));
	modify1Value(i + 1, j + 1, m(1,1));
	modify1Value(i + 1, j + 2, m(1,2));

	modify1Value(i + 2, j + 0, m(2,0));
	modify1Value(i + 2, j + 1, m(2,1));
	modify1Value(i + 2, j + 2, m(2,2));
}

bool CSparseMatrix::IsSymmetric(){
	double thresh = 0.00000000001f;
	for(int rowIndex = 0;  rowIndex < numRows; rowIndex++){ 
		for(CMatrixElement *theElem = rowList[rowIndex]; theElem != NULL; theElem = theElem->rowNext){
			double x = GetValue(theElem->j, theElem->i);
			if(abs(theElem->value - x) > thresh) 
				return false;
		}
	}
	return true;
}

void
CSparseMatrix::add1Mat3(int i, int j, XMat33 m, double sign){
	//modify1Value(i + 0, j + 0, GetValue(i + 0,j + 0) + m.v[0].n[0] * sign);
	//modify1Value(i + 0, j + 1, GetValue(i + 0,j + 1) + m.v[0].n[1] * sign);
	//modify1Value(i + 0, j + 2, GetValue(i + 0,j + 2) + m.v[0].n[2] * sign);

	//modify1Value(i + 1, j + 0, GetValue(i + 1,j + 0) + m.v[1].n[0] * sign);
	//modify1Value(i + 1, j + 1, GetValue(i + 1,j + 1) + m.v[1].n[1] * sign);
	//modify1Value(i + 1, j + 2, GetValue(i + 1,j + 2) + m.v[1].n[2] * sign);

	//modify1Value(i + 2, j + 0, GetValue(i + 2,j + 0) + m.v[2].n[0] * sign);
	//modify1Value(i + 2, j + 1, GetValue(i + 2,j + 1) + m.v[2].n[1] * sign);
	//modify1Value(i + 2, j + 2, GetValue(i + 2,j + 2) + m.v[2].n[2] * sign);
	modify1Value(i + 0, j + 0, GetValue(i + 0,j + 0) + m(0,0) * sign);
	modify1Value(i + 0, j + 1, GetValue(i + 0,j + 1) + m(0,1) * sign);
	modify1Value(i + 0, j + 2, GetValue(i + 0,j + 2) + m(0,2) * sign);

	modify1Value(i + 1, j + 0, GetValue(i + 1,j + 0) + m(1,0) * sign);
	modify1Value(i + 1, j + 1, GetValue(i + 1,j + 1) + m(1,1) * sign);
	modify1Value(i + 1, j + 2, GetValue(i + 1,j + 2) + m(1,2) * sign);

	modify1Value(i + 2, j + 0, GetValue(i + 2,j + 0) + m(2,0) * sign);
	modify1Value(i + 2, j + 1, GetValue(i + 2,j + 1) + m(2,1) * sign);
	modify1Value(i + 2, j + 2, GetValue(i + 2,j + 2) + m(2,2) * sign);
}

void
CSparseMatrix::add1Mat3Trans(int i, int j,XMat33 m, double sign){
	//modify1Value(i + 0, j + 0, GetValue(i + 0,j + 0) + m.v[0].n[0] * sign);
	//modify1Value(i + 0, j + 1, GetValue(i + 0,j + 1) + m.v[1].n[0] * sign);
	//modify1Value(i + 0, j + 2, GetValue(i + 0,j + 2) + m.v[2].n[0] * sign);

	//modify1Value(i + 1, j + 0, GetValue(i + 1,j + 0) + m.v[0].n[1] * sign);
	//modify1Value(i + 1, j + 1, GetValue(i + 1,j + 1) + m.v[1].n[1] * sign);
	//modify1Value(i + 1, j + 2, GetValue(i + 1,j + 2) + m.v[2].n[1] * sign);

	//modify1Value(i + 2, j + 0, GetValue(i + 2,j + 0) + m.v[0].n[2] * sign);
	//modify1Value(i + 2, j + 1, GetValue(i + 2,j + 1) + m.v[1].n[2] * sign);
	//modify1Value(i + 2, j + 2, GetValue(i + 2,j + 2) + m.v[2].n[2] * sign);

	modify1Value(i + 0, j + 0, GetValue(i + 0,j + 0) + m(0,0) * sign);
	modify1Value(i + 0, j + 1, GetValue(i + 0,j + 1) + m(1,0) * sign);
	modify1Value(i + 0, j + 2, GetValue(i + 0,j + 2) + m(2,0) * sign);

	modify1Value(i + 1, j + 0, GetValue(i + 1,j + 0) + m(0,1) * sign);
	modify1Value(i + 1, j + 1, GetValue(i + 1,j + 1) + m(1,1) * sign);
	modify1Value(i + 1, j + 2, GetValue(i + 1,j + 2) + m(2,0) * sign);

	modify1Value(i + 2, j + 0, GetValue(i + 2,j + 0) + m(0,2) * sign);
	modify1Value(i + 2, j + 1, GetValue(i + 2,j + 1) + m(1,2) * sign);
	modify1Value(i + 2, j + 2, GetValue(i + 2,j + 2) + m(2,2) * sign);
}

void CSparseMatrix::set1Value(int i, int j, double val){
   // Insertion in rows
   CMatrixElement *theElem = new CMatrixElement(i,j,val);
   theElem->rowNext = rowList[i];
   rowList[i] = theElem;

   // Insertion in columns
   theElem->colNext = colList[j];
   colList[j] = theElem;

   // If on the diagonal, store it for fast access
   if(i==j) diagonal[i] = val;
   
}

void
CSparseMatrix::modify1Value(int i, int j, double val) {
	CMatrixElement *theElem = GetElement(i,j);

	if(theElem == NULL) set1Value(i,j,val);
	else theElem->value = val;
}

void
CSparseMatrix::add1Value(int i, int j, double val)
{
	CMatrixElement *theElem = GetElement(i,j);

	if(theElem == NULL)
		set1Value(i,j,val);
	else
		theElem->value += val;
}

void
CSparseMatrix::setRow(int i, CMatrixElement *head)
{
    // Set it in the row
	rowList[i] = head;
	// And in the column (and diagonal)
	CMatrixElement *theElem;
	for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext){
		theElem->colNext = colList[theElem->j];
		colList[theElem->j] = theElem;
		if(i == theElem->j)
			diagonal[i] = theElem->value;
	}
}

void CSparseMatrix::setDimensions(int nRows){
	// Clean up anyway. Safer, since life is a jungle.
    CleanUp();
    numRows = nRows;
    numCols = nRows;
		/*
		if(rowList)
			delete [] rowList;
		if(colList)
			delete [] colList;
		if(diagonal)
			delete [] diagonal;*/
    rowList = new CMatrixElement*[numRows];
    colList = new CMatrixElement*[numCols];
    diagonal = new double[numRows];
    for(int k = 0; k < numRows; k++){
	   diagonal[k] = 0.;
	   rowList[k] = colList[k] = NULL;
    }	

}

void
CSparseMatrix::setDimensions(int nRows, int nCols)
{
	// Clean up anyway. Safer, since life is a jungle.
  CleanUp();
  numRows = nRows;
  numCols = nCols;
  rowList = new CMatrixElement*[numRows];
  colList = new CMatrixElement*[numCols];
  diagonal = new double[numRows];
  for(int k = 0; k < numRows; k++)
  {
   diagonal[k] = 0.;
   rowList[k] = NULL;
  }
	for(int l = 0; l < numCols; l++)
		colList[l] = NULL;
}

CMatrixElement*
CSparseMatrix::GetElement(int i, int j)
{
	CMatrixElement *theElem;
	for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
		if(theElem->j == j)
			return theElem;
	return NULL;
}

double
CSparseMatrix::GetValue(int i, int j)
{
	CMatrixElement *theElem;
	for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
		if(theElem->j == j)
			return theElem->value;
	return 0.0;
}

double
CSparseMatrix::diagonalElement(int i)
{
    assert(i < numRows);
	return diagonal[i];
}

void
CSparseMatrix::Print()
{
   CMatrixElement *theElem;
   for(int i = 0; i < numRows; i++)
	   for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
		   printf("i=%d, j=%d: %f\n", theElem->i, theElem->j, theElem->value);
}

void
CSparseMatrix::PrintMathematica(FILE *fp)
{
   int i, j;
   fprintf(fp,"m = {");
   for(i = 0; i < numRows; i++)
   {
	   fprintf(fp,"\n{");
	   for(j = 0; j < numCols; j++)
	   {
		   fprintf(fp,"%f",GetValue(i,j));
		   if(j != numCols-1) fprintf(fp,", ");
	   }
	   fprintf(fp,"}");
	   if(i != numRows-1) fprintf(fp,",");
   }
   fprintf(fp,"}\n\n");
}
void 
CSparseMatrix::PrintMathematica_wyz(FILE *fp)
{
	int i, j;
	fprintf(fp,"m = {");
	for(i = 0; i < numRows; i++)
	{
		fprintf(fp,"\n{");
		for(j = 0; j < numCols; j++)
		{
			double a = GetValue(i,j);
			if( abs(a) > 0.001 )
			{
				fprintf(fp,"%d,%d %f",i,j,GetValue(i,j));
				if(j != numCols-1) fprintf(fp,", ");
			}
		}
		fprintf(fp,"}");
		if(i != numRows-1) fprintf(fp,",");
	}
	fprintf(fp,"}\n\n");
}

void
PrintVectorMathematica(FILE *fp, double theVec[], int n)
{
   int i;
   fprintf(fp,"v = {");
   for(i = 0; i < n; i++)
   {
	   fprintf(fp,"%f",theVec[i]);
	   if(i != n-1) fprintf(fp,", ");
   }
   fprintf(fp,"}\n\n");
}


void
CSparseMatrix::multMatVec(double *src,
													double *dest)
{
	assert(src && dest);
	CMatrixElement *theElem = NULL;
	for(int i = 0; i < numRows; i++)
	{
		double sum = 0;
		for(theElem = rowList[i];
		    theElem != NULL;
		    theElem = theElem->rowNext)
			sum += theElem->value * src[theElem->j];
		dest[i] = sum;
	}
}

void
CSparseMatrix::multTransMatVec(double *src,
															 double *dest)
{
	assert(src && dest);
	double sum;

	CMatrixElement *theElem = NULL;
	for(int j = 0; j < numCols; j++)
	{
		sum = 0.0;
		for(theElem = colList[j]; theElem != NULL; theElem = theElem->colNext)
			sum += theElem->value * src[theElem->i];
		dest[j] = sum;
	}
}

////row version by xiaojun
//void
//CSparseMatrix::multTransMatVec2(float *src,
//															 float *dest)
//{
//	assert(src && dest);
//	float sum;
//
//	CMatrixElement *theElem = NULL;
//	for(int i = 0; i < numRows; i++)
//	{
//		sum = 0.0;
//		for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
//			sum += theElem->value * src[theElem->j];
//		dest[i] = sum;
//	}
//}

//void CSparseMatrix::multTransMatMat(){
//  // M = transpose(M)*M
//  CSparseMatrix* tempMat = new CSparseMatrix(numCols);
//  double *colVec = new double[numRows];
//  double *des_colVec = new double[numCols];
//
//  int i, j;
//  CMatrixElement *theElem;
//
//  for(j = 0; j < numCols; j++)
//  {
//		// initialize the result
//		for(i = 0; i < numCols; i++)
//			des_colVec[i] = 0.0;
//
//		// grab column j
//		for(i = 0; i < numRows; i++)
//			colVec[i] = 0.0;
//		for(theElem = colList[j]; theElem != NULL; theElem = theElem->colNext)
//			colVec[theElem->i] = theElem->value;
//
//		// compute a column of M
//		multTransMatVec(colVec, des_colVec);
//		
//		// store in tempMat
//		for(i = 0; i < numCols; i++)
//	     if(des_colVec[i] != 0.0)
//				 tempMat->set1Value(i,j,des_colVec[i]);
//  }
//	
//  // copy tempMat to this
//  setDimensions(numCols);
//  rowList = tempMat->rowList;
//  colList = tempMat->colList;
//  diagonal = tempMat->diagonal;
//  delete [] colVec;
//  delete [] des_colVec;
//	// delete tempMat; // must fix that ! -> memory leaks here
//}

void CSparseMatrix::multTransMatMat_yz(){
	// M = transpose(M)*M  <-> A * B
	CSparseMatrix* tempMat = new CSparseMatrix(numCols);

	int j;
	CMatrixElement *theElem;

#ifdef ___DEBUG_YZ
	FILE * fp;
	fp = fopen("mat_mup_log.txt", "w");
	fclose(fp);
#endif

	for (j = 0; j < numCols; j++)
	{
		theElem = colList[j];
		if (theElem == NULL)
			continue;

		for (CMatrixElement * theElem_e = colList[j]; theElem_e != NULL; theElem_e = theElem_e->colNext)  //D^T: col means A.row
		{
			int r_id = theElem_e->i;
			int tem_rid = theElem_e->j;
			for (CMatrixElement *theElem_r = rowList[r_id]; theElem_r != NULL; theElem_r = theElem_r->rowNext)
			{
#ifdef ___DEBUG_YZ
				fp = fopen("mat_mup_log.txt", "a+");
#endif

				int tem_cid = theElem_r->j;

#ifdef ___DEBUG_YZ
				fprintf(fp, "j=%d r_id=%d tem_rid=%d tem_cid=%d\n", j, r_id, tem_rid, tem_cid);
#endif

				double value = theElem_e->value * theElem_r->value;

#ifdef ___DEBUG_YZ
				fprintf(fp, "A[%d,%d] = %f, B[%d,%d] = %f, C[%d,%d] += %f\n", tem_rid, r_id, theElem_e->value, r_id, tem_cid, theElem_r->value, tem_rid, tem_cid, value);
#endif

				tempMat->add1Value(tem_rid, tem_cid, value);

#ifdef ___DEBUG_YZ
				fclose(fp);
#endif

			}
		}

	}
	// copy tempMat to this
	//Cleanup();
	setDimensions(numCols);
	rowList = tempMat->rowList;
	colList = tempMat->colList;
	diagonal = tempMat->diagonal;
	// delete tempMat; // must fix that ! -> memory leaks here
}

void
CSparseMatrix::AddMatrix(CSparseMatrix *mat)
{
	int i;
	CMatrixElement *theElem, *matElem;
	for(i = 0; i < numRows; i++)
	{
		for(matElem = mat->rowList[i]; matElem != NULL; matElem = matElem->rowNext)
		{
			theElem = GetElement(matElem->i,matElem->j);
			if(theElem == NULL)
			{
				theElem = new CMatrixElement(matElem->i,matElem->j,matElem->value);
				theElem->rowNext = rowList[i];
				rowList[i] = theElem;
				theElem->colNext = colList[theElem->j];
				colList[theElem->j] = theElem;
			}
			else
				theElem->value += matElem->value;
		}
	}
	for(i = 0; i < numRows; i++)
		diagonal[i] += mat->diagonal[i];
}

void
CSparseMatrix::ScaleRow(int i, double s)
{
   CMatrixElement *theElem;
   for(theElem = rowList[i]; theElem != NULL; theElem = theElem->rowNext)
	   theElem->value *= s;
   diagonal[i] *= s;
}

//***************************************
// preconditionedBiConjugateGradient
//***************************************
unsigned int 
preconditionedBiConjugateGradient(double x[],
																	double b[],
																	double tol,
																	const unsigned int iter_max)
{
	double *dr = new double[numRows];
	double *drb = new double[numRows];
	double *dp = new double[numRows];
	double *dpb = new double[numRows];
	double *dz = new double[numRows];
	double *dzb = new double[numRows];
	double *dAp = new double[numRows];
	double *dATpb = new double[numRows];

	assert(dr && drb && dp && dpb && dz && dzb && dAp && dATpb);
	double mag_r, mag_rOld, mag_pbAp, mag_Residual, Residual0, alpha, beta;

	multMatVec(x,dAp);
	mag_r = mag_Residual = Residual0 = 0.;
	int i = 0;
	for(i = 0; i < numRows; i++)
	{

		dr[i] = drb[i] = b[i] - dAp[i];
		dp[i] = dpb[i] = dz[i] = dzb[i] = dr[i]/diagonalElement(i);	// Simple preconditioning
		mag_r += drb[i] * dz[i];
		mag_Residual += dz[i] * dz[i];
		Residual0 += b[i]*b[i]/(diagonalElement(i)*diagonalElement(i));
	}
	
	mag_Residual = Residual0*100; // Force the first iteration anyway.
	if(Residual0 == 0)
		Residual0 = 1.;	// To make it work even if ||b|| = 0
	unsigned int nbIter = 0;
	while(mag_Residual > tol && nbIter < iter_max)
	{
		nbIter++;
		multMatVec(dp,dAp);
		multTransMatVec(dpb,dATpb);
		mag_pbAp = 0.0;
		for(i = 0; i < numRows; i++)
			mag_pbAp += dpb[i] * dAp[i];
		
		if(mag_pbAp == 0)
		{
			//fprintf(stderr,"OOOOOCH!!! (mag_pbAp==0)\n");
			//return 0;
		}

		if(mag_r == 0 && mag_pbAp == 0)
			alpha = 1;
		else
			alpha = mag_r / mag_pbAp;
		mag_rOld = mag_r;
		mag_r = 0.0;
		for(i = 0; i < numRows; i++)
		{
			x[i] += alpha * dp[i];
			dr[i] -= alpha * dAp[i];
			drb[i] -= alpha * dATpb[i];
			dz[i] = dr[i]/diagonalElement(i);
			dzb[i] = drb[i]/diagonalElement(i);
			mag_r += drb[i] * dz[i];
		}
		
		if(mag_rOld == 0)
		{
			//fprintf(stderr,"OOOOOCH!!! (mag_rOld==0)\n");
			//return 0;
		}
		
		if(mag_r == 0 && mag_rOld == 0)
			beta = 1.0;
		else
			beta = mag_r / mag_rOld;
		mag_Residual = 0.;
		for(i = 0; i < numRows; i++)
		{
			dp[i] = dz[i] + beta * dp[i];
			dpb[i] = dzb[i] + beta * dpb[i];
			mag_Residual += dz[i] * dz[i];
		}
	}
	delete [] dr;
	delete [] drb;
	delete [] dp;
	delete [] dpb;
	delete [] dz;
	delete [] dzb;
	delete [] dAp;
	delete [] dATpb;

	return nbIter;
}


void
CSparseMatrix::Debug()
{
	int i, j;
	fprintf(stderr,"M [ %i %i ] =\n",numRows,numCols);
	for(i = 0; i < numRows; i++)
	{
		for(j = 0; j < numCols; j++)
			fprintf(stderr,"%f ",GetValue(i,j));
		fprintf(stderr,"\n");
	}
}


int CSparseMatrix::GetElementNumber(){
	int count = 0;
	CMatrixElement *theElem = NULL;
	for (int j = 0; j < numCols; j++)
	{
		for (theElem = colList[j]; theElem != NULL; theElem = theElem->colNext)
			count++;
	}
	return count;
}


int CSparseMatrix::GetMaxColNumer(){
	int count = 0;
	int temp = 0;
	CMatrixElement *theElem = NULL;
	for (int j = 0; j < numCols; j++)
	{
		count = 0;
		for (theElem = colList[j]; theElem != NULL; theElem = theElem->colNext)
			count++;
		if (count > temp)
			temp = count;
	}
	return temp;
}

};

#endif // _SPARSE_MATRIX_

