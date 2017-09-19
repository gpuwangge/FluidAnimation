// MatlabEngine.cpp: implementation of the CMatlabEng class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
#include "MatlabEng.h"
#pragma warning(disable : 4267) 
//////////////////////////////////////////////////////////////////////
// Construction/Destru=ction
//////////////////////////////////////////////////////////////////////

CMatlabEng::CMatlabEng()
{
    pEng=NULL;
}

CMatlabEng::~CMatlabEng()
{
    if (pEng!=NULL)
        Close();
}

void CMatlabEng::Open(const char *StartCmd)
{
    pEng=engOpen(StartCmd);
}

int CMatlabEng::Close()
{
    int Result=engClose(pEng);
    if (Result==0)	//Success
        pEng=NULL;

    return Result;
}

int CMatlabEng::EvalString(const char *string)
{
    return (engEvalString(pEng, string));
}

mxArray* CMatlabEng::GetVariable(const char *name)
{
    return (engGetVariable(pEng, name));
}


int CMatlabEng::GetVisible(bool* value)
{
    return (engGetVisible(pEng, value));
}

void CMatlabEng::OpenSingleUse(const char *startcmd, void *dcom, int *retstatus)
{
    pEng=engOpenSingleUse(startcmd, dcom, retstatus);
}

int CMatlabEng::OutputBuffer(char *p, int n)
{
    return (engOutputBuffer(pEng, p, n));
}

int CMatlabEng::PutVariable(const char *name, const mxArray *mp)
{
    return (engPutVariable(pEng, name, mp));
}

int CMatlabEng::SetVisible(bool value)
{
    return (engSetVisible(pEng, value));
}





void CMatlabEng::PutVar_sparse(const char *name, CSparseMatrix *ms) {
	int col = ms->numCols;
	int row = ms->numRows;
	int num = ms->GetElementNumber();
	int max_col_num = ms->GetMaxColNumer();

	mxArray *A = mxCreateSparse(row, col, num, mxREAL);
	double *pA = mxGetPr(A);
	mwIndex *adr_i = mxGetIr(A);
	mwIndex *adr_j = mxGetJc(A);

	int index = 0;//
	CMatrixElement * theElm;
	double * col_elem;
	int *  col_id;
	int col_count;
	col_elem = new double[max_col_num + 10];
	col_id = new int[max_col_num + 10];
	for (int i = 0; i < col; i++)
	{

		adr_j[i] = index;
		col_count = 0;
		for (theElm = ms->colList[i]; theElm != NULL; theElm = theElm->colNext)
		{
			col_elem[col_count] = theElm->value;
			col_id[col_count] = theElm->i;//
			col_count++;
		} //end of

		for (int j = 0; j < col_count; j++)
		{
			int test = col_id[j];
			pA[index] = col_elem[col_count - j - 1];
			adr_i[index] = col_id[col_count - j - 1];
			//pA[index] = col_elem[j];
			//adr_i[index] = col_id[j];
			index++;
		} //end of
	}
	adr_j[col] = index;
	delete[] col_id;
	delete[] col_elem;
	PutVariable(name, A);
	mxDestroyArray(A);

	////////////////////////////////////////////////
	char makeitsparse[1024];
	sprintf_s(makeitsparse, "tmp__ = sparse(%s); %s = tmp__;", name, name);
	this->EvalString(makeitsparse);
}

void CMatlabEng::GetVar_d(const char* name, double * sm) {
	mxArray *A = GetVariable(name);
	double* pA = mxGetPr(A);
	int col = mxGetN(A); // should be 1
	int row = mxGetM(A); // should be n
						 //CSparseMatrix *sm = new CSparseMatrix(row,col);

	int index = 0;
	for (int i = 0; i < col; i++)
		for (int j = 0; j < row; j++) {
			sm[index] = (double)pA[index];
			index++;
		}
	mxDestroyArray(A);
}

int CMatlabEng::PutVar_d(const char *name, double *ms, int row) {
	mxArray *A = NULL;
	int col = 1;
	A = mxCreateDoubleMatrix(row, col, mxREAL);
	double* pA = mxGetPr(A);
	int index = 0;
	for (int i = 0; i < col; i++)
		for (int j = 0; j < row; j++)
			pA[index++] = ms[j];

	PutVariable(name, A);
	int ret = PutVariable(name, A);
	mxDestroyArray(A);
	return ret;
}