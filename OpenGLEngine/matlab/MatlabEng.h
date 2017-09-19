#ifndef _MATLAB_ENGINE_H_
#define _MATLAB_ENGINE_H_

//To use this project, first make sure lib is compatible with the win32/x64
//Second in Debuging, Environment: Path=address_of_the dll;

//  Open 			Start up MATLAB engine
//  Close 			Shut down MATLAB engine
//  GetVariable 	Get a MATLAB array from the MATLAB engine
//  PutVariable 	Send a MATLAB array to the MATLAB engine
//  EvalString 		Execute a MATLAB command
//  OutputBuffer 	Create a buffer to store MATLAB text output
//  OpenSingleUse 	Start a MATLAB engine session for single, nonshared use 
//  GetVisible 		Determine visibility of MATLAB engine session
//  SetVisible 		Show or hide MATLAB engine session

#ifndef FALSE
#define FALSE	0
#endif

#ifndef TRUE
#define TRUE	1
#endif

#include "MatlabInclude\matrix.h"
#include "MatlabInclude\engine.h"

#include "..\mathutilities\SparseMatrix.h"

#if defined(_M_X64) || defined(__amd64__)
#define WINDOWS_X64
#endif

#ifdef WINDOWS_X64
#pragma comment(lib, "./win64version/libeng.lib")
#pragma comment(lib, "./win64version/libmx.lib")
#pragma comment(lib, "./win64version/libmat.lib")
#pragma comment(lib, "./win64version/glfw3.lib")
#else
#pragma comment(lib, "./win32version/libeng.lib")
#pragma comment(lib, "./win32version/libmx.lib")
#pragma comment(lib, "./win32version/libmat.lib")
#pragma comment(lib, "./win32version/glfw3.lib")
#endif 


class CMatlabEng  
{
public:
	
	int OutputBuffer(char *p, int n);
	/*
	Purpose			Specify buffer for MATLAB output
	Arguments		n
					Length of buffer p.
					p
					Pointer to character buffer of length n.
	Description		OutputBuffer defines a character buffer for engEvalString to 
					return any output that ordinarily appears on the screen.
					The default behavior of EvalString is to discard any standard 
					output caused by the command it is executing. 
					OutputBuffer(ep,p,n) tells any subsequent calls to 
					EvalString to save the first n characters of output in the 
					character buffer pointed to by p. 
					To turn off output buffering, use OutputBuffer(ep,NULL,0);
	*/

	void OpenSingleUse(const char *startcmd, void *dcom, int *retstatus);
	/*
	Purpose			Start a MATLAB engine session for single, nonshared use
	Arguments		startcmd
					String to start MATLAB process. 
					On Windows, the startcmd string must be NULL.
					dcom
					Reserved for future use; must be NULL.
					retstatus
					Return status; possible cause of failure.
	Description		This routine allows you to start multiple MATLAB processes for 
					the purpose of using MATLAB as a computational engine. 
					OpenSingleUse starts a MATLAB process, establishes a connection, 
					and returns a unique engine identifier, or NULL if the 
					open fails. OpenSingleUse starts a new MATLAB process each time 
					it is called.
					OpenSingleUse opens a COM channel to MATLAB. This starts the
					MATLAB that was registered during installation. 
					If you did not register during installation, on the command line 
					you can enter the command: 
						matlab /regserver
					OpenSingleUse allows single-use instances of a MATLAB engine 
					server. OpenSingleUse differs from Open, which allows multiple 
					users to use the same MATLAB engine server.
	*/

	int GetVisible(bool* value);
	/*
	Purpose			Copy a variable from a MATLAB engine’s workspace
	Arguments		name
					Name of mxArray to get from MATLAB.
	Description		return status of the window for the MATLAB engine session, 
					is visible or invisible on the Windows desktop
	Returns			SetVisible returns 0 on success, and 1 otherwise.
	*/

	int SetVisible(bool value);
	/*
	Purpose			Show or hide MATLAB engine session
	Arguments		value
					Value to set the Visible property to. Set value to 1 to make 
					the engine window visible, or to 0 to make it invisible.
	Description		SetVisible makes the window for the MATLAB engine session, 
					either visible or invisible on the Windows desktop. 
					You can use this function to enable or disable user interaction 
					with the MATLAB engine session.
	
	Returns			SetVisible returns 0 on success, and 1 otherwise.
	*/
	
	mxArray* GetVariable(const char* name);
	/*
	Purpose			Copy a variable from a MATLAB engine’s workspace 
	Arguments		name
					Name of mxArray to get from MATLAB.
	Description		reads the named mxArray from the MATLAB engine session
					associated with ep and returns a pointer to a newly allocated 
					mxArray structure, or NULL if the attempt fails. GetVariable 
					fails if the named variable does not exist.
					Be careful in your code to free the mxArray created by this 
					routine when you are finished with it.
	*/

	int PutVariable(const char *name, const mxArray *mp);
	/*
	Purpose			Put variables into a MATLAB engine’s workspace
	Arguments		name
					Name given to the mxArray in the engine’s workspace.
					mp
					mxArray pointer.
	Description		PutVariable writes mxArray mp to the engine ep, giving it 
					the variable name, name. If the mxArray does not exist in the 
					workspace, it is created. If an mxArray with the same name 
					already exists in the workspace, the existing mxArray is 
					replaced with the new mxArray.
					
	Returns			PutVariable returns 0 if successful and 1 if an error occurs.
	*/

	int EvalString(const char* string);
	/*
	Purpose			Evaluate expression in string
	Arguments		string 
					String to execute.
	Description		evaluates the expression contained in string for the MATLAB
					engine session, previously started by Open. 
	Returns			It returns a nonzero value if the MATLAB session is no 
					longer running, and zero otherwise.
	*/

	void Open(const char* StartCmd);
	/*
	Purpose			Start a MATLAB engine session
	Arguments		Startcmd
					String to start MATLAB process. On Windows, the startcmd 
					string must be NULL.
	Description		This routine allows you to start a MATLAB process for the 
					purpose of using MATLAB as a computational engine.
	*/

	int Close();
	/*
	Close			Quit a MATLAB engine session
	Arguments		-
	Description		This routine allows you to quit a MATLAB engine session.
	Returns			Close sends a quit command to the MATLAB engine session and 
					closes the connection. It returns 0 on success, and 1 otherwise. 
					Possible failure includes attempting to terminate a MATLAB 
					engine session that was already terminated.
	*/

	CMatlabEng();
	virtual ~CMatlabEng();
protected:
	Engine* pEng;

public:
	void PutVar_sparse(const char *name, CSparseMatrix *ms);
	void GetVar_d(const char* name, double * sm);
	int PutVar_d(const char *name, double *ms, int row);

};

#endif // _MATLAB_ENGINE_H_

/*
	//test of output matrix from matlab
	//matlab.EvalString("T=[1 4; 3 0; 0 5]");
	//SMatrix *m1 = NULL;
	//m1 = matlab.GetVar("T");
	//m1->ShowDenseMatrix();

	//test of input matrix to matlab
	//matlab.PutVar("U",m1);
	//mxArray *U = NULL;
	//U = mxCreateDoubleMatrix(3, 2, mxREAL);	
	//U = matlab.GetVariable("U");
	//double* pU = mxGetPr(U);
	//for(int i = 0; i < 6; i++)
	//	cout<<pU[i]<<endl;

	//SPAESE - test of output matrix from matlab
	matlab.EvalString("A=[0 11 22 33; 987 0 0 0; 3 2 0 3]");
	matlab.EvalString("T=sparse(A)");
	SMatrix *m1 = NULL;
	m1 = matlab.GetSparseVar("T");
	m1->ShowDenseMatrix();
	m1->ShowMatrix();


	//SPAESE - test of input matrix to matlab
	matlab.PutSparseVar("U",m1);
	SMatrix *m2 = NULL;
	m2 = matlab.GetSparseVar("U");
	m2->ShowDenseMatrix();
	m2->ShowMatrix();

	//SPAESE - test of calculation in matlab
	matlab.EvalString("B=[0 11 0; 987 0 0; 3 2 0;3 2 5]");
	matlab.EvalString("S=sparse(B)");
	matlab.EvalString("C=T*S");
	SMatrix *m3 = NULL;
	m3 = matlab.GetSparseVar("C");
	m3->ShowDenseMatrix();
	m3->ShowMatrix();

	delete m1,m2;
	*/
