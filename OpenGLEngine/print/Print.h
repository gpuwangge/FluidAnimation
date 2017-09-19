/***Created by Xiaojun***/
#ifndef __WXJ_PRI_H__
#define __WXJ_PRI_H__
//#include "Auxiliary.h"
#include "fstream"
//#include "App.h"


class CPrint
{	
//private:
public:
	std::fstream fs;
public:
	CPrint(void);
	~CPrint(void);
	void CreateFile_ASCII(std::string filename);
	void CreateFile(std::string filename);
	void CreateFileAdd(std::string filename);

	void PrintToFile(double x);
	void PrintToFile_Space(double x);
	void PrintToFile_Enter(double x);

	void PrintToFile_ASCII(char x);

	void PrintToFile(char *x);
	void PrintToFile_Space(char *x);
	void PrintToFile_Enter(char *x);

	void PrintToFile(double x, int p);
	void PrintToFile_Space(double x, int p);
	void PrintToFile_Enter(double x, int p);

	void PrintToFile(int x);
	void PrintToFile_Space(int x);
	void PrintToFile_Enter(int x);


	void PrintToFile_Compact(int x);

	void PrintToFile_Bracket(double x);
	//void PrintToFile_Enter(vec3 x);
	
	//void PrintToFile(CSparseMatrix *x);
	void ReadFile(std::string filename, double **x);
	int ReadFile_Frame(std::string filename, double **x, int frameNum);
	void SaveFile(std::string filename, int size, double *x);
	void SaveFile_Add(int size, double *x, int frameNum);
	void ClearFile(std::string filename);
	void SaveFile_EndMark(std::string filename);
	void OpenForRead(std::string filename);
	void OpenForRead_AscII(std::string filename);
	void OpenForSave_Binary(std::string filename, bool isApp);

	void SaveFile_ASCII(std::string filename, int size, double *x, std::string *s);
	void ReadFile_ASCII(std::string filename, double *x);

	bool ReadFile_One_Frame(double **pos, double **posPrev,double **posPrev2,double **acc);
	
	void CloseFile();

	void WriteBinaryFile();
	double ReadBinaryFile();
	void OpenBinaryFile(std::string filename);
	void ClearBinaryFile(std::string filename);
	void WriteBinaryFile_Add(int size, double *x, int frameNum);
	void WriteBinaryFile_Add(int size, double *x);
	bool ReadBinaryFile_One_Frame(double **pos);
	double ReadBinaryFile_One_Number();
	double ReadASCIIFile_One_Number();
};
#endif