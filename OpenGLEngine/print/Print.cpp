/***Created by Xiaojun***/
//#include "StdAfx.h"
#include "Print.h"

#include <iomanip>
using namespace std;

CPrint::CPrint(void){}
CPrint::~CPrint(void){}

void CPrint::ClearFile(std::string filename){
	fs.open(filename.c_str(), std::fstream::out|std::fstream::binary);
}

void CPrint::SaveFile_EndMark(std::string filename){
	fs.open(filename.c_str(), std::fstream::app|std::fstream::binary);
	fs<<"-99999\n";
}

void CPrint::OpenForSave_Binary(std::string filename, bool isApp){
	//fs.open(filename.c_str(), std::fstream::in);
	if(isApp)
		fs.open(filename.c_str(), std::fstream::app|std::fstream::binary);
	else
		fs.open(filename.c_str(), std::fstream::out|std::fstream::binary);
}

void CPrint::SaveFile_Add(int size, double *x, int frameNum){////////
	fs<<frameNum;
	fs<<" ";
	fs<<size;
	fs<<"\n";
	for(int i = 0; i < size; i++){
		fs<<x[i];
		fs<<"\n";
	}
}

void CPrint::WriteBinaryFile_Add(int size, double *x, int frameNum){////////
	fs.write((char*)&frameNum, sizeof(int));
	fs.write((char*)&size, sizeof(int));
	for(int i = 0; i < size; i++){
		fs.write((char*)&(x[i]), sizeof(double));
	}
}

void CPrint::WriteBinaryFile_Add(int size, double *x){////////
	for(int i = 0; i < size; i++){
		fs.write((char*)&(x[i]), sizeof(double));
	}
}



bool CPrint::ReadFile_One_Frame(double **pos, double **posPrev, double **posPrev2, double **acc){
	int frameNum;
	int size;

	fs>>frameNum;
	fs>>size;

	//if(size == EOF)//how to know if this is the end of file?
		//return false;

	if((*pos)!= NULL)
		delete *pos;
	*pos = new double[size];

	if((*posPrev)!= NULL)
		delete *posPrev;
	*posPrev = new double[size];

	if((*posPrev2)!= NULL)
		delete *posPrev2;
	*posPrev2 = new double[size];

	if((*acc)!= NULL)
		delete *acc;
	*acc = new double[size];

	for(int i = 0; i < (size / 3); i++){
			fs>>(*pos)[3 * i + 0];
			fs>>(*pos)[3 * i + 1];
			fs>>(*pos)[3 * i + 2];
			//fs>>(*acc)[3 * i + 0];
			//fs>>(*acc)[3 * i + 1];
			//fs>>(*acc)[3 * i + 2];

			/*if(i > 0){
				(*posPrev)[3 * i + 0] = (*pos)[3 * (i-1) + 0];
				(*posPrev)[3 * i + 1] = (*pos)[3 * (i-1) + 1];
				(*posPrev)[3 * i + 2] = (*pos)[3 * (i-1) + 2];
			}else{
				(*posPrev)[3 * i + 0] = 0;
				(*posPrev)[3 * i + 1] = 0;
				(*posPrev)[3 * i + 2] = 0;
			}

			if(i > 1){
				(*posPrev2)[3 * i + 0] = (*pos)[3 * (i-2) + 0];
				(*posPrev2)[3 * i + 1] = (*pos)[3 * (i-2) + 1];
				(*posPrev2)[3 * i + 2] = (*pos)[3 * (i-2) + 2];
			}else{
				(*posPrev2)[3 * i + 0] = 0;
				(*posPrev2)[3 * i + 1] = 0;
				(*posPrev2)[3 * i + 2] = 0;
			}*/
	}

	return true;

}


bool CPrint::ReadBinaryFile_One_Frame(double **pos){
	int frameNum;
	int size;

	fs.read((char*)&frameNum, sizeof(int));
	fs.read((char*)&size, sizeof(int));

	if((*pos)!= NULL) delete *pos;
	*pos = new double[size];

	for(int i = 0; i < (size / 3); i++){
		fs.read((char*)&((*pos)[3 * i + 0]), sizeof(double));
		fs.read((char*)&((*pos)[3 * i + 1]), sizeof(double));
		fs.read((char*)&((*pos)[3 * i + 2]), sizeof(double));
	}

	return true;
}


double CPrint::ReadASCIIFile_One_Number(){
	double data;
	//fs.read((char*)&data, sizeof(double));
	fs>>data;
	return data;
}

double CPrint::ReadBinaryFile_One_Number(){
	double data;
	fs.read((char*)&data, sizeof(double));
	return data;
}


int CPrint::ReadFile_Frame(std::string filename, double **x, int maxFrameNum){
	fs.open(filename.c_str(), std::fstream::in|std::fstream::binary);
	int frameNum;
	int size;
	int returnNum = 0;

	fs>>frameNum;
	fs>>size;

	if((*x)!=NULL)
		delete *x;
	*x = new double[size * maxFrameNum];


	for(int j  = 0; j < maxFrameNum; j++){

		for(int i = 0; i < size; i++){
			fs>>(*x)[j * size + i];
		}
		fs>>frameNum;
		returnNum++;
		if(frameNum == -99999)
			break;

		fs>>size;

	}

	return returnNum;
}

void CPrint::ReadFile(std::string filename, double **x){
	fs.open(filename.c_str(), std::fstream::in|std::fstream::binary);
	int size;
	fs>>size;
	if((*x)!=NULL)
		delete *x;
	*x = new double[size];
	for(int i = 0; i < size; i++){
		fs>>(*x)[i];
	}
}

void CPrint::SaveFile(std::string filename, int size, double *x){
	fs.open(filename.c_str(), std::fstream::out|std::fstream::binary);
	fs<<size;
	for(int i = 0; i < size; i++){
		fs<<"\n";
		fs<<x[i];
	}
}

void CPrint::CreateFile_ASCII(std::string filename){
	fs.open(filename.c_str(), std::fstream::out);
}

void CPrint::SaveFile_ASCII(std::string filename, int size, double *x, std::string *s){
	fs.open(filename.c_str(), std::fstream::out);
	for(int i = 0; i < size; i++){
		for(int j = 0; j < s[i].length(); j++) fs<<s[i][j];
		fs<<" ";
		fs<<x[i];
		fs<<"\n";
	}
	CloseFile();
}

void CPrint::ReadFile_ASCII(std::string filename, double *x){
	fs.open(filename.c_str(), std::fstream::in|std::fstream::binary);
	char c;
	int i = 0;
	while(!fs.eof()){
		do{ fs>>c;}while(c != ':');
		fs>>x[i++];
	}
	CloseFile();
}


void CPrint::CreateFile(std::string filename){
	fs.open(filename.c_str(), std::fstream::out|std::fstream::binary);
}

void CPrint::CreateFileAdd(std::string filename){
	fs.open(filename.c_str(), std::fstream::app|std::fstream::out|std::fstream::binary);
}



//void CPrint::PrintToFile(CSparseMatrix *x){
//	CMatrixElement *theElem;
//	for(int i = 0; i < x->numRows; i++){
//		PrintToFile("(");
//		PrintToFile(i);
//		PrintToFile_Space(")");
//		for(theElem = x->rowList[i]; theElem != NULL; theElem = theElem->rowNext)
//			if(theElem->value != 0){
//				PrintToFile(theElem->value);
//				PrintToFile("(");
//				PrintToFile(theElem->j);
//				PrintToFile_Space(")");
//			}
//		PrintToFile_Enter("");
//	}
//}

void CPrint::PrintToFile(double x){
	fs<<setiosflags(ios::fixed)<<setw(25)<<setprecision(15)<<x;
}
void CPrint::PrintToFile_Space(double x){
	PrintToFile(x);
	fs<<"  ";
}
void CPrint::PrintToFile_Enter(double x){
	PrintToFile(x);
	fs<<"\r\n";
	/*Binary*/
	//char result[sizeof(double)];
	//memcpy(result,&x,sizeof(x));
	//fs.write(result, sizeof(result));
}

void CPrint::PrintToFile(double x, int p){
	fs<<setiosflags(ios::fixed)<<setw(p*3)<<setprecision(p)<<x;
}
void CPrint::PrintToFile_Space(double x, int p){
	PrintToFile(x,p);
	fs<<"  ";
}
void CPrint::PrintToFile_Enter(double x, int p){
	PrintToFile(x,p);
	fs<<"\r\n";
}

void CPrint::PrintToFile_ASCII(char x){
	fs<<x;
}

void CPrint::PrintToFile(char *x){
	//fs<<setiosflags(ios::fixed)<<x;
	fs<<x;
}
void CPrint::PrintToFile_Space(char *x){
	PrintToFile(x);
	fs<<"  ";
}
void CPrint::PrintToFile_Enter(char *x){
	//fs.write(x, sizeof(x));
	PrintToFile(x);
	fs<<"\r\n";
}

void CPrint::PrintToFile_Compact(int x){
	fs<<x;
}

void CPrint::PrintToFile(int x){
	fs<<setiosflags(ios::fixed)<<setw(10)<<setprecision(5)<<x;
}
void CPrint::PrintToFile_Space(int x){
	PrintToFile(x);
	fs<<"  ";
}
void CPrint::PrintToFile_Enter(int x){
	PrintToFile(x);
	fs<<"\r\n";
}
	


void CPrint::PrintToFile_Bracket(double x){
	fs<<"(";
	fs<<x;
	fs<<") ";
}

//void CPrint::PrintToFile_Space(vec3 x){
//	fs<<x.n[0];
//	fs<< " ";
//	fs<<x.n[1];
//	fs<< " ";
//	fs<<x.n[2];
//	fs<<"  ";
//}
//void CPrint::PrintToFile_Enter(vec3 x){
//	fs<<x.n[0];
//	fs<< " ";
//	fs<<x.n[1];
//	fs<< " ";
//	fs<<x.n[2];
//	fs<<"\n";
//}

void CPrint::CloseFile(){
	fs.close();
}

void CPrint::OpenForRead(std::string filename){
	fs.open(filename.c_str(), std::fstream::in|std::fstream::binary);
}

void CPrint::OpenForRead_AscII(std::string filename){
	fs.open(filename.c_str(), std::fstream::in);
}

//-------------------


void CPrint::WriteBinaryFile(){
	double nNum = 111.34595645643;
	//std::string str("123456");
 
	fs.write((char*)&nNum, sizeof(double));
	//fs.write(str.c_str(), sizeof(char) * (str.size()));
}

void CPrint::OpenBinaryFile(std::string filename){
	fs.open(filename.c_str(), std::fstream::in|std::fstream::binary);
}


double CPrint::ReadBinaryFile(){
	double nNum;
	fs.read((char*)&nNum, sizeof(double));
	return nNum;
}

void CPrint::ClearBinaryFile(std::string filename){
	fs.open(filename.c_str(), std::fstream::out|std::fstream::binary);
}