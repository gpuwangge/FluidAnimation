#pragma once

#include "..\cmathhelper.h"
#include "..\mathutilities\BasicDataStructure.h"
#include "..\mathutilities\SparseMatrix.h"
#include "..\matlab\MatlabEng.h"
#include "..\print\Print.h"
#include "cmarchingtriangle.h"

class CLCS {
public:
	float integrat_T;
	float integrat_dt;

	//============================================
	glm::vec2 gridSize;
	vector<vector<glm::vec3>> phi;//same size of the grid
	float maxPhi; float minPhi;
	//============================================

	glm::vec2 FTLESize;
	vector<vector<glm::mat2x2>> gradPhi;//(-1,-1)
	vector<vector<float>> FTLE; float maxFTLE; float minFTLE;

	vector < vector<float>> testFTLE;

	//============================================
	glm::vec2 LCSSize;
	vector<vector<glm::vec2>> gradFTLE;//(-3,-3)
	vector<vector<glm::mat2x2>> hessian;
	vector<vector<glm::mat2x2>> eigenvectorOfHessian;
	float maxEigenvectorHessian1; float minEigenvectorHessian1;
	float maxEigenvectorHessian2; float minEigenvectorHessian2;

	vector<vector<glm::mat2x2>> testEigenvectorOfHessian;

	vector<vector<glm::vec2>> eigenvaluesOFHessian;

	vector<vector<float>> S; float maxS; float minS;
	vector<vector<float>> F; float maxF; float minF;
	vector<vector<float>> testScaler;

	vector<vector<glm::mat2x2>> tensor_vertex;
	vector<vector<glm::mat2x2>> eigenvectorOfTensor;
	float maxEigenvectorTensor1; float minEigenvectorTensor1;
	float maxEigenvectorTensor2; float minEigenvectorTensor2;
	vector<vector<glm::vec2>> eigenvaluesOfTensor;
	//============================================

	glm::vec2 rectSize;
	vector<CTriangle> triangles;
	//unordered_map<int, CTriangle> triangles;
	vector<CHalfedge> halfedges;

	CSparseMatrix *L_Alpha; //laplacian with new tensor (A)
	CSparseMatrix *L_B; //B = L_Id' * L_Id + w*D (B)
	CSparseMatrix *L_Id;
	//laplacian
	//CSparseMatrix *L_IdwD; 

	//CSparseMatrix *L_S;//for smooth S
	//============================================
	float weight_b;
	float smooth_h;
	float smooth_n;
	string currentDirectory;

	CMarchingTriangle marchingTriangle;

	vector<pair<float, float>> C;

	CLCS(int width, int height) {
		integrat_T = 0;
		integrat_dt = 0;

		//example: 13x7 grid
		gridSize.x = height; gridSize.y = width;
		phi.resize(gridSize.x, vector<glm::vec3>(gridSize.y)); //13x7: 0~12, 0~6

		//============================================

		FTLESize.x = height - 1; FTLESize.y = width - 1;
		gradPhi.resize(FTLESize.x, vector<glm::mat2x2>(FTLESize.y));//12x6: 0~11, 0~5
		FTLE.resize(FTLESize.x, vector<float>(FTLESize.y));//12x6

		testFTLE.resize(FTLESize.x, vector<float>(FTLESize.y));//12x6

		//for (int i = 0; i < width - 1; i++) testFTLE[2][i] = 100;
		float centerW = FTLESize.y / 2;
		float centerH = FTLESize.x / 2;
		for (int i = 0; i < FTLESize.x; i++) {
			for (int j = 0; j < FTLESize.y; j++) {
				float w = j - centerW;
				float h = i - centerH;
				float radius = sqrt(w*w + h*h);
				if (radius > 24 && radius < 26) testFTLE[i][j] = 100;
				//if (radius < 1) testFTLE[i][j] = 100;
			}
		}
		//To smooth the testFTLE

		//============================================
		LCSSize.x = height - 3; LCSSize.y = width - 3;

		gradFTLE.resize(LCSSize.x, vector<glm::vec2>(LCSSize.y));//11x5(10x4 valid): 1~10, 1~4 (use median integral, ignore the first row and column)
		hessian.resize(LCSSize.x, vector<glm::mat2x2>(LCSSize.y));  //11x5(10x4 valid)
		
		eigenvectorOfHessian.resize(LCSSize.x, vector<glm::mat2x2>(LCSSize.y)); //11x5(10x4 valid)

		testEigenvectorOfHessian.resize(LCSSize.x, vector<glm::mat2x2>(LCSSize.y)); //11x5(10x4 valid)
		for (int i = 0; i < LCSSize.x; i++) {
			for (int j = 0; j < LCSSize.y; j++) {
				float w = j - centerW;
				float h = i - centerH;
				float radius = sqrt(w*w + h*h);
				//if (radius > 26 && radius < 28) testFTLE[i][j] = 100;
				//if (radius < 1) testFTLE[i][j] = 100;
				if (radius < 0.0001) radius = 0.0001;
				if (radius > 24 && radius < 26) {
					testEigenvectorOfHessian[i][j] = glm::mat2x2(w / radius, h / radius,
																 h / radius, -w / radius);
				}
				else {
					testEigenvectorOfHessian[i][j] = glm::mat2x2(-w / radius, -h / radius,
																 -h / radius, w / radius);
				}
			}
		}

		eigenvaluesOFHessian.resize(LCSSize.x, vector<glm::vec2>(LCSSize.y)); //11x5(10x4 valid)

		S.resize(LCSSize.x, vector<float>(LCSSize.y)); //11x5(10x4 valid)
		F.resize(LCSSize.x, vector<float>(LCSSize.y)); //11x5(10x4 valid)
		
		testScaler.resize(LCSSize.x, vector<float>(LCSSize.y));

		tensor_vertex.resize(LCSSize.x, vector<glm::mat2x2>(LCSSize.y));//11x5(10x4 valid)
		eigenvectorOfTensor.resize(LCSSize.x, vector<glm::mat2x2>(LCSSize.y));
		eigenvaluesOfTensor.resize(LCSSize.x, vector<glm::vec2>(LCSSize.y));
		//============================================

		rectSize.x = height - 4; rectSize.y = width - 4;
		//each row has k=(width-4)*2 triangles.
		//rect: 9x3 (grid is still 11x5-10x4valid-ignore first row and column)
		triangles.resize(rectSize.x * rectSize.y * 2, CTriangle());  //9x3x2 (10x4 valid grid points, 9x3 rects, 2 triangle for each rect)
		halfedges.resize(rectSize.x * rectSize.y * 2 * 3, CHalfedge()); //9x3x2x3 (each triangle has 3 halfedge)	
		//============================================

		weight_b = 0;
		smooth_h = 0;
		smooth_n = 0;
		//char *buffer;
		//buffer = _getcwd(NULL, 0);
		currentDirectory = string(_getcwd(NULL, 0));

		

		ResetMinMax();
	}

	void ResetMinMax() {
		maxPhi = -99999;
		minPhi = 99999;

		maxEigenvectorHessian1 = -99999;  minEigenvectorHessian1 = 99999;
		maxEigenvectorHessian2 = -99999;  minEigenvectorHessian2 = 99999;

		maxEigenvectorTensor1 = -99999;  minEigenvectorTensor1 = 99999;
		maxEigenvectorTensor2 = -99999;  minEigenvectorTensor2 = 99999;

		maxS = -99999;
		minS = 99999;

		maxF = -99999;
		minF = 99999;

		maxFTLE = -99999;
		minFTLE = 99999;
	}

	/*void InitialTriangles(CGrid *grid) {
		for (int i = 0; i < triangles.size(); i++) {

		}
	}*/

	void InitialCotan_fromAngle(CGrid *grid) {
#ifdef DEBUG
		CPrint printer_cotan;
		string filename_cotan = "cotan_from_angle.txt";
		printer_cotan.CreateFile(OUTPUT_DATA + filename_cotan);
#endif
		//////////////////////////////////////////////////
		for (int i = 0; i < triangles.size(); i++) {
			float x0, y0, x1, y1, x2, y2;
			GetAllPointsOfTriangle(i, rectSize.y, x0, y0, x1, y1, x2, y2);  //Order is start from the up-left point, counter-clock wise

			//int tx = grid->location.size();//65
			//int ty = grid->location[0].size();//129

			glm::vec3 p0 = grid->location[x0][y0];
			glm::vec3 p1 = grid->location[x1][y1];
			glm::vec3 p2 = grid->location[x2][y2];

			triangles[i].AssignPosition(p0, p1, p2);

			glm::vec3 line10 = p1 - p0;
			glm::vec3 line21 = p2 - p1;
			glm::vec3 line02 = p0 - p2;
			
			float e10 = getLength(line10);//for e10, this is x
			float e21 = getLength(line21);
			float e02 = getLength(line02);

			float e10_oppositeAngle = acos((e21*e21 + e02*e02 - e10*e10) / (2 * e21*e02));
			float e21_oppositeAngle = acos((e10*e10 + e02*e02 - e21*e21) / (2 * e10*e02));
			float e02_oppositeAngle = acos((e21*e21 + e10*e10 - e02*e02) / (2 * e21*e10));


			halfedges[i * 3 + 0].start.x = x0;
			halfedges[i * 3 + 0].start.y = y0;
			halfedges[i * 3 + 0].end.x = x1;
			halfedges[i * 3 + 0].end.y = y1;
			halfedges[i * 3 + 0].cotan = 1.0 / glm::tan(e10_oppositeAngle);

			halfedges[i * 3 + 1].start.x = x1;
			halfedges[i * 3 + 1].start.y = y1;
			halfedges[i * 3 + 1].end.x = x2;
			halfedges[i * 3 + 1].end.y = y2;
			halfedges[i * 3 + 1].cotan = 1.0 / glm::tan(e21_oppositeAngle);

			halfedges[i * 3 + 2].start.x = y2;
			halfedges[i * 3 + 2].start.y = y2;
			halfedges[i * 3 + 2].end.x = x0;
			halfedges[i * 3 + 2].end.y = y0;
			halfedges[i * 3 + 2].cotan = 1.0 / glm::tan(e02_oppositeAngle);
		}

		///////////////////////////////////////////////////
#ifdef DEBUG
		for (int i = 0; i < halfedges.size(); i++) {
			printer_cotan.PrintToFile("Halfedge (");
			printer_cotan.PrintToFile(halfedges[i].start.x, 0);
			printer_cotan.PrintToFile(",");
			printer_cotan.PrintToFile(halfedges[i].start.y, 0);
			printer_cotan.PrintToFile(")->(");
			printer_cotan.PrintToFile(halfedges[i].end.x, 0);
			printer_cotan.PrintToFile(",");
			printer_cotan.PrintToFile(halfedges[i].end.y, 0);
			printer_cotan.PrintToFile("):        Cotan");
			printer_cotan.PrintToFile_Space(halfedges[i].cotan);
			printer_cotan.PrintToFile_Enter("");
		}
		printer_cotan.CloseFile();
#endif
		
	}

	void InitialCotan(CGrid *grid) {//bugged
#ifdef DEBUG
		CPrint printer_cotan;
		string filename_cotan = "cotan.txt";
		printer_cotan.CreateFile(OUTPUT_DATA + filename_cotan);
#endif
		///////////////////////////////////////////////////////////////

		for (int i = 0; i < triangles.size(); i++) {
			float x0, y0, x1, y1, x2, y2;
			GetAllPointsOfTriangle(i, rectSize.y, x0, y0, x1, y1, x2, y2);  //Order is start from the up-left point, counter-clock wise

			glm::vec3 p0 = grid->location[x0][y0];
			glm::vec3 p1 = grid->location[x1][y1];
			glm::vec3 p2 = grid->location[x2][y2];

			XVec2 e01(p1.x - p0.x, p1.y - p0.y);
			XVec2 e10(p0.x - p1.x, p0.y - p1.y);
			XVec2 e12(p2.x - p1.x, p2.y - p1.y);
			XVec2 e21(p1.x - p2.x, p1.y - p2.y);
			XVec2 e20(p0.x - p2.x, p0.y - p2.y);
			XVec2 e02(p2.x - p0.x, p2.y - p0.y);

			float a = e01.GetLength();
			float b = e12.GetLength();
			float c = e20.GetLength();
			float s = (a + b + c) / 2.0;
			float area = sqrt(s * (s - a)*(s - b)*(s - c));

			triangles[i].area = area;

			halfedges[i * 3 + 0].start.x = x0;
			halfedges[i * 3 + 0].start.y = y0;
			halfedges[i * 3 + 0].end.x = x1;
			halfedges[i * 3 + 0].end.y = y1;
			halfedges[i * 3 + 0].cotan = e21.Dot(e20) / 2.0f / area;

			halfedges[i * 3 + 1].start.x = x1;
			halfedges[i * 3 + 1].start.y = y1;
			halfedges[i * 3 + 1].end.x = x2;
			halfedges[i * 3 + 1].end.y = y2;
			halfedges[i * 3 + 1].cotan = e01.Dot(e02) / 2.0f / area;

			halfedges[i * 3 + 2].start.x = y2;
			halfedges[i * 3 + 2].start.y = y2;
			halfedges[i * 3 + 2].end.x = x0;
			halfedges[i * 3 + 2].end.y = y0;
			halfedges[i * 3 + 2].cotan = e10.Dot(e12) / 2.0f / area;
		}

		////////////////////////////////////////////////////////////
#ifdef DEBUG
		for (int i = 0; i < halfedges.size(); i++) {
			printer_cotan.PrintToFile("Halfedge (");
			printer_cotan.PrintToFile(halfedges[i].start.x, 0);
			printer_cotan.PrintToFile(",");
			printer_cotan.PrintToFile(halfedges[i].start.y, 0);
			printer_cotan.PrintToFile(")->(");
			printer_cotan.PrintToFile(halfedges[i].end.x, 0);
			printer_cotan.PrintToFile(",");
			printer_cotan.PrintToFile(halfedges[i].end.y, 0);
			printer_cotan.PrintToFile("):        Cotan");
			printer_cotan.PrintToFile_Space(halfedges[i].cotan);
			printer_cotan.PrintToFile_Enter("");
		}
		printer_cotan.CloseFile();
#endif
		
	}


	float getCotan(int n, int i, bool useNew) {
		if (n < 0 || n >= triangles.size()) return 0;
		int hIdx = n * 3 + i;
		if (useNew) return halfedges[hIdx].cotan_new;
		return halfedges[hIdx].cotan;
	}

	void AssembleLaplace(CSparseMatrix *sparseMatrix, bool useNewCotan, float h, string filename) {
#ifdef DEBUG
		CPrint printer;
		//string filename_LId = "L.txt";
		printer.CreateFile(OUTPUT_DATA + filename);
#endif

		int k = (LCSSize.y - 1) * 2; //how many triangles in a row

		for (int i = 0; i < LCSSize.x; i++) { 
			for (int j = 0; j < LCSSize.y; j++) {

				int x, y;
				int n0 = (i - 1)  * k + 2 * (j - 1); //if the neighbor triangle exists, this is it's triangle index
				int n1 = n0 + 1;
				int n2 = n1 + 1;
				int n3 = n0 + k + 1;
				int n4 = n3 + 1;
				int n5 = n4 + 1;
				float cotan_sum = 0;

				x = i - 1; y = j - 1; //neighbor point
				if (x >= 0 && x < LCSSize.x && y >= 0 && y < LCSSize.y) {
					float cotan = (getCotan(n0, 2, useNewCotan) + getCotan(n1, 0, useNewCotan)) / 2.0f * h;
					sparseMatrix->add1Value(GetIndexOfPoint(i, j, LCSSize.y), GetIndexOfPoint(x, y, LCSSize.y), -cotan); //return i * w + j;
					cotan_sum += cotan;
				}

				x = i - 1; y = j;
				if (x >= 0 && x < LCSSize.x && y >= 0 && y < LCSSize.y) {
					float cotan = (getCotan(n1, 1, useNewCotan) + getCotan(n2, 0, useNewCotan)) / 2.0f * h;
					sparseMatrix->add1Value(GetIndexOfPoint(i, j, LCSSize.y), GetIndexOfPoint(x, y, LCSSize.y), -cotan);
					cotan_sum += cotan;
				}

				x = i; y = j + 1;
				if (x >= 0 && x < LCSSize.x && y >= 0 && y < LCSSize.y) {
					float cotan = (getCotan(n2, 1, useNewCotan) + getCotan(n5, 2, useNewCotan)) / 2.0f * h;
					sparseMatrix->add1Value(GetIndexOfPoint(i, j, LCSSize.y), GetIndexOfPoint(x, y, LCSSize.y), -cotan);
					cotan_sum += cotan;
				}
				
				x = i + 1; y = j + 1;
				if (x >= 0 && x < LCSSize.x && y >= 0 && y < LCSSize.y) {
					float cotan = (getCotan(n5, 0, useNewCotan) + getCotan(n4, 2, useNewCotan)) / 2.0f * h;
					sparseMatrix->add1Value(GetIndexOfPoint(i, j, LCSSize.y), GetIndexOfPoint(x, y, LCSSize.y), -cotan);
					cotan_sum += cotan;
				}
				
				x = i + 1; y = j;
				if (x >= 0 && x < LCSSize.x && y >= 0 && y < LCSSize.y) {
					//float cotan = getCotan(n4, 0, n3, 1, useNewCotan) * h;
					float cotan = (getCotan(n4, 0, useNewCotan) + getCotan(n3, 1, useNewCotan)) / 2.0f * h;
					sparseMatrix->add1Value(GetIndexOfPoint(i, j, LCSSize.y), GetIndexOfPoint(x, y, LCSSize.y), -cotan);
					cotan_sum += cotan;
				}

				x = i; y = j - 1;
				if (x >= 0 && x < LCSSize.x && y >= 0 && y < LCSSize.y) {
					//float cotan = getCotan(n3, 2, n0, 1, useNewCotan) * h;
					float cotan = (getCotan(n3, 2, useNewCotan) + getCotan(n0, 1, useNewCotan)) / 2.0f * h;
					sparseMatrix->add1Value(GetIndexOfPoint(i, j, LCSSize.y), GetIndexOfPoint(x, y, LCSSize.y), -cotan);
					cotan_sum += cotan;
				}

				sparseMatrix->add1Value(GetIndexOfPoint(i, j, LCSSize.y), GetIndexOfPoint(i, j, LCSSize.y), cotan_sum);

				//float dualArea = 
			}
		}

#ifdef DEBUG
		CMatrixElement *theElem;
		for (int i = 0; i < sparseMatrix->numRows; i++)
			for (theElem = sparseMatrix->rowList[i]; theElem != NULL; theElem = theElem->rowNext) {
				printer.PrintToFile("i = ");
				printer.PrintToFile(theElem->i, 0);
				printer.PrintToFile(", j = ");
				printer.PrintToFile(theElem->j, 0);
				printer.PrintToFile(": ");
				printer.PrintToFile_Enter(theElem->value);
			}
		printer.CloseFile();
#endif
		
	}


	float MatlabEigenSolver(CMatlabEng &matlabEng, CSparseMatrix *A, CSparseMatrix *B, CSparseMatrix *L, int stepNum) { //, CSparseMatrix *L_alpha, CSparseMatrix *L_Id
		matlabEng.PutVar_sparse("A", A);
		matlabEng.PutVar_sparse("B", B);
		matlabEng.PutVar_sparse("L", L);

		std::string command;
		int sampleSize = 1;
		std::ostringstream E; E << sampleSize;
		std::ostringstream N; N << LCSSize.x * LCSSize.y;


		cout << "-Smoothing S..." << endl;
		double *array_S = new double[LCSSize.x * LCSSize.y];
		for (int i = 0; i < LCSSize.x; i++)
			for (int j = 0; j < LCSSize.y; j++)
				array_S[i * (int)LCSSize.y + j] = S[i][j];
		matlabEng.PutVar_d("S", array_S, LCSSize.x * LCSSize.y);

		//Explicit: M * S^(t+1) = (M + h * L) * S^t
		//command = "L_smooth = 0.01 * L + speye(" + N.str() + ")"; matlabEng.EvalString(command.c_str());
		//command = "R = L_smooth * S"; matlabEng.EvalString(command.c_str());

		//smooth_h = 0.025;
		//smooth_n = 1;

		std::ostringstream SH; SH << smooth_h;
		
		//Implicit: (M - h * L) * S^(t+1) = M * S^t
		command = "L_smooth = speye(" + N.str() + ") - " + SH.str() + " * L"; matlabEng.EvalString(command.c_str());
		//command = "R = linsolve(L_smooth, S)"; matlabEng.EvalString(command.c_str()); //not work 
		for (int i = 0; i < smooth_n; i++) {
			command = "S = L_smooth\\S"; matlabEng.EvalString(command.c_str());
		}
		double* array_tmp = new double[LCSSize.x * LCSSize.y];
		if (smooth_n > 0) {
			matlabEng.GetVar_d("S", array_tmp);
			for (int i = 0; i < LCSSize.x * LCSSize.y; i++) {
				//S[i / LCSSize.y][i % (int)LCSSize.y] = array_tmp[i];
				//if (array_tmp[i] > maxS) maxS = array_tmp[i];
				//if (array_tmp[i] < minS) minS = array_tmp[i];
			}
		}

		cout << "-Running eigensolver..." << endl;
		//command = "[V,D] = eigs(A, B+speye(" + N.str() + ")*1e-6, " + E.str() + ",\'lm\')"; matlabEng.EvalString(command.c_str());
		command = "[V,D] = eigs(A, B, " + E.str() + ",\'lm\')"; matlabEng.EvalString(command.c_str());

		//错误使用/ 内存不足。请键入HELP MEMORY查看选项
		//出错 eigs/checkInputs (line 747) B = B / scaleB;
		//出错 eigs (line 96) checkInputs(varargin{:});

		matlabEng.GetVar_d("V", array_tmp);
		for (int i = 0; i < LCSSize.x * LCSSize.y; i++){
			//testScaler[i / LCSSize.y][i % (int)LCSSize.y] = tmp[i];//F
			F[i / LCSSize.y][i % (int)LCSSize.y] = array_tmp[i];
			if (array_tmp[i] > maxF) maxF = array_tmp[i];
			if (array_tmp[i] < minF) minF = array_tmp[i];
		}

		//command = "R = B*V"; matlabEng.EvalString(command.c_str());
		//matlabEng.GetVar_d("R", tmp); 
		//for (int i = 0; i < LCSSize.x * LCSSize.y; i++) {
		//	testScaler[i / LCSSize.y][i % (int)LCSSize.y] = tmp[i];
		//}

		int count = 0;
		float sum = 0;
		for (int i = 0; i < LCSSize.x; i++)
			for (int j = 0; j < LCSSize.y; j++) {
				if (S[i][j] >= 5) {
					count++;
					sum += F[i][j];
				}
				
			}
		float aver = sum / count;

		marchingTriangle.MatlabContour(matlabEng, F, aver, stepNum, currentDirectory);

		//std::ostringstream WIDTH; WIDTH << LCSSize.y;
		//std::ostringstream HEIGHT; HEIGHT << LCSSize.x;
		//command = "F = reshape(V, " + WIDTH.str() + ", " + HEIGHT.str() + ");"; matlabEng.EvalString(command.c_str());
		//command = "F = F';"; matlabEng.EvalString(command.c_str());
		//std::ostringstream AVER; AVER << aver;
		//command = "[C,h] = contour(F, [" + AVER.str() + " " + AVER.str() + "])"; matlabEng.EvalString(command.c_str());

		//float ratio = LCSSize.y / LCSSize.x;
		//std::ostringstream RATIO; RATIO << ratio;
		//command = "pbaspect(["+RATIO.str()+" 1 1])"; matlabEng.EvalString(command.c_str());

		//std::ostringstream STEPNUM; STEPNUM << stepNum;
		//command = "saveas(h, '" + currentDirectory + "\\" + OUTPUT_DATA + "imgs\\LCS\\LCS" + STEPNUM.str() + ".png')";
		//matlabEng.EvalString(command.c_str());

		delete array_S;
		delete array_tmp;

		return aver;
	}

	~CLCS() {
		if (L_Alpha) {
			L_Alpha->CleanUp();
			delete L_Alpha;
		}
		if (L_B) {
			L_B->CleanUp();
			delete L_B;
		}

		if (L_Id) {
			L_Id->CleanUp();
			delete L_Id;
		}
	}

};