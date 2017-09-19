#pragma once
#include "cgrid.h"
#include "clcs.h"
#include "cdisplaycontroller.h"
#include "..\cmathhelper.h"
//#include "cmarchingtriangle.h"

//https://www.mathworks.com/help/matlab/calling-matlab-engine-from-c-c-and-fortran-programs.html
#include "..\matlab\MatlabEng.h"

#include "..\mathutilities\SparseMatrix.h"

#include "..\print\Print.h"

class CFieldSystem {
public:
	CGrid *grid;
	vector<CNodes> nodes;
	CLCS *lcs;
	CDisplayController *displayController;
	int m_width;
	int m_height;
	
	glm::vec2 centerIdx1;
	glm::vec2 centerIdx2;
	glm::vec3 centerLocation;

	float t; //current time
	
	//float T; //backward/forward track time interval
	//float dt; //backward/forward track time

	float offline_t; //current offline time
	float offline_dt;

	CMatlabEng matlabEng;

	bool enableTestFTLE;
	bool enableTestEigenvectorOfHessian;
	bool enableTestS;

	float epsilon_s;
	float thresh_s;
	float epsilon_tensor;
	float threshold_FTLE_s;

	int m_maxSteps;

	//int integral_step;
	//CMarchingTriangle marchingTriangle;
	

	CFieldSystem(int width, int height, int maxSteps) {
		enableTestFTLE = 0;
		enableTestEigenvectorOfHessian = 0;
		enableTestS = 0;

		t = 0;

		//T = 15;
		//dt = 0.05;

		offline_t = 0;
		offline_dt = 0.1;

		epsilon_s = 0;
		thresh_s = 0;
		epsilon_tensor = 0;
		threshold_FTLE_s = 0;
		
		m_width = width;
		m_height = height;

		m_maxSteps = maxSteps;

		grid = new CGrid(width, height);

		//for(int i = 0; i < maxSteps; i++)
			//nodes.push_back(CNodes(width, height, i));

		lcs = new CLCS(width, height);
		displayController = new CDisplayController(width, height);

		cout << "Opening Matlab..." << endl;
		matlabEng.Open(NULL);

		//lcs->InitialTriangles(grid);
		//lcs->InitialCotan(grid);
		lcs->InitialCotan_fromAngle(grid);
	
	}

	void Update(float secondsElapsed, int stepNum, int integral_step) {
		lcs->ResetMinMax();

		cout << "Update: Current time: " << t << endl;
		cout << "Update: Current offline time: " << offline_t << endl;
	
		cout << "Update: Computing Phi..." << endl;
		ComputePhi(stepNum, integral_step);

		cout << "Update: Computing FTLE..." << endl;
		ComputeFTLE();
		
		cout << "Update: Computing S..." << endl;
		ComputeS();

		cout << "Update: Computing triange tensors..." << endl;
		ComputeTensorForTriangle();

		cout << "Update: Computing cotan for halfedges..." << endl;
		//UpdateCotanForHalfedge();
		UpdateCotanForHalfedge_fromAngle();


		cout << "Update: Assemble L_Alpha..." << endl;
		if (lcs->L_Alpha) {
			lcs->L_Alpha->CleanUp();
			delete lcs->L_Alpha;
		}
		lcs->L_Alpha = new CSparseMatrix(lcs->LCSSize.x * lcs->LCSSize.y);
		lcs->AssembleLaplace(lcs->L_Alpha, true, 1.0, "L_Alpha.txt");

		matlabEng.PutVar_sparse("L_alpha", lcs->L_Alpha);//for test

		cout << "-Construct B matrix..." << endl;
		if (lcs->L_B) {
			lcs->L_B->CleanUp();
			delete lcs->L_B;
		}
		lcs->L_B = new CSparseMatrix(lcs->LCSSize.x * lcs->LCSSize.y);
		lcs->AssembleLaplace(lcs->L_B, true, 1.0, "L_B.txt");//?
		//B = L'L + wD
		lcs->L_B->multTransMatMat_yz();
		//LB->multTransMatMat();

		matlabEng.PutVar_sparse("DL_alpha", lcs->L_B);//for test

		int scount = 0;
		for (int i = 0; i < lcs->LCSSize.x; i++) {
			for (int j = 0; j < lcs->LCSSize.y; j++) {
				int idx = i * lcs->LCSSize.y + j;
				if (lcs->S[i][j] >= 5) {
					lcs->L_B->add1Value(idx, idx, lcs->weight_b);
					scount++;
				}
				lcs->L_B->add1Value(idx, idx, 0.00001);
			}
		}

		cout << "Update: Assemble L_Id..." << endl;
		if (lcs->L_Id) {
			lcs->L_Id->CleanUp();
			delete lcs->L_Id;
		}
		lcs->L_Id = new CSparseMatrix(lcs->LCSSize.x * lcs->LCSSize.y);
		lcs->AssembleLaplace(lcs->L_Id, false, 1.0, "L_Id.txt");
		
		//Make B matrix, B = L_Id + weight * D, while Dij = 1 if i=j in S
		//cout << "Update: Assemble L_Id + w * D..." << endl;
		//if (lcs->L_IdwD) delete lcs->L_IdwD;
		//lcs->L_IdwD = new CSparseMatrix(lcs->LCSSize.x * lcs->LCSSize.y);
		//lcs->AssembleLaplace(lcs->L_IdwD, false, 1.0, "L_IdwD.txt");
		//float weight = 0.0000001;
		//for (int i = 0; i < lcs->LCSSize.x; i++)
		//	for (int j = 0; j < lcs->LCSSize.y; j++)
		//		if (lcs->S[i][j] > 1000)
		//			lcs->L_IdwD->add1Value(i *  lcs->LCSSize.y + j, i * lcs->LCSSize.y + j, weight);


		//smooth S using (M+h*L_Id)S
		//cout << "Update: Smoothing S..." << endl;
		//float h = 1;
		//float mass = 0.0f;
		//if (lcs->L_S) delete lcs->L_S;
		//lcs->L_S = new CSparseMatrix(lcs->LCSSize.x * lcs->LCSSize.y);
		//lcs->AssembleLaplace(lcs->L_S, false, h, "L_S.txt");

		//for (int i = 0; i < lcs->LCSSize.x; i++) 
		//	for (int j = 0; j < lcs->LCSSize.y; j++) 
		//		lcs->L_S->add1Value(i * lcs->LCSSize.y + j, i * lcs->LCSSize.y + j, mass);
		//lcs->SmoothSinMatlab(matlabEng);
		//lcs->TestMatlab(matlabEng);
	
		cout << "Update: Solve for L_Alpha * f = lambda_max * L_IdwD * f..." << endl;
		float isovalue = lcs->MatlabEigenSolver(matlabEng, lcs->L_Alpha, lcs->L_B, lcs->L_Id, stepNum);


		cout << "Update contour..." << endl;
		lcs->marchingTriangle.Contour(lcs->F, isovalue, lcs->C, lcs->triangles, "Countour.txt");


		cout << "Update: Summary" << endl;
		cout << "Min FTLE: " << lcs->minFTLE << ", Max FTLE: " << lcs->maxFTLE << endl;
		cout << "Min S: " << lcs->minS << ", Max S: " << lcs->maxS << endl;
		cout << "Min F: " << lcs->minF << ", Max F: " << lcs->maxF << endl;

		t += secondsElapsed;
		offline_t += offline_dt;
	}

	//void ParseVelocity(float px, float py, float *vx, float *vy, float delta) {
	//	//input: px: 0 ~ width-1; py: 0 ~ height-1
	//	//rk4 variables
	//	glm::vec2 p1, p2, p3, p4;
	//	glm::vec2 v1, v2, v3, v4;

	//	p1.x = px;
	//	p1.y = py;
	//	v1.x = float(-PI * sin(2 * p1.x * PI) * cos(p1.y * PI));
	//	v1.y = float(PI * cos(2 * p1.x * PI) * sin(p1.y * PI));

	//	p2.x = float(p1.x + delta / 2.0f * v1.x);
	//	p2.y = float(p1.y + delta / 2.0f * v1.y);
	//	//GridInterpolation(p2.X, p2.Y, &v2.X, &v2.Y);//interpolate the velocity
	//	v2.x = float(-PI * sin(2 * p2.x * PI) * cos(p2.y * PI));
	//	v2.y = float(PI * cos(2 * p2.x * PI) * sin(p2.y * PI));

	//	p3.x = float(p1.x + delta / 2.0f * v2.x);
	//	p3.y = float(p1.y + delta / 2.0f * v2.y);
	//	//GridInterpolation(p3.X, p3.Y, &v3.X, &v3.Y);//interpolate the velocity
	//	v3.x = float(-PI * sin(2 * p3.x * PI) * cos(p3.y * PI));
	//	v3.y = float(PI * cos(2 * p3.x * PI) * sin(p3.y * PI));

	//	p4.x = float(p1.x + delta * v3.x);
	//	p4.y = float(p1.y + delta * v3.y);
	//	//GridInterpolation(p4.X, p4.Y, &v4.X, &v4.Y);//interpolate the velocity
	//	v4.x = float(-PI * sin(2 * p4.x * PI) * cos(p4.y * PI));
	//	v4.y = float(PI * cos(2 * p4.x * PI) * sin(p4.y * PI));

	//	*vx = (v1.x + 2.0f * v2.x + 2.0f * v3.x + v4.x) / 6.0f;
	//	*vy = (v1.y + 2.0f * v2.y + 2.0f * v3.y + v4.y) / 6.0f;
	//}



	//void AdvectPosition(float &px, float &py, float t, float dt) {
	//	//http://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html#x1-1200511
	//	//float A = 0.1;
	//	//float dfdx;
	//	//float fx = f(px, t, dfdx);
	//	//float vx = (float)(-PI * A * sin(PI * fx) * cos(PI * py));
	//	//float vy = (float)(PI * A * cos(PI * fx) * sin(PI * py) * dfdx);

	//	glm::vec2 v = GetVelocity(px, py, t);

	//	px += v.x * dt;
	//	py += v.y * dt;

	//	//return glm::vec2(vx, vy);
	//}


	void ComputePhi(int stepNum, int integral_step) { //phi(x) = x + h * v(x)
		if (m_maxSteps == 0)
			nodes.push_back(CNodes(m_width, m_height, stepNum));
		
		for (int i = 0; i < lcs->gridSize.x; i++) {
			for (int j = 0; j < lcs->gridSize.y; j++) {
				float px = grid->location[i][j].x;
				float py = grid->location[i][j].y;

				//only for visualization
				if (m_maxSteps == 0) {
					glm::vec2 v = GetVelocity(px, py, offline_t);
					nodes[stepNum].SetVelocity(i, j, v.x, v.y, 0);
				}

				for (int k = 0; k < integral_step; k++) {//backward track or forward track for k steps
					glm::vec2 v;
					if (m_maxSteps == 0)	//use generated data
						v = GetVelocity(px, py, offline_t + k * lcs->integrat_dt);
						//AdvectPosition(px, py, offline_t + k * lcs->integrat_dt, lcs->integrat_dt);
					else //use data read from file
						v = GetVelocity_fromFile(px, py, k, stepNum, grid, nodes, i, j);

					if (lcs->integrat_T >= 0) {
						px += v.x * lcs->integrat_dt;
						py += v.y * lcs->integrat_dt;
					}
					else {
						px -= v.x * lcs->integrat_dt;
						py -= v.y * lcs->integrat_dt;
					}
				}
			


				lcs->phi[i][j].x = px;
				lcs->phi[i][j].y = py;
				float magnitude = sqrt(px*px + py*py);
				if (magnitude > lcs->maxPhi) lcs->maxPhi = magnitude;
				if (magnitude < lcs->minPhi) lcs->minPhi = magnitude;
			}
		}
	}

	void ComputeFTLE() {
		//float max = -9999;
		//float min = 9999;
		//float max_lambda = -9999;
		for (int i = 0; i < lcs->FTLESize.x; i++) {
			for (int j = 0; j < lcs->FTLESize.y; j++) {
				/***gradPhi = [ dphi/dx | dphi/dy]***/
				//float test1 = grid->dx;
				//float test2 = grid->dy;
				lcs->gradPhi[i][j] = glm::mat2x2(
					(lcs->phi[i][j + 1][0] - lcs->phi[i][j][0]) / grid->dx, (lcs->phi[i + 1][j][0] - lcs->phi[i][j][0]) / grid->dy,
					(lcs->phi[i][j + 1][1] - lcs->phi[i][j][1]) / grid->dx, (lcs->phi[i + 1][j][1] - lcs->phi[i][j][1]) / grid->dy);

				/***Delta = gradPhi * (gradPhi)^T***/
				glm::mat2x2 gradPhi_transpose(
					lcs->gradPhi[i][j][0][0], lcs->gradPhi[i][j][1][0],
					lcs->gradPhi[i][j][0][1], lcs->gradPhi[i][j][1][1]);
				glm::mat2x2 Delta = lcs->gradPhi[i][j] * gradPhi_transpose; 

				/***FTLE = ln(maxEigenvalue)^2 / h***/
				glm::vec2 lambda = GetEigenvalue(Delta);
				if (!Trivial(lambda[0], TRIVIAL_NUMBER))
					//lcs->FTLE[i][j] = log(sqrt(lambda[0])) / abs(secondsElapsed);
					lcs->FTLE[i][j] = log(sqrt(lambda[0])) / abs(lcs->integrat_T);
				else 
					lcs->FTLE[i][j] = 0;

				//float test_lambda_max = lambda[0];
				//float test_lambda_min = lambda[1];
				//float test_FTLE = lcs->FTLE[i][j];

				//int test = 1;
				//test++;

				//if (lambda[0] > max_lambda) max_lambda = lambda[0];
				//if (lcs->FTLE[i][j] > max) max = lcs->FTLE[i][j];
				//if (lcs->FTLE[i][j] < min) min = lcs->FTLE[i][j];

				if (lcs->FTLE[i][j] > lcs->maxFTLE) lcs->maxFTLE = lcs->FTLE[i][j];
				if (lcs->FTLE[i][j] < lcs->minFTLE) lcs->minFTLE = lcs->FTLE[i][j];
				
			}
		}

		//cout << "Max Lambda of grad(phi): " << max_lambda << endl;
		//cout << "Max FTLE: " << max << endl;
		//cout << "Min FTLE: " << min << endl;

		if (enableTestFTLE) {
			cout << "!!!Enable test FTLE..." << endl;
			lcs->FTLE = lcs->testFTLE;
		}
	}

	void ComputeS() {
		//The two creterias: http://shaddenlab.berkeley.edu/uploads/LCS-tutorial/LCSdef.html
#ifdef DEBUG
		CPrint printer_s;
		string filename_s = "S.txt";
		printer_s.CreateFile(OUTPUT_DATA + filename_s);
#endif
		/***calculate the gradient of the exponent***/
		//float max = -9999;
		//float min = 9999;
		//cout << "Creating gradient of FTLE..." << endl;


		for (int i = 0; i < lcs->LCSSize.x; i++) {
			for (int j = 0; j < lcs->LCSSize.y; j++) {
				lcs->gradFTLE[i][j].x = (lcs->FTLE[i + 1][j + 2] - lcs->FTLE[i + 1][j]) / (2 * grid->dx);
				lcs->gradFTLE[i][j].y = (lcs->FTLE[i + 2][j + 1] - lcs->FTLE[i][j + 1]) / (2 * grid->dy);
			}
		}

		if (enableTestEigenvectorOfHessian)
			cout << "!!!Enable test eigenvector of hessian..." << endl;

		//calculate the hessian matrix of exponent and min and max eigenvectors of the hessian
		for (int i = 0; i < lcs->LCSSize.x; i++) {
			for (int j = 0; j < lcs->LCSSize.y; j++) {

				//f00 correspond to (i+1,j+1)
				float f0_1 = lcs->FTLE[i + 1][j];
				float f_10 = lcs->FTLE[i][j + 1];
				float f_1_1 = lcs->FTLE[i][j];

				float f00 = lcs->FTLE[i + 1][j + 1];

				float f01 = lcs->FTLE[i + 1][j + 2];
				float f10 = lcs->FTLE[i + 2][j + 1];
				float f11 = lcs->FTLE[i + 2][j + 2];

				float f1_1 = lcs->FTLE[i + 2][j];
				float f_11 = lcs->FTLE[i][j + 2];

				float value1 = (f01 + f0_1 - 2 * f00) / grid->dx / grid->dx;
				float value2 = (f11 - f1_1 - f_11 + f_1_1) / 4 / grid->dx / grid->dy;
				float value4 = (f10 + f_10 - 2 * f00) / grid->dy / grid->dy;

				lcs->hessian[i][j][0][0] = value1; lcs->hessian[i][j][0][1] = value2;
				lcs->hessian[i][j][1][0] = value2; lcs->hessian[i][j][1][1] = value4;

				GetEigenInfo(lcs->hessian[i][j], lcs->eigenvectorOfHessian[i][j], lcs->eigenvaluesOFHessian[i][j]);
				glm::vec2 lambda = lcs->eigenvaluesOFHessian[i][j];
				if (enableTestEigenvectorOfHessian)lcs->eigenvectorOfHessian[i][j] = lcs->testEigenvectorOfHessian[i][j];//for test purpose
				glm::vec2 e0(lcs->eigenvectorOfHessian[i][j][0][0], lcs->eigenvectorOfHessian[i][j][1][0]);
				glm::vec2 e1(lcs->eigenvectorOfHessian[i][j][0][1], lcs->eigenvectorOfHessian[i][j][1][1]);

				float magnitude = sqrt(e0.x*e0.x + e0.y*e0.y);
				if (magnitude > lcs->maxEigenvectorHessian1) lcs->maxEigenvectorHessian1 = magnitude;
				if (magnitude < lcs->minEigenvectorHessian1) lcs->minEigenvectorHessian1 = magnitude;
				magnitude = sqrt(e1.x*e1.x + e1.y*e1.y);
				if (magnitude > lcs->maxEigenvectorHessian2) lcs->maxEigenvectorHessian2 = magnitude;
				if (magnitude < lcs->minEigenvectorHessian2) lcs->minEigenvectorHessian2 = magnitude;

				float norm = sqrt(lcs->gradFTLE[i][j].x * lcs->gradFTLE[i][j].x + lcs->gradFTLE[i][j].y * lcs->gradFTLE[i][j].y);
				glm::vec2 normed_gradFTLE(lcs->gradFTLE[i][j].x / norm, lcs->gradFTLE[i][j].y / norm);

				float t = dotProduct(normed_gradFTLE, e1);
				lcs->S[i][j] = 1.0f / (t*t + epsilon_s);

				if (lambda[1] > thresh_s) lcs->S[i][j] = 1;
				//if (lcs->FTLE[i + 1][j + 1] < (9.0f / 30.0f)) 
				//if (lcs->FTLE[i + 1][j + 1] < (20.0f / 30.0f))
				if (lcs->FTLE[i + 1][j + 1] < threshold_FTLE_s)
					lcs->S[i][j] = 1;
				else 
					lcs->S[i][j] = 5;

				if (lcs->S[i][j] > lcs->maxS) lcs->maxS = lcs->S[i][j];
				if (lcs->S[i][j] < lcs->minS) lcs->minS = lcs->S[i][j];
			}
		}

		if (enableTestS) {
			cout << "!!!Enable test S..." << endl;
			float centerW = lcs->FTLESize.y / 2;
			float centerH = lcs->FTLESize.x / 2;
			for (int i = 0; i < lcs->LCSSize.x; i++) {
				for (int j = 0; j < lcs->LCSSize.y; j++) {
					float w = j - centerW;
					float h = i - centerH;
					float radius = sqrt(w*w + h*h);
					//if (radius < 0.0001) radius = 0.0001;
					if (radius > 24 && radius < 26)
						lcs->S[i][j] = 5;
					else
						lcs->S[i][j] = 1;

					if (lcs->S[i][j] > lcs->maxS) lcs->maxS = lcs->S[i][j];
					if (lcs->S[i][j] < lcs->minS) lcs->minS = lcs->S[i][j];
				}
			}
			enableTestS = false;
		}

#ifdef DEBUG
		for (int i = 0; i < lcs->S.size(); i++) {
			for (int j = 0; j < lcs->S[0].size(); j++) {
				printer_s.PrintToFile("i = ");
				printer_s.PrintToFile(i, 0);
				printer_s.PrintToFile(", j = ");
				printer_s.PrintToFile(j, 0);
				printer_s.PrintToFile(": ");
				printer_s.PrintToFile_Enter(lcs->S[i][j]);
			}
		}

		printer_s.CloseFile();
#endif
	}

	void UpdateCotanForHalfedge_fromAngle() {
#ifdef DEBUG
		CPrint printer_cotan;
		string filename_cotan = "cotan_new_from_angle.txt";
		printer_cotan.CreateFile(OUTPUT_DATA + filename_cotan);
#endif
		/////////////////////////////////////////
		for (int i = 0; i < lcs->triangles.size(); i++) {
			float x0, y0, x1, y1, x2, y2;
			GetAllPointsOfTriangle(i, lcs->rectSize.y, x0, y0, x1, y1, x2, y2);  //Order is start from the up-left point, counter-clock wise

			glm::vec3 p0 = grid->location[x0][y0];
			glm::vec3 p1 = grid->location[x1][y1];
			glm::vec3 p2 = grid->location[x2][y2];

			glm::vec3 line10 = p1 - p0;
			glm::vec3 line21 = p2 - p1;
			glm::vec3 line02 = p0 - p2;

			//---
			XMat22 t(lcs->triangles[i].tensor_triangle[0][0], lcs->triangles[i].tensor_triangle[0][1],
				lcs->triangles[i].tensor_triangle[1][0], lcs->triangles[i].tensor_triangle[1][1]);
			//t = XMat22(1, 0, 0, 1);//test
			XVec2 xline10(line10.x, line10.y);
			XVec2 xline21(line21.x, line21.y);
			XVec2 xline02(line02.x, line02.y);

			float e10 = sqrt(xline10.Dot(t * xline10));//for e10, this is x
			float e21 = sqrt(xline21.Dot(t * xline21));
			float e02 = sqrt(xline02.Dot(t * xline02));
			//---

			float e10_oppositeAngle = acos((e21*e21 + e02*e02 - e10*e10) / (2 * e21*e02));
			float e21_oppositeAngle = acos((e10*e10 + e02*e02 - e21*e21) / (2 * e10*e02));
			float e02_oppositeAngle = acos((e21*e21 + e10*e10 - e02*e02) / (2 * e21*e10));


			float epsilon = 0.000001;
			e10_oppositeAngle = glm::max(e10_oppositeAngle, epsilon);
			e21_oppositeAngle = glm::max(e21_oppositeAngle, epsilon);
			e02_oppositeAngle = glm::max(e02_oppositeAngle, epsilon);

			//halfedges[i * 3 + 0].start.x = x0;
			//halfedges[i * 3 + 0].start.y = y0;
			//halfedges[i * 3 + 0].end.x = x1;
			//halfedges[i * 3 + 0].end.y = y1;
			
			lcs->halfedges[i * 3 + 0].cotan_new = 1.0 / glm::tan(e10_oppositeAngle);

			//halfedges[i * 3 + 1].start.x = x1;
			//halfedges[i * 3 + 1].start.y = y1;
			//halfedges[i * 3 + 1].end.x = x2;
			//halfedges[i * 3 + 1].end.y = y2;
			
			lcs->halfedges[i * 3 + 1].cotan_new = 1.0 / glm::tan(e21_oppositeAngle);

			//halfedges[i * 3 + 2].start.x = y2;
			//halfedges[i * 3 + 2].start.y = y2;
			//halfedges[i * 3 + 2].end.x = x0;
			//halfedges[i * 3 + 2].end.y = y0;
			lcs->halfedges[i * 3 + 2].cotan_new = 1.0 / glm::tan(e02_oppositeAngle);


			if (lcs->halfedges[i*3+1].start.x == 61 &&
				lcs->halfedges[i*3+1].start.y == 6 &&
				lcs->halfedges[i*3+1].end.x == 61 &&
				lcs->halfedges[i*3+1].end.y == 7) {
				cout << lcs->halfedges[i*3+1].cotan_new << endl;
			}

			if (lcs->halfedges[i * 3 + 2].start.x == 61 &&
				lcs->halfedges[i * 3 + 2].start.y == 6 &&
				lcs->halfedges[i * 3 + 2].end.x == 61 &&
				lcs->halfedges[i * 3 + 2].end.y == 7) {
				cout << lcs->halfedges[i * 3 + 2].cotan_new << endl;
			}
		}

		////////////////////////////////////////////
#ifdef DEBUG
		for (int i = 0; i < lcs->halfedges.size(); i++) {
			printer_cotan.PrintToFile("Halfedge (");
			printer_cotan.PrintToFile(lcs->halfedges[i].start.x, 0);
			printer_cotan.PrintToFile(",");
			printer_cotan.PrintToFile(lcs->halfedges[i].start.y, 0);
			printer_cotan.PrintToFile(")->(");
			printer_cotan.PrintToFile(lcs->halfedges[i].end.x, 0);
			printer_cotan.PrintToFile(",");
			printer_cotan.PrintToFile(lcs->halfedges[i].end.y, 0);
			printer_cotan.PrintToFile("):        Cotan");
			printer_cotan.PrintToFile_Space(lcs->halfedges[i].cotan_new);
			printer_cotan.PrintToFile_Enter("");
		}
		printer_cotan.CloseFile();
#endif
		
	}

	void UpdateCotanForHalfedge() {//bugged
#ifdef DEBUG
		CPrint printer_cotan;
		string filename_cotan = "cotan_new.txt";
		printer_cotan.CreateFile(OUTPUT_DATA + filename_cotan);
#endif
		for (int i = 0; i < lcs->triangles.size(); i++) {
			float x0, y0, x1, y1, x2, y2;
			GetAllPointsOfTriangle(i, lcs->rectSize.y, x0, y0, x1, y1, x2, y2);  //Order is start from the up-left point, counter-clock wise

			glm::vec3 p0 = grid->location[x0][y0];
			glm::vec3 p1 = grid->location[x1][y1];
			glm::vec3 p2 = grid->location[x2][y2];

			XVec2 e01(p1.x - p0.x, p1.y - p0.y);
			XVec2 e10(p0.x - p1.x, p0.y - p1.y);
			XVec2 e12(p2.x - p1.x, p2.y - p1.y);
			XVec2 e21(p1.x - p2.x, p1.y - p2.y);
			XVec2 e20(p0.x - p2.x, p0.y - p2.y);
			XVec2 e02(p2.x - p0.x, p2.y - p0.y);

			//XMat22 t(1, 0, 0, 1);//test

			XMat22 t(lcs->triangles[i].tensor_triangle[0][0], lcs->triangles[i].tensor_triangle[0][1], 
					 lcs->triangles[i].tensor_triangle[1][0], lcs->triangles[i].tensor_triangle[1][1]);

			lcs->halfedges[i * 3 + 0].cotan_new = e21.Dot(t * e20) / 2.0f / lcs->triangles[i].area;
			lcs->halfedges[i * 3 + 1].cotan_new = e01.Dot(t * e02) / 2.0f / lcs->triangles[i].area;
			lcs->halfedges[i * 3 + 2].cotan_new = e10.Dot(t * e12) / 2.0f / lcs->triangles[i].area;
		}

#ifdef DEBUG
		for (int i = 0; i < lcs->halfedges.size(); i++) {
			printer_cotan.PrintToFile("Halfedge (");
			printer_cotan.PrintToFile(lcs->halfedges[i].start.x, 0);
			printer_cotan.PrintToFile(",");
			printer_cotan.PrintToFile(lcs->halfedges[i].start.y, 0);
			printer_cotan.PrintToFile(")->(");
			printer_cotan.PrintToFile(lcs->halfedges[i].end.x, 0);
			printer_cotan.PrintToFile(",");
			printer_cotan.PrintToFile(lcs->halfedges[i].end.y, 0);
			printer_cotan.PrintToFile("):        Cotan");
			printer_cotan.PrintToFile_Space(lcs->halfedges[i].cotan);
			printer_cotan.PrintToFile_Enter("");
		}
		printer_cotan.CloseFile();
#endif
		
	}

	void ComputeTensorForTriangle() {
		//float epsilon = 0.00000001f;

		//Step1 compute tensor_vertex = S*ee^T+episilon*I in grid points
		for (int i = 0; i < lcs->LCSSize.x; i++) {
			for (int j = 0; j < lcs->LCSSize.y; j++) {
				XVec2 e(lcs->eigenvectorOfHessian[i][j][0][0], lcs->eigenvectorOfHessian[i][j][1][0]);
				XMat22 eeT(e(0) * e(0), e(0) * e(1), e(0) * e(1), e(1) * e(1));
				XMat22 I = XMat22(1, 0, 0, 1);
				XMat22 t = lcs->S[i][j] * lcs->S[i][j] * eeT + epsilon_tensor * I;
				lcs->tensor_vertex[i][j][0][0] = t(0, 0); lcs->tensor_vertex[i][j][0][1] = t(0, 1);
				lcs->tensor_vertex[i][j][1][0] = t(1, 0); lcs->tensor_vertex[i][j][1][1] = t(1, 1);

				GetEigenInfo(lcs->tensor_vertex[i][j], lcs->eigenvectorOfTensor[i][j], lcs->eigenvaluesOfTensor[i][j]);

				glm::vec2 e0(lcs->eigenvectorOfTensor[i][j][0][0], lcs->eigenvectorOfTensor[i][j][1][0]);
				glm::vec2 e1(lcs->eigenvectorOfTensor[i][j][0][1], lcs->eigenvectorOfTensor[i][j][1][1]);

				float magnitude = sqrt(e0.x*e0.x + e0.y*e0.y);
				if (magnitude > lcs->maxEigenvectorTensor1) lcs->maxEigenvectorTensor1 = magnitude;
				if (magnitude < lcs->minEigenvectorTensor1) lcs->minEigenvectorTensor1 = magnitude;
				magnitude = sqrt(e1.x*e1.x + e1.y*e1.y);
				if (magnitude > lcs->maxEigenvectorTensor2) lcs->maxEigenvectorTensor2 = magnitude;
				if (magnitude < lcs->minEigenvectorTensor2) lcs->minEigenvectorTensor2 = magnitude;

			}
		}

		//Step2 average tensor_facet in all triangles
		for (int i = 0; i < lcs->triangles.size(); i++) {
			float x0, y0, x1, y1, x2, y2;
			GetAllPointsOfTriangle(i, lcs->rectSize.y, x0, y0, x1, y1, x2, y2);

			lcs->triangles[i].tensor_triangle += lcs->tensor_vertex[x0][y0];
			lcs->triangles[i].tensor_triangle += lcs->tensor_vertex[x1][y1];
			lcs->triangles[i].tensor_triangle += lcs->tensor_vertex[x2][y2];
			lcs->triangles[i].tensor_triangle /= 3;
		}
	}

	

	void CreateInstance(ModelAsset &LineAsset, ModelAsset &RectAsset, std::list<ModelInstance> &gInstances) {
		displayController->CreateInstance(LineAsset, RectAsset, gInstances);
	}
	void UpdateInstance(std::list<ModelInstance> &gInstances, int stepNum) {
		displayController->UpdateInstance(gInstances, grid, nodes, lcs, offline_t, stepNum);
	}

	~CFieldSystem() {
		//if (nodes) delete nodes;
		if (grid)  delete grid;
		if (lcs) delete lcs;
		if (displayController) delete displayController;
	}
};