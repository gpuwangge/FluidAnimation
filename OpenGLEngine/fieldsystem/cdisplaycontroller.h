#pragma once
#include "..\cmathhelper.h"

class CDisplayController {
	int m_width;
	int m_height;
	//float dx, dy, dm;
public:
	VisualScalarType visualScalarType;
	VisualVectorType visualVectorType;
	int downSample;
	bool showUnitVector;

	int testNum;
	int testStep;

	int contourSize;
	//int contourStep;

	CDisplayController(int width, int height) {
		m_width = width;
		m_height = height;

		//float x_min = 0; float x_max = 2;
		//float y_min = 0; float y_max = 1;
		//dx = (x_max - x_min) / (width - 1);
		//dy = (y_max - y_min) / (height - 1);
		//dm = glm::min(dx, dy);

		visualScalarType = SHOW_F;
		visualVectorType = SHOW_NONE;

		downSample = 1;
		showUnitVector = false;

		testNum = 0;
		testStep = 300;

		contourSize = 1000;
		//contourStep = 0;
	}

	void CreateInstance(ModelAsset &LineAsset, ModelAsset &RectAsset, std::list<ModelInstance> &gInstances) {
		for (int i = 0; i < m_height; i++)
			for (int j = 0; j < m_width; j++)
				gInstances.push_back(ModelInstance(&RectAsset, glm::mat4(), glm::vec4(1, 1, 1, 1))); //grid point circles
		for (int i = 0; i < m_height; i++)
			for (int j = 0; j < m_width; j++)
				gInstances.push_back(ModelInstance(&LineAsset, glm::mat4(), glm::vec4(0, 0, 0, 1))); //vector data
		//for (int j = 0; j < testNum; j++)  //aux lines Part I
		//	for (int i = 0; i < testStep; i++) //each grid point has 300 steps
		//		gInstances.push_back(ModelInstance(&LineAsset, glm::mat4(), glm::vec4(0, 0, 0, 1))); //
		for (int j = 0; j < contourSize; j++)  //Contour
			gInstances.push_back(ModelInstance(&LineAsset, glm::mat4(), glm::vec4(0, 0, 0, 1))); //
	}

	glm::vec4 InterpolateColor(float x1, float x, float x2, glm::vec4 c1, glm::vec4 c2) { //x1 < x < x2
		glm::vec4 r = c1 * (x2 - x) / (x2 - x1) + c2 * (x - x1) / (x2 - x1);
		return r;
		//return c1 * (x2 - x) / (x2 - x1) + c2 * (x - x1) / (x2 - x1);
	}


	glm::vec4 GetColor(float x, float min, float max) {
		//float xRed = max; glm::vec4 red(1, 0, 0, 1);
		//float xYellow = (max + min) * 0.75; glm::vec4 yellow(1, 1, 0, 1);
		//float xGreen = (max + min) * 0.5; glm::vec4 green(0, 1, 0, 1);
		//float xSkyBlue = (max + min) * 0.25; glm::vec4 skyBlue(0, 0.75, 1, 1);
		//float xBlue = min; glm::vec4 blue(0, 0, 1, 1);

		//if (x <= xBlue) return blue;
		//if (x > xBlue && x <= xSkyBlue) return InterpolateColor(xBlue, x, xSkyBlue, blue, skyBlue);
		//if (x > xSkyBlue && x <= xGreen) return InterpolateColor(xSkyBlue, x, xGreen, skyBlue, green);
		//if (x > xGreen && x <= xYellow) return InterpolateColor(xGreen, x, xYellow, green, yellow);
		//if (x > xYellow && x <= xRed) return InterpolateColor(xYellow, x, xRed, yellow, red);
		//return red;

		glm::vec4 red(1, 0, 0, 1);
		glm::vec4 green(0, 1, 0, 1);
		glm::vec4 blue(0, 0, 1, 1);


		//float middleValue = (min + max) / 2.0f;
		float valueBlue = min;
		float valueGreen = (min + max) / 2;
		float valueRed = max;

		//float valueRed = 0;
		//float valueGreen = 0;
		//float valueBlue = 0;

		//if (min > 0) { //min and max > 0
		//	valueBlue = min;
		//	valueGreen = (min + max) / 2;
		//	valueRed = max;
		//}
		//else if (max < 0) { //min and max < 0
		//	valueBlue = min;
		//	valueGreen = (min + max) / 2;
		//	valueRed = max;
		//}
		//else { //min < 0 and max > 0
		//	valueBlue = min;
		//	valueGreen = 0;
		//	valueRed = max;
		//}



		if (x <= valueBlue) return blue;
		if (x > valueBlue && x <= valueGreen) return InterpolateColor(valueBlue, x, valueGreen, blue, green);
		if (x > valueGreen && x <= valueRed) return InterpolateColor(valueGreen, x, valueRed, green, red);
		return red;
	}

	//glm::vec4 GetColor(float x) {
	//	//float xRed = 7; glm::vec4 red(1, 0, 0, 1);
	//	//float xYellow = 5; glm::vec4 yellow(1, 1, 0, 1);
	//	//float xGreen = 2; glm::vec4 green(0, 1, 0, 1);
	//	//float xSkyBlue = 0; glm::vec4 skyBlue(0, 0.75, 1, 1);
	//	//float xBlue = -2; glm::vec4 blue(0, 0, 1, 1);

	//	//float xRed = 9; glm::vec4 red(1, 0, 0, 1);
	//	//float xYellow = 7; glm::vec4 yellow(1, 1, 0, 1);
	//	//float xGreen = 5; glm::vec4 green(0, 1, 0, 1);
	//	//float xSkyBlue = 2; glm::vec4 skyBlue(0, 0.75, 1, 1);
	//	//float xBlue = 0; glm::vec4 blue(0, 0, 1, 1);

	//	//float xRed = 500; glm::vec4 red(1, 0, 0, 1);
	//	//float xYellow = 250; glm::vec4 yellow(1, 1, 0, 1);
	//	//float xGreen = 0; glm::vec4 green(0, 1, 0, 1);
	//	//float xSkyBlue = -250; glm::vec4 skyBlue(0, 0.75, 1, 1);
	//	//float xBlue = -500; glm::vec4 blue(0, 0, 1, 1);

	//	//float xRed = 500; glm::vec4 red(1, 0, 0, 1);
	//	//float xYellow = 70; glm::vec4 yellow(1, 1, 0, 1);
	//	//float xGreen = 50; glm::vec4 green(0, 1, 0, 1);
	//	//float xSkyBlue = 20; glm::vec4 skyBlue(0, 0.75, 1, 1);
	//	//float xBlue = 0; glm::vec4 blue(0, 0, 1, 1);

	//	float xRed = 200; glm::vec4 red(1, 0, 0, 1);
	//	float xYellow = 100; glm::vec4 yellow(1, 1, 0, 1);
	//	float xGreen = 0; glm::vec4 green(0, 1, 0, 1);
	//	float xSkyBlue = -100; glm::vec4 skyBlue(0, 0.75, 1, 1);
	//	float xBlue = -200; glm::vec4 blue(0, 0, 1, 1);

	//	if (x <= xBlue) return blue;
	//	if (x > xBlue && x <= xSkyBlue) return InterpolateColor(xBlue, x, xSkyBlue, blue, skyBlue);
	//	if (x > xSkyBlue && x <= xGreen) return InterpolateColor(xSkyBlue, x, xGreen, skyBlue, green);
	//	if (x > xGreen && x <= xYellow) return InterpolateColor(xGreen, x, xYellow, green, yellow);
	//	if (x > xYellow && x <= xRed) return InterpolateColor(xYellow, x, xRed, yellow, red);
	//	return red;
	//}

	void UpdateScalar(int i, int j, std::list<ModelInstance>::iterator &iter, CGrid *grid, CLCS *lcs) {
		switch (visualScalarType) {
		case SHOW_FTLE:
			iter->m_transform = translate(
				grid->location[i][j].x,
				grid->location[i][j].y,
				grid->location[i][j].z)*
				scale(grid->dx, grid->dy, 1);
			if (i < lcs->FTLE.size() && j < lcs->FTLE[0].size())
				iter->m_color = GetColor(lcs->FTLE[i][j], lcs->minFTLE, lcs->maxFTLE);
			else iter->m_color = glm::vec4(1, 1, 1, 1);
			break;
		case SHOW_S:
			iter->m_transform = translate(
				grid->location[i][j].x,
				grid->location[i][j].y,
				grid->location[i][j].z)*
				scale(grid->dx, grid->dy, 1);
			if (i < lcs->S.size() && j < lcs->S[0].size())
				iter->m_color = GetColor(lcs->S[i][j], lcs->minS, lcs->maxS);
			else iter->m_color = glm::vec4(1, 1, 1, 1);
			break;
		case SHOW_F:
			iter->m_transform = translate(
				grid->location[i][j].x,
				grid->location[i][j].y,
				grid->location[i][j].z)*
				scale(grid->dx, grid->dy, 1);
			if (i < lcs->F.size() && j < lcs->F[0].size())
				iter->m_color = GetColor(lcs->F[i][j], lcs->minF, lcs->maxF);
			else iter->m_color = glm::vec4(1, 1, 1, 1);
			break;
			//case SHOW_TESTSCALER:
			//	iter->m_transform = translate(
			//		grid->location[i][j].x,
			//		grid->location[i][j].y,
			//		grid->location[i][j].z)*
			//		scale(grid->dx, grid->dy, 1);
			//	if (i < lcs->testScaler.size() && j < lcs->testScaler[0].size())
			//		iter->m_color = GetColor(lcs->testScaler[i][j] * displaycoff5 * 100.0);
			//	break;
			//case EIGENVALUE_H_MAX:
			//	iter->m_transform = translate(
			//		grid->location[i][j].x,
			//		grid->location[i][j].y,
			//		grid->location[i][j].z)*
			//		scale(grid->dx, grid->dy, 1);
			//	if (i < lcs->eigenvaluesOFHessian.size() && j < lcs->eigenvaluesOFHessian[0].size())
			//		iter->m_color = GetColor(lcs->eigenvaluesOFHessian[i][j][0] * displaycoff4);
			//	break;
			//case EIGENVALUE_H_MIN:
			//	iter->m_transform = translate(
			//		grid->location[i][j].x,
			//		grid->location[i][j].y,
			//		grid->location[i][j].z)*
			//		scale(grid->dx, grid->dy, 1);
			//	if (i < lcs->eigenvaluesOFHessian.size() && j < lcs->eigenvaluesOFHessian[0].size())
			//		iter->m_color = GetColor(lcs->eigenvaluesOFHessian[i][j][1] * displaycoff4);
			//	break;
		default:
			break;
		}
	}

	void UpdateVector(int i, int j, std::list<ModelInstance>::iterator &iter, CGrid *grid, vector<CNodes> &nodes, CLCS *lcs, int stepNum) {
		if (i % downSample != 0 || j % downSample != 0) return;

		glm::vec3 v;
		glm::vec2 values;
		double magnitude = 0;
		switch (visualVectorType) {
		case SHOW_NONE:
			values.x = 0;
			values.y = 0;
			magnitude = 0;
			break;
		case SHOW_VELOCITY:
			values.x = nodes[stepNum].GetVelocity(i, j).x;
			values.y = nodes[stepNum].GetVelocity(i, j).y;
			magnitude = sqrt(values.x * values.x + values.y * values.y) / glm::max(nodes[stepNum].maxVelocity, abs(nodes[stepNum].minVelocity));
			break;
		case SHOW_PHI:
			values.x = lcs->phi[i][j].x - grid->location[i][j].x;
			values.y = lcs->phi[i][j].y - grid->location[i][j].y;
			magnitude = sqrt(values.x * values.x + values.y * values.y) / glm::max(lcs->maxPhi, abs(lcs->minPhi));
			break;
		case SHOW_MAX_EIGENVECTOR_HASSIAN:
			if (i >= lcs->eigenvectorOfHessian.size() || j >= lcs->eigenvectorOfHessian[0].size()) {
				values.x = 0;
				values.y = 0;
				magnitude = 0;
			}
			else {
				values.x = lcs->eigenvectorOfHessian[i][j][0][0];
				values.y = lcs->eigenvectorOfHessian[i][j][1][0];
				magnitude = sqrt(values.x * values.x + values.y * values.y) / glm::max(lcs->maxEigenvectorHessian1, abs(lcs->minEigenvectorHessian1));
			}
			break;
		case SHOW_MIN_EIGENVECTOR_HASSIAN:
			if (i >= lcs->eigenvectorOfHessian.size() || j >= lcs->eigenvectorOfHessian[0].size()) {
				values.x = 0;
				values.y = 0;
				magnitude = 0;
			}
			else {
				values.x = lcs->eigenvectorOfHessian[i][j][0][1];
				values.y = lcs->eigenvectorOfHessian[i][j][1][1];
				magnitude = sqrt(values.x * values.x + values.y * values.y) / glm::max(lcs->maxEigenvectorHessian2, abs(lcs->minEigenvectorHessian2));
			}
			break;
		case SHOW_MAX_EIGENVECTOR_TENSOR:
			if (i >= lcs->eigenvectorOfTensor.size() || j >= lcs->eigenvectorOfTensor[0].size()) {
				values.x = 0;
				values.y = 0;
				magnitude = 0;
			}
			else {
				values.x = lcs->eigenvectorOfTensor[i][j][0][0];
				values.y = lcs->eigenvectorOfTensor[i][j][1][0];
				magnitude = sqrt(values.x * values.x + values.y * values.y) / glm::max(lcs->maxEigenvectorTensor1, abs(lcs->minEigenvectorTensor1));
			}
			break;
		case SHOW_MIN_EIGENVECTOR_TENSOR:
			if (i >= lcs->eigenvectorOfTensor.size() || j >= lcs->eigenvectorOfTensor[0].size()) {
				values.x = 0;
				values.y = 0;
				magnitude = 0;
			}
			else {
				values.x = lcs->eigenvectorOfTensor[i][j][0][1];
				values.y = lcs->eigenvectorOfTensor[i][j][1][1];
				magnitude = sqrt(values.x * values.x + values.y * values.y) / glm::max(lcs->maxEigenvectorTensor2, abs(lcs->minEigenvectorTensor2));
			}
			break;
		case SHOW_CONTOUR:
			values.x = 0;
			values.y = 0;
			magnitude = 0;
			break;
			/*
			case VELOCITY:
			values.x = nodes->velocity[i][j].x;
			values.y = nodes->velocity[i][j].y;
			break;
			case GPHIX:
			if (i >= lcs->gradPhi.size() || j >= lcs->gradPhi[0].size()) break;
			values.x = lcs->gradPhi[i][j][0][0];
			values.y = lcs->gradPhi[i][j][1][0];
			break;
			case GPHIY:
			if (i >= lcs->gradPhi.size() || j >= lcs->gradPhi[0].size()) break;
			values.x = lcs->gradPhi[i][j][0][1];
			values.y = lcs->gradPhi[i][j][1][1];
			break;
			case GH:
			if (i >= lcs->eigenvectorOfHessian.size() || j >= lcs->eigenvectorOfHessian[0].size()) break;
			if (i % 8 != 0 || j % 8 != 0) break;
			values.x = lcs->eigenvectorOfHessian[i][j][0][0];
			values.y = lcs->eigenvectorOfHessian[i][j][1][0];
			break;
			case SHOW_GRAD_FTLE:
			if (i >= lcs->gradFTLE.size() || j >= lcs->gradFTLE[0].size()) break;
			if (i % 8 != 0 || j % 8 != 0) break;
			values.x = lcs->gradFTLE[i][j][0];
			values.y = lcs->gradFTLE[i][j][1];
			break;*/
		default:
			break;
		}

		float dm = glm::min(grid->dx, grid->dy);
		if (showUnitVector && visualVectorType != SHOW_NONE) magnitude = dm * downSample;
		else magnitude = magnitude * dm * downSample;

		glm::mat4 transform = glm::mat4();
		if (!Trivial(magnitude, 0.0000001))
			transform = GetTransform(values.x, values.y);

		iter->m_transform =
			translate(
				grid->location[i][j].x,
				grid->location[i][j].y,
				grid->location[i][j].z) *
			transform *
			scale(magnitude, 1, 1);
	}

	void UpdateContour(int i, std::list<ModelInstance>::iterator &iter, CLCS *lcs) {
		pair<float, float> start(0, 0);
		float magnitude = 0;
		glm::mat4 transform = glm::mat4();
		if (visualVectorType == SHOW_CONTOUR && (i + 1) < lcs->C.size() && i % 2 == 0) {
			start = lcs->C[i];
			pair<float, float> end = lcs->C[i + 1];
			glm::vec2 values(end.first - start.first, end.second - start.second);

			magnitude = sqrt(values.x * values.x + values.y * values.y);
			transform = GetTransform(values.x, values.y);
		}

		iter->m_transform =
			translate(
				start.first,
				start.second,
				0) *
			transform *
			scale(magnitude, 1, 1);
	}

	void UpdateInstance(std::list<ModelInstance> &gInstances, CGrid *grid, vector<CNodes> &nodes, CLCS *lcs,
		float offline_t, int stepNum) {

		int selectIdx = 0;

		//int sy = 30;
		//int sx = 20;

		//float px = grid->location[sy][sx].x;
		//float py = grid->location[sy][sx].y;
		//float pz = grid->location[sy][sx].z;
		//float vx, vy;
		//float fx = px;
		//float fy = py;

		//glm::vec3 selectPos;
		//selectPos.x = grid->location[selectY][selectX].x;
		//selectPos.y = grid->location[selectY][selectX].y;
		//selectPos.z = grid->location[selectY][selectX].z;

		float px = 0;
		float py = 0;
		float pz = 0;

		int idx = 0;
		int block1 = m_height * m_width;
		int block2 = 2 * m_height * m_width;
		for (std::list<ModelInstance>::iterator iter = gInstances.begin(); iter != gInstances.end(); iter++, idx++) {
			int i = idx / m_width;
			int j = idx % m_width;

			if (idx < block1)  UpdateScalar(i, j, iter, grid, lcs);
			else if (idx >= block1 && idx < block2) UpdateVector(i - m_height, j, iter, grid, nodes, lcs, stepNum);
			else UpdateContour(idx - block2, iter, lcs);
			

			//else { //aux lines Part II
			//	i = i - m_height;

			//	//selectIdx: 0~100*300
			//	int vID = selectIdx / testStep; //vertex id: 0~8
			//	int k = selectIdx % testStep; //step num: 0~299

			//	int sy = vID / testNum * 3 + 20; //position index of the vertex
			//	int sx = vID % testNum * 3 + 60;

			//	if (k == 0) {
			//		px = grid->location[sy][sx].x;
			//		py = grid->location[sy][sx].y;
			//		pz = grid->location[sy][sx].z;
			//	}

			//	//glm::vec2 v = GetVelocity(px, py, offline_t + k * lcs->integrat_dt);
			//	glm::vec2 v = GetVelocity_fromFile(px, py, k, stepNum, grid, nodes, sy, sx);

			//	float deltaPx = v.x * lcs->integrat_dt;
			//	float deltaPy = v.y * lcs->integrat_dt;
			//	double magnitude = sqrt(deltaPx * deltaPx + deltaPy * deltaPy);

			//	px += deltaPx;
			//	py += deltaPy;

			//	if (!Trivial(magnitude, 0.0000001)) {
			//		iter->m_transform =
			//			translate(px, py, pz) *
			//			GetTransform(
			//				deltaPx,
			//				deltaPy) *
			//			scale(magnitude, 1, 1);
			//	}
			//	else iter->m_transform = glm::mat4();
			//	/////////////////

			//	selectIdx++;
			//}
			
		}
	}

};