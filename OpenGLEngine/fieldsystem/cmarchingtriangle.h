#pragma once
#include "..\matlab\MatlabEng.h"

class CMarchingTriangle {
public:
	CMarchingTriangle() {

	}

	bool MarchingCubeIntersect(float isovalue, float f1, float f2, float &d1, float &d2) {
		d1 = abs(f1 - isovalue);
		d2 = abs(f2 - isovalue);
		if ((f1 < isovalue) && (f2 > isovalue)) return true;
		if ((f1 > isovalue) && (f2 < isovalue)) return true;
		return false;
	}

	void Contour(vector<vector<float>> &F, float isovalue, vector<pair<float,float>> &C, vector<CTriangle> &triangles, string filename) { //, CGrid *grid
#ifdef DEBUG
		CPrint printer;
		printer.CreateFile(OUTPUT_DATA + filename);
#endif

		C.clear();

		int height = F.size();//62
		int width = F[0].size();//126

		//int s = triangles.size();//61*125*2=15250

		//Assign F to triangle
		int rectSizey = width - 1;
		for (int i = 0; i < triangles.size(); i++) {
			float x0, y0, x1, y1, x2, y2;
			GetAllPointsOfTriangle(i, rectSizey, x0, y0, x1, y1, x2, y2);  //Order is start from the up-left point, counter-clock wise

			triangles[i].AssignF(F[x0][y0], F[x1][y1], F[x2][y2]);//need double check
		}

		//Update C
		float d1 = 0;
		float d2 = 0;
		for (int i = 0; i < triangles.size(); i++) {
			//edge 01
			if (MarchingCubeIntersect(isovalue, triangles[i].f[0], triangles[i].f[1], d1, d2))
				C.push_back(make_pair<float, float>(
					(triangles[i].p[0].x * d2 + triangles[i].p[1].x * d1) / (d1 + d2),
					(triangles[i].p[0].y * d2 + triangles[i].p[1].y * d1) / (d1 + d2)));

			//edge 12
			if (MarchingCubeIntersect(isovalue, triangles[i].f[1], triangles[i].f[2], d1, d2))
				C.push_back(make_pair<float, float>(
					(triangles[i].p[1].x * d2 + triangles[i].p[2].x * d1) / (d1 + d2),
					(triangles[i].p[1].y * d2 + triangles[i].p[2].y * d1) / (d1 + d2)));

			//edge 20
			if (MarchingCubeIntersect(isovalue, triangles[i].f[2], triangles[i].f[0], d1, d2))
				C.push_back(make_pair<float, float>(
					(triangles[i].p[2].x * d2 + triangles[i].p[0].x * d1) / (d1 + d2),
					(triangles[i].p[2].y * d2 + triangles[i].p[0].y * d1) / (d1 + d2)));
		}

#ifdef DEBUG
		for (int i = 0; i < C.size(); i++){
			printer.PrintToFile("i = ");
			printer.PrintToFile(i, 0);
			printer.PrintToFile(": ");
			printer.PrintToFile(C[i].first);
			printer.PrintToFile(", ");
			printer.PrintToFile_Enter(C[i].second);
		}

		printer.CloseFile();
#endif
	}


	void MatlabContour(CMatlabEng &matlabEng, vector<vector<float>> &F, float isovalue, int stepNum, string currentDirectory) {
		int height = F.size();
		int width = F[0].size();
		std:string command;

		std::ostringstream WIDTH; WIDTH << width;
		std::ostringstream HEIGHT; HEIGHT << height;
		command = "F = reshape(V, " + WIDTH.str() + ", " + HEIGHT.str() + ");"; matlabEng.EvalString(command.c_str());
		command = "F = F';"; matlabEng.EvalString(command.c_str());
		std::ostringstream ISOVALUE; ISOVALUE << isovalue;
		command = "[C,h] = contour(F, [" + ISOVALUE.str() + " " + ISOVALUE.str() + "])"; matlabEng.EvalString(command.c_str());

		float ratio = width / height;
		std::ostringstream RATIO; RATIO << ratio;
		command = "pbaspect([" + RATIO.str() + " 1 1])"; matlabEng.EvalString(command.c_str());

		std::ostringstream STEPNUM; STEPNUM << stepNum;
		command = "saveas(h, '" + currentDirectory + "\\" + OUTPUT_DATA + "imgs\\MatlabContour\\LCS" + STEPNUM.str() + ".png')";
		matlabEng.EvalString(command.c_str());
	}

};