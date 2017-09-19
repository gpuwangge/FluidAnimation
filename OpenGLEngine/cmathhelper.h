#pragma once
#include "app.h"

// convenience function that returns a translation matrix
glm::mat4 translate(GLfloat x, GLfloat y, GLfloat z) { //Tutorial#5
	return glm::translate(glm::mat4(), glm::vec3(x, y, z));
}

// convenience function that returns a scaling matrix
glm::mat4 scale(GLfloat x, GLfloat y, GLfloat z) { //Tutorial#5
	return glm::scale(glm::mat4(), glm::vec3(x, y, z));
}

glm::mat4 rotate(GLfloat angle, glm::vec3 axis) { 
	return glm::rotate(glm::mat4(), angle, axis);
}

bool Trivial(double x, float parameter) {
	if (abs(x) < parameter)
		return true;
	return false;
}

glm::mat4 GetTransform(double x, double y) {
	if (Trivial(x, 0.00001) && Trivial(y, 0.00001))
		return glm::mat4();

	if(Trivial(x, 0.00001)) {
		if(y > 0) return glm::rotate(glm::mat4(), (GLfloat)PI/2, glm::vec3(0, 0, 1));
		else return glm::rotate(glm::mat4(), -(GLfloat)PI / 2, glm::vec3(0, 0, 1));
	}

	if (Trivial(y, 0.00001)) {
		if (x > 0) return glm::rotate(glm::mat4(), (GLfloat)0, glm::vec3(0, 0, 1));
		else return glm::rotate(glm::mat4(), (GLfloat)PI, glm::vec3(0, 0, 1));
	}

	if((x > 0 && y > 0) || (x > 0 && y < 0))
		return glm::rotate(glm::mat4(), (GLfloat)(glm::atan(y / x)), glm::vec3(0, 0, 1));
	else 
		return glm::rotate(glm::mat4(), (GLfloat)(PI + glm::atan(y / x)), glm::vec3(0, 0, 1));
}


bool QuadraticEq(double a, double b, double c, float &x1, float &x2) {
	if (Trivial(a, 0.00001)) { cout << "Quadratic faild! (a = 0)" << endl; return false; }
	float delta = b * b - 4 * a * c;
	if (Trivial(delta, 0.00001)) delta = 0;//??????????????????
	if (delta < 0) {
		cout << "Quadratic faild! (delta = " << delta << ")" << endl;;
		return false; 
	}else if (Trivial(delta, 0.00001)) {
		x1 = -1 * b / (2.0f * a);
		x2 = x1;
	}
	else {
		x1 = (-1 * b + sqrt(delta)) / (2.0f * a);
		x2 = (-1 * b - sqrt(delta)) / (2.0f * a);
	}
	return true;
}

//Get all the eigenvalue
// r is from bigger to smaller
glm::vec2 GetEigenvalue(glm::mat2x2 &x){
	//double testb00 = x[0][0];
	//double testb01 = x[0][1];
	//double testb10 = x[1][0];
	//double testb11 = x[1][1];

	double a = 1.0f;
	double b = -(x[0][0] + x[1][1]);
	double c = x[0][0] * x[1][1] - (x[0][1] * x[1][0]);

	glm::vec2 r(0,0);

	if (QuadraticEq( a,b,c, r[0], r[1])){
		//if (abs(r[0]) < abs(r[1])){
		if (r[0] < r[1]) {
			float t = r[0];
			r[0] = r[1];
			r[1] = t;
		}

	}
	else cout << "Error: Eigenvalue false." << endl;
	
	return r;

}


bool GetEigenInfo(glm::mat2x2 &x, glm::mat2x2 &e, glm::vec2 &lambda) {
	//! right now only works if m is a 2 X 2 matrix
	//since eigenvalues are in decrease order, the first eigenvector is the larger, the second is smaller
	//double *eigenvalue;
	//eigenvalue = new double[2];
	//SMatrix *tmp = new SMatrix(2, 2);
	lambda = GetEigenvalue(x);

	for (int i = 0; i < 2; i++) {
		float a = float(lambda[i] -x[0][0]);
		float b = float(0 - x[0][1]);
		float c = float(0 - x[1][0]);
		float d = float(lambda[i] - x[1][1]);

		float smallValue = 0.1f;

		if (abs(a * d - b * c) <= smallValue) { //singular, non-invertible, has non-zero result x
			if (abs(a) > smallValue) {
				e[0][i] = -1 * b / a;
				e[1][i] = 1;
			}else if (abs(b)> smallValue) {
				e[0][i] = 1;
				e[1][i] = -1 * a / b;
			}else if (abs(c)> smallValue) {
				e[0][i] = -1 * d / c;
				e[1][i] = 1;
			}else if (abs(d)> smallValue) {
				e[0][i] = 1;
				e[1][i] = -1 * c / d;
			}else{
				e[0][i] = 1;
				e[1][i] = 1;
			}
			float norm = sqrt(e[0][i] * e[0][i] + e[1][i] * e[1][i]);
			e[0][i] = e[0][i] / norm;
			e[1][i] = e[1][i] / norm;
		}
		else { //non-singular, invertible, has zero result x
			e[0][i] = 0;
			e[1][i] = 0;
		}
	}

	//if(eigenvalue[i] < 0){ //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//	r->SetValue(0,i, 0);
	//	r->SetValue(1,i, 0);	
	//}

	//delete eigenvalue;
	return true;
	//}
	//else {
	//	delete tmp;
		//delete eigenvalue;
		//return false;
	//}
}



float dotProduct(glm::vec2 x, glm::vec2 y) {
	return x[0] * y[0] + x[1] * y[1];
}


float f(float x, float t, float &dfdx) {
	float epsilon =  0.25f; //if epsilon (perturbation strength) is zero, velocity is time independant
	float w = 2 * PI / 10.0f;

	float a = epsilon * sin(w * t);
	float b = 1 - 2 * epsilon * sin(w * t);

	dfdx = 2 * a * x + b;

	return a * x * x + b * x;
}




//Conclusion: Given grid start from (si,sj), for triangle[i], look for these:
//The first point is (si+m,sj+n)
//if(i%2==0) look for (si+1+m,sj+n)
//else look for (si+m,sj+1+n)
//The last point is (si+1+m,sj+1+n)
//Given m=i/k, n=(i%k)/2, k=(width-4)*2 
pair<int, int> GetFirstPointOfTriangle(int i, int rectSize) {
	int k = rectSize * 2;
	int m = i / k;
	int n = (i % k) / 2;
	return make_pair<int, int>(m+0, n+0);
}
//Order is start from the up-left point, counter-clock wise
void GetAllPointsOfTriangle(int i, int rectSize, float &x0, float &y0, float &x1, float &y1, float &x2, float &y2) {
	pair<int, int> pt = GetFirstPointOfTriangle(i, rectSize);
	x0 = pt.first;
	y0 = pt.second;
	
	if (i % 2 == 0) { //down
		x1 = x0 + 1;
		y1 = y0;

		x2 = x0 + 1;
		y2 = y0 + 1;
	}
	else {//up
		x2 = x0;
		y2 = y0 + 1;

		x1 = x0 + 1;
		y1 = y0 + 1;
	}
}


//int GetFirstOneRingNeighbor(int i, int j, int si, int sj, int k) {
//	return (i - si) * k + 2 * (j - sj);
//	
//}

int GetIndexOfPoint(int i, int j, int width) {
	return i * width + j;
}

void GetPointByIndex(int idx, int width, int &i, int &j) {
	i = idx / width;
	j = idx % width;
}

inline float getLength(glm::vec3 &v) {
	return sqrt(v.x * v.x + v.y * v.y + v.z*v.z);
}


glm::vec2 GetVelocity(float px, float py, float t) {
	glm::vec2 v;

	float A = 0.1;
	float dfdx;
	float fx = f(px, t, dfdx);
	v.x = (float)(-PI * A * sin(PI * fx) * cos(PI * py));
	v.y = (float)(PI * A * cos(PI * fx) * sin(PI * py) * dfdx);

	return v;
}

glm::vec3 BilinearInterpolation(glm::vec3 f11, glm::vec3 f12, glm::vec3 f21, glm::vec3 f22,
	float x1, float x2, float y1, float y2, float x, float y) {

	float dx = x2 - x1;
	float dy = y2 - y1;
	float rx = x - x1;
	float ry = y - y1;

	float coffx_1 = (dx - rx) / dx;
	float coffx_2 = rx / dx;
	glm::vec3 r1 = f11 * coffx_1 + f12 * coffx_2;
	glm::vec3 r2 = f21 * coffx_1 + f22 * coffx_2;

	float coffy_1 = (dy - ry) / dy;
	float coffy_2 = ry / dy;
	glm::vec3 r = r1 * coffy_1 + r2 * coffy_2;

	return r;
}

glm::vec2 GetVelocity_fromFile(float px, float py, int k, int step, CGrid *grid, vector<CNodes> &nodes, int i, int j) {
	glm::vec2 v;

	//no location interpolation
	//v.x = nodes[step - k].GetVelocity(i, j).x;
	//v.y = nodes[step - k].GetVelocity(i, j).y;

	
	//int x = px / grid->dx;//width  - x
	//int y = py / grid->dy;//height  - y
	//if (x >= 0 && y >= 0 && x < grid->location[0].size() - 1 && y < grid->location.size() - 1) {
	//	v.x = nodes[step - k].GetVelocity(y, x).x;
	//	v.y = nodes[step - k].GetVelocity(y, x).y;
	//}
	
	/*Interpolation*/
	int tmp_i = py / grid->dy;
	int tmp_j = px / grid->dx;
	
	if (tmp_j >= 0 && tmp_i >= 0 && tmp_j < grid->location[0].size() - 1 && tmp_i < grid->location.size() - 1) {
		glm::vec3 f11 = nodes[step - k].GetVelocity(tmp_i, tmp_j);			glm::vec3 f12 = nodes[step - k].GetVelocity(tmp_i, tmp_j + 1);
		glm::vec3 f21 = nodes[step - k].GetVelocity(tmp_i + 1, tmp_j); 		glm::vec3 f22 = nodes[step - k].GetVelocity(tmp_i + 1, tmp_j + 1); ;

		glm::vec3 r = BilinearInterpolation(
			f11, f12, f21, f22,
			tmp_j*grid->dx, tmp_j*grid->dx + grid->dx, tmp_i*grid->dy, tmp_i*grid->dy + grid->dy,
			px, py
		);
		v.x = r.x;
		v.y = r.y;
	}



	return v;
}


