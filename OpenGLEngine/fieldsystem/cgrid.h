#pragma once
//#include "cnodes.h"

class CNodes {
private:
	vector<vector<glm::vec3>> velocity;
	//float m_t;
	int m_step;
public:
	float minVelocity, maxVelocity;

	CNodes(int width, int height, int step) {
		velocity.resize(height, vector<glm::vec3>(width));
		m_step = step;
		//m_t = 0;

		maxVelocity = -99999;
		minVelocity = 99999;
	}
	//void SetTime(float t) {
	//	m_t = t;
	//}

	void SetVelocity(int j, int i, float u, float v, float w) {
		velocity[j][i].x = u;
		velocity[j][i].y = v;
		velocity[j][i].z = w;

		float magnitude = sqrt(u*u + v*v + w*w);
		if (magnitude > maxVelocity) maxVelocity = magnitude;
		if (magnitude < minVelocity) minVelocity = magnitude;
	}
	glm::vec3 GetVelocity(int j, int i) {
		return velocity[j][i];
	}
};

class CTriangle {
public:
	glm::mat2x2 tensor_triangle;
	float area;
	//CTriangle(glm::mat2x2 t) : tensor_triangle(t){
	//	//tensor_triangle[0][0] = 0; tensor_triangle[0][1] = 0;
	//	//tensor_triangle[1][0] = 0; tensor_triangle[1][1] = 0;
	//}

	//to use in marching cube:
	glm::vec3 p[3];
	float f[3];

	CTriangle() {}

	void AssignPosition(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
	}

	void AssignF(float f0, float f1, float f2) {
		f[0] = f0;
		f[1] = f1;
		f[2] = f2;
	}
};

class CHalfedge {
public:
	glm::vec2 start;
	glm::vec2 end;

	float cotan_new;
	float cotan;

	CHalfedge() {}
};


class CGrid {
public:

	vector<vector<glm::vec3>> location;
	double x_min, x_max, y_min, y_max, z_min, z_max;
	double dx, dy, dz;

	CGrid(int width, int height){
		//x_min = -1; x_max = 1;
		//y_min = -1; y_max = 1;
		//z_min = -1; z_max = 1;
		//dx = (x_max - x_min) / (width-1);
		//dy = (y_max - y_min) / (height-1);
		//dz = 0;

		/*
		dx = 0.1;
		dy = 0.1;
		dz = 0;
		float w = (width - 1)* dx;
		float h = (height - 1) * dy;
		float d = 0;
		x_min = -w/2; x_max = w/2;
		y_min = -h/2; y_max = h/2;
		z_min = -d/2; z_max = d/2;*/

		x_min = 0; x_max = 2;
		y_min = 0; y_max = 1;
		z_min = 0; z_max = 0;
		dx = (x_max - x_min) / (width - 1);
		dy = (y_max - y_min) / (height - 1);
		dz = 0;

		location.resize(height, vector<glm::vec3>(width));

		double current_x = x_min, current_y = y_min, current_z = 0;
		for (int i = 0; i < location.size(); i++) {
			current_x = x_min;
			for (int j = 0; j < location[0].size(); j++) {
				location[i][j] = glm::vec3(current_x, current_y, current_z);
				current_x += dx;
			}
			current_y += dy;
		}
	}

	~CGrid(){}

};


