#pragma once

// standard C++ libraries
#include <cassert>
#include <stdexcept>
#include <cmath>
#include <list>
#include <sstream>
#include <iostream>
#include <vector>
#include <unordered_map>

// third-party libraries
#include "..\OpenGLEngine\common\glew\include\GL\glew.h"
#include "..\OpenGLEngine\common\glfw\include\GLFW\glfw3.h"
#include "..\OpenGLEngine\common\glm\glm.hpp"
#include "..\OpenGLEngine\common\glm\gtc\matrix_transform.hpp"

// tdogl classes
#include "..\OpenGLEngine\utilities\Program.h"
#include "..\OpenGLEngine\utilities\Texture.h"
#include "..\OpenGLEngine\utilities\Camera.h"

#include "..\OpenGLEngine\containers.h"

//#include "cmathhelper.h"
#include <math.h>

#include <direct.h>

#define PI 3.1415926
#define TRIVIAL_NUMBER 0.0000001

#define OUTPUT_DATA "OutputData\\"
#define INPUT_DATA "InputData\\"

//#define DEBUG

enum VisualScalarType {
	SHOW_FTLE,
	SHOW_S,
	SHOW_F,
	SHOW_TESTSCALER
};

enum VisualVectorType {
	SHOW_NONE, 
	SHOW_VELOCITY,
	SHOW_PHI,
	SHOW_MAX_EIGENVECTOR_HASSIAN,
	SHOW_MIN_EIGENVECTOR_HASSIAN,
	SHOW_MAX_EIGENVECTOR_TENSOR,
	SHOW_MIN_EIGENVECTOR_TENSOR,
	SHOW_CONTOUR,

	VELOCITY,  GPHIX, GPHIY, GH, EIGENVALUE_H_MAX, EIGENVALUE_H_MIN,
	SHOW_GRAD_FTLE,
};


using namespace std;

