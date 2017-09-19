#define GLM_FORCE_RADIANS

#include "app.h"

#include "loadassets.h"
#include "fieldsystem\cfieldsystem.h"

#include <windows.h>

#include "print\ImageGenerator.h"


// constants
const glm::vec2 SCREEN_SIZE(650, 1440);

// globals
GLFWwindow* gWindow = NULL;
//GLfloat gDegreesRotated = 0.0f;

tdogl::Camera gCamera;
//double gScrollY = 0.0;

ModelAsset gWoodenCrate;
ModelAsset gTriangle;
ModelAsset gFilledTriangle;
ModelAsset gFilledRect;
ModelAsset gLine;
ModelAsset gCircle;
std::list<ModelInstance> gInstances;

//Light gLight;
std::vector<Light> gLights; //Tutorial#8: Even more lights

CFieldSystem *gFieldSystem;

int gStepNum;

//HANDLE hMutex;
//int gOption_disp;


template <typename T>
void SetLightUniform(tdogl::Program* shaders, const char* propertyName, size_t lightIndex, const T& value) {
	std::ostringstream ss;
	ss << "allLights[" << lightIndex << "]." << propertyName;
	std::string uniformName = ss.str();

	shaders->setUniform(uniformName.c_str(), value);
}

//renders a single 'ModelInstance'
static void RenderInstance(const ModelInstance& inst) { //Tutorial#5
	ModelAsset* asset = inst.m_asset;
	tdogl::Program* shaders = asset->shaders;

	//bind the shaders
	shaders->use();

	//set the shader uniforms
	shaders->setUniform("camera", gCamera.matrix());
	shaders->setUniform("model", inst.m_transform);
	shaders->setUniform("color", inst.m_color);

	if (asset->useTexture) {
		shaders->setUniform("materialTex", 0); //set to 0 because the texture will be bound to GL_TEXTURE0, Tutorial#7: More lighting
		shaders->setUniform("materialShininess", asset->shininess);
		shaders->setUniform("materialSpecularColor", asset->specularColor);
		shaders->setUniform("cameraPosition", gCamera.position());
		shaders->setUniform("numLights", (int)gLights.size());

		for (size_t i = 0; i < gLights.size(); ++i) { //Tutorial#8: Even more lighting
			SetLightUniform(shaders, "position", i, gLights[i].position);
			SetLightUniform(shaders, "intensities", i, gLights[i].intensities);
			SetLightUniform(shaders, "attenuation", i, gLights[i].attenuation);
			SetLightUniform(shaders, "ambientCoefficient", i, gLights[i].ambientCoefficient);
			SetLightUniform(shaders, "coneAngle", i, gLights[i].coneAngle);
			SetLightUniform(shaders, "coneDirection", i, gLights[i].coneDirection);
		}
	}

	//bind the texture
	if (asset->useLight) {
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, asset->texture->object());
	}

	//bind VAO and draw
	glBindVertexArray(asset->vao);
	glDrawArrays(asset->drawType, asset->drawStart, asset->drawCount);

	//unbind everything
	glBindVertexArray(0);
	glBindTexture(GL_TEXTURE_2D, 0);
	shaders->stopUsing();
}

static void SaveVariable(string variableName, int scalarType, int vectorType) {
	// clear everything
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	gFieldSystem->displayController->visualScalarType = (VisualScalarType)scalarType;
	gFieldSystem->displayController->visualVectorType = (VisualVectorType)vectorType;
	cout << "Render " + variableName + "..." << endl;
	gFieldSystem->UpdateInstance(gInstances, gStepNum);
	for (std::list<ModelInstance>::const_iterator it = gInstances.begin(); it != gInstances.end(); ++it)
		RenderInstance(*it);
	string filename = "imgs\\" + variableName + "\\" + variableName + to_string(gStepNum) + ".bmp";
	filename = OUTPUT_DATA + filename;
	SaveBMP(filename, SCREEN_SIZE.y, SCREEN_SIZE.x);
}

// draws a single frame
static void Render() {
	glClearColor(1, 1, 1, 1);

	// render all the instances
	SaveVariable("FTLE", SHOW_FTLE, SHOW_NONE);
	SaveVariable("S", SHOW_S, SHOW_NONE);
	SaveVariable("F", SHOW_F, SHOW_NONE);

	SaveVariable("Velocity", SHOW_F, SHOW_VELOCITY);
	SaveVariable("Phi", SHOW_F, SHOW_PHI);
	SaveVariable("Eigenvector_Max_Hassian", SHOW_F, SHOW_MAX_EIGENVECTOR_HASSIAN);
	SaveVariable("Eigenvector_Min_Hassian", SHOW_F, SHOW_MIN_EIGENVECTOR_HASSIAN);
	SaveVariable("Eigenvector_Max_Tensor", SHOW_F, SHOW_MAX_EIGENVECTOR_TENSOR);
	SaveVariable("Eigenvector_Min_Tensor", SHOW_F, SHOW_MIN_EIGENVECTOR_TENSOR);

	SaveVariable("Contour", SHOW_F, SHOW_CONTOUR);

    //swap the display buffers (displays what was just drawn)
    glfwSwapBuffers(gWindow);
}


//void UpdateUI() {
	//move position of camera based on WASD keys, and XZ keys for up and down
	//const float moveSpeed = 2.0; //units per second //Tutorial#4
	//float speed = secondsElapsed * moveSpeed;

	//if (glfwGetKey(gWindow, 'Z'))
	//	gCamera.offsetPosition( speed * -gCamera.forward());
	//else if (glfwGetKey(gWindow, 'X'))
	//	gCamera.offsetPosition(speed * gCamera.forward());
	//if (glfwGetKey(gWindow, 'A'))
	//	gCamera.offsetPosition(speed * -gCamera.right());
	//else if (glfwGetKey(gWindow, 'D'))
	//	gCamera.offsetPosition(speed * gCamera.right());
	//if (glfwGetKey(gWindow, 'S'))
	//	gCamera.offsetPosition(speed * -glm::vec3(0, 1, 0));
	//else if (glfwGetKey(gWindow, 'W'))
	//	gCamera.offsetPosition(speed * glm::vec3(0, 1, 0));

	//if (glfwGetKey(gWindow, '0'))  gFieldSystem->visualDataType = NONE;
	//else if (glfwGetKey(gWindow, '1'))  gFieldSystem->visualDataType = VELOCITY;
	//else if (glfwGetKey(gWindow, '2')) gFieldSystem->visualDataTyp0e = PHI;
	//else if (glfwGetKey(gWindow, '3')) gFieldSystem->visualDataType = GPHIX;
	//else if (glfwGetKey(gWindow, '4')) gFieldSystem->visualDataType = GPHIY;
	//else if (glfwGetKey(gWindow, '5')) gFieldSystem->visualDataType = VFTLE;
	//else if (glfwGetKey(gWindow, '6')) gFieldSystem->visualDataType = SHOW_S;

	//if (glfwGetKey(gWindow, 'Z'))
	//	cout << "z" << endl;

	//int option_disp = 0;
	//cout << "Enter what to render: " << endl;
	//cout << "0: FTLE" << endl;
	//cout << "1: 6x12" << endl;
	//cout << "2: 32x64" << endl;
	//cout << "3: 64x128" << endl;
	//cout << "4: 128x256" << endl;
	//cout << "5: 256x512" << endl;
//}


// update the scene based on the time elapsed since last update
void Update(float secondsElapsed, int integral_step) {
	cout << "Processing step " << gStepNum << endl;
	//************************//
	cout << "Update system..." << endl;
	gFieldSystem->Update(secondsElapsed, gStepNum, integral_step);
	//************************//

	//gLights[0].position = glm::vec4(gCamera.position(), 1.0);
	//gLights[0].coneDirection = gCamera.forward();
	//gLights[0].intensities = glm::vec3(1, 0, 0); //red
	//gLights[0].intensities = glm::vec3(0, 1, 0); //green
	//gLights[0].intensities = glm::vec3(1, 1, 1); //white

	//rotate camera based on mouse movement
	//const float mouseSensitivity = 0.1f; //Tutorial#4
	//double mouseX, mouseY;
	//glfwGetCursorPos(gWindow, &mouseX, &mouseY);
	//gCamera.offsetOrientation(mouseSensitivity * (float)mouseY, mouseSensitivity * (float)mouseX);
	//glfwSetCursorPos(gWindow, 0, 0); //reset the mouse, so it doesn't go out of the window

	//increase or decrease field of view based on mouse wheel
	//const float zoomSensitivity = -0.2f; //Tutorial#4
	//float fieldOfView = gCamera.fieldOfView() + zoomSensitivity * (float)gScrollY;
	//if (fieldOfView < 5.0f) fieldOfView = 5.0f;
	//if (fieldOfView > 130.0f) fieldOfView = 130.0f;
	//gCamera.setFieldOfView(fieldOfView);
	//gScrollY = 0;


}


//DWORD WINAPI SubThread(LPVOID lpParamter) {
//	while (1) {
//		//WaitForSingleObject(hMutex, INFINITE);
//
//		//gThreadData += 1000;
//		//cout << "Sub Thread!" << gThreadData << endl;
//
//		Sleep(500);
//		//cout << "2";
//		UpdateUI();
//
//		//ReleaseMutex(hMutex);
//	}
//}


void OnError(int errorCode, const char* msg) { throw std::runtime_error(msg); }

// the program starts here
void AppMain() {
	// initialise GLFW
	glfwSetErrorCallback(OnError);
	if (!glfwInit())
		throw std::runtime_error("glfwInit failed");

	// open a window with GLFW
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	gWindow = glfwCreateWindow((int)SCREEN_SIZE.y, (int)SCREEN_SIZE.x, "LCS Generator", NULL, NULL);
	glfwSetWindowPos(gWindow, 700, 100);
	if (!gWindow)
		throw std::runtime_error("glfwCreateWindow failed. Can your hardware handle OpenGL 3.2?");

	// GLFW settings
	//glfwSetInputMode(gWindow, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	//glfwSetCursorPos(gWindow, 0, 0);
	//glfwSetScrollCallback(gWindow, OnScroll);
	glfwMakeContextCurrent(gWindow);

	// initialise GLEW
	glewExperimental = GL_TRUE; //stops glew crashing on OSX :-/
	if (glewInit() != GLEW_OK)
		throw std::runtime_error("glewInit failed");

	// GLEW throws some errors, so discard all the errors so far
	while (glGetError() != GL_NO_ERROR) {} //Tutorial#3

	// print out some info about the graphics drivers
	std::cout << "OpenGL version: " << glGetString(GL_VERSION) << std::endl;
	std::cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
	std::cout << "Vendor: " << glGetString(GL_VENDOR) << std::endl;
	std::cout << "Renderer: " << glGetString(GL_RENDERER) << std::endl;
	cout << endl;

	// make sure OpenGL version 3.2 API is available
	if (!GLEW_VERSION_3_2)
		throw std::runtime_error("OpenGL 3.2 API is not available.");

	// OpenGL settings
	glEnable(GL_DEPTH_TEST); //Tutorial#3
	glDepthFunc(GL_LESS);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//************************//
	gStepNum = 0;
	string filename;

	filename = OUTPUT_DATA; _mkdir(filename.c_str());

	filename = OUTPUT_DATA + (string)"imgs"; _mkdir(filename.c_str());

	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"MatlabContour"; _mkdir(filename.c_str());

	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"FTLE"; _mkdir(filename.c_str());
	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"S"; _mkdir(filename.c_str());
	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"F"; _mkdir(filename.c_str());

	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"Velocity"; _mkdir(filename.c_str());
	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"Phi"; _mkdir(filename.c_str());
	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"Eigenvector_Max_Hassian"; _mkdir(filename.c_str());
	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"Eigenvector_Min_Hassian"; _mkdir(filename.c_str());
	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"Eigenvector_Max_Tensor"; _mkdir(filename.c_str());
	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"Eigenvector_Min_Tensor"; _mkdir(filename.c_str());


	filename = OUTPUT_DATA + (string)"imgs\\" + (string)"Contour"; _mkdir(filename.c_str());
	//************************//
	CPrint printer;
	printer.OpenForRead_AscII(INPUT_DATA + (string)"config.txt");
	int option_res = printer.ReadASCIIFile_One_Number();
	int option_maxstep = printer.ReadASCIIFile_One_Number();
	bool option_pause = printer.ReadASCIIFile_One_Number();
	//************************//
	int width = 0, height = 0;
	switch (option_res) {
	case 0:
		width = 8 + 1; height = 4 + 1;
		//gFieldSystem = new CFieldSystem(8 + 1, 4 + 1, option_maxstep);
		break;
	case 1:
		width = 12 + 1; height = 6 + 1;
		//gFieldSystem = new CFieldSystem(12 + 1, 6 + 1, option_maxstep);
		break;
	case 2:
		width = 64 + 1; height = 32 + 1;
		//gFieldSystem = new CFieldSystem(64 + 1, 32 + 1, option_maxstep);
		break;
	case 3:
		width = 128 + 1; height = 64 + 1;
		//gFieldSystem = new CFieldSystem(128 + 1, 64 + 1, option_maxstep);
		break;
	case 4:
		width = 256 + 1; height = 128 + 1;
		//gFieldSystem = new CFieldSystem(256 + 1, 128 + 1, option_maxstep);
		break;
	case 5:
		width = 512 + 1; height = 256 + 1;
		//gFieldSystem = new CFieldSystem(512 + 1, 256 + 1, option_maxstep);
		break;
	case 6:
		width = 50 + 1; height = 25 + 1;
		//gFieldSystem = new CFieldSystem(50 + 1, 25 + 1, option_maxstep);
		break;
	}

	gFieldSystem = new CFieldSystem(width, height, option_maxstep);

	gFieldSystem->displayController->downSample = printer.ReadASCIIFile_One_Number();
	gFieldSystem->displayController->showUnitVector = printer.ReadASCIIFile_One_Number();
	printer.CloseFile();
	//************************//
	filename = "Parameters.txt";
	printer.OpenForRead_AscII(INPUT_DATA + filename);
	gFieldSystem->lcs->integrat_T = printer.ReadASCIIFile_One_Number();
	gFieldSystem->lcs->integrat_dt = printer.ReadASCIIFile_One_Number();
	gFieldSystem->epsilon_s = printer.ReadASCIIFile_One_Number();
	gFieldSystem->thresh_s = printer.ReadASCIIFile_One_Number();
	gFieldSystem->epsilon_tensor = printer.ReadASCIIFile_One_Number();
	gFieldSystem->lcs->weight_b = printer.ReadASCIIFile_One_Number();
	gFieldSystem->lcs->smooth_h = printer.ReadASCIIFile_One_Number();
	gFieldSystem->lcs->smooth_n = printer.ReadASCIIFile_One_Number();
	gFieldSystem->threshold_FTLE_s = printer.ReadASCIIFile_One_Number();
	printer.CloseFile();
	//************************//
	if (option_maxstep > 0) {
		filename = "airfoil-velocity.dat";
		printer.OpenForRead_AscII(INPUT_DATA + filename);
		string test_s;
		for (int i = 0; i < 10; i++) //Overview header
			printer.fs >> test_s;

		float input_x;
		float input_y;
		float input_u;
		float input_v;
		float start_t = 1.0;

		int h = 26;
		int w = 51;
		

		for (int k = 0; k < option_maxstep; k++) {
			for (int i = 0; i < 6; i++) printer.fs >> test_s; //Frame header
			
			vector<vector<glm::vec3>> tmpVel;
			tmpVel.resize(h, vector<glm::vec3>(w));

			for (int i = 0; i < h * w; i++) {
				printer.fs >> input_x;//not used
				printer.fs >> input_y;//not used
				printer.fs >> input_u;
				printer.fs >> input_v;

				//gFieldSystem->nodes[k].SetTime(start_t + k * 0.1);
				//gFieldSystem->nodes[k].SetVelocity(i / 51, i % 51, input_u, input_v, 0);
				tmpVel[i / w][i % w].x = input_u;
				tmpVel[i / w][i % w].y = input_v;
			}

			float tmp_di = 1.0 / (h - 1);
			float tmp_dj = 2.0 / (w - 1);
			
			float di = 1.0 / (height - 1);
			float dj = 2.0 / (width - 1);
			

			gFieldSystem->nodes.push_back(CNodes(width, height, k));
			for (int i = 0; i < height - 1; i++)
				for (int j = 0; j < width - 1; j++) {
					float pi = i * di;
					float pj = j * dj;

					int tmp_i = pi / tmp_di;
					int tmp_j = pj / tmp_dj;
					//float tmp_ri = pi - tmp_i * tmp_di;
					//float tmp_rj = pj - tmp_j * tmp_dj;

					//Bilinear interpolation
					glm::vec3 f11 = tmpVel[tmp_i][tmp_j];			glm::vec3 f12 = tmpVel[tmp_i][tmp_j + 1];
					glm::vec3 f21 = tmpVel[tmp_i + 1][tmp_j];		glm::vec3 f22 = tmpVel[tmp_i + 1][tmp_j + 1];

					glm::vec3 r = BilinearInterpolation(
						f11, f12, f21, f22,
						tmp_j*tmp_dj, tmp_j*tmp_dj + tmp_dj, tmp_i*tmp_di, tmp_i*tmp_di + tmp_di,
						pj, pi
					);

					float coff = 2.0f; //because the grid height and length is half in the example
					gFieldSystem->nodes[k].SetVelocity(i, j, r.x * coff, r.y * coff, 0);

				}



		}

		printer.CloseFile();
	}
	//************************//

	// initialise the gWoodenCrate asset
	//LoadWoodenCrateAsset(gWoodenCrate); // Tutorial#5
	//LoadTriangleAsset(gTriangle);
	//LoadFilledTriangleAsset(gFilledTriangle);
	LoadLineAsset(gLine);
	LoadCircleAsset(gCircle);
	LoadFilledRectAsset(gFilledRect);
	//LoadFieldSystemAsset(gLine, fieldSystem);

	//************************//
	gFieldSystem->CreateInstance(gLine, gFilledRect, gInstances);

	//gInstances.push_back(ModelInstance(&gWoodenCrate, glm::mat4(), glm::vec4(0, 0, 0, 1))); //test crate
	//gInstances.push_back(ModelInstance(&gTriangle, glm::mat4(), glm::vec4(0, 0, 0, 1))); //test triangle
	//gInstances.push_back(ModelInstance(&gFilledTriangle, glm::mat4(), glm::vec4(1, 0, 0, 1))); //test triangle
	//gInstances.push_back(ModelInstance(&gCircle, translate(1, 1, 0), glm::vec4(1, 0, 0, 1))); //test circle
	//gInstances.push_back(ModelInstance(&gFilledRect, translate(0.5, 0.5, 0), glm::vec4(1, 0, 0, 1))); //test rect
	//************************//

	// setup gCamera
	gCamera.setPosition(glm::vec3(1.0, 0.5, 1.25)); //Tutorial#4: Camera
	gCamera.setViewportAspectRatio(SCREEN_SIZE.y / SCREEN_SIZE.x);
	gCamera.setNearAndFarPlanes(0.5f, 100.0f); //Tutorial#7: More lighting

	// setup gLight
	Light spotlight; //Tutorial#8: Even more lighting
	spotlight.position = glm::vec4(-4, 0, 10, 1);
	spotlight.intensities = glm::vec3(2, 2, 2); //strong white light
	spotlight.attenuation = 0.1f;
	spotlight.ambientCoefficient = 0.0f; //no ambient light
	spotlight.coneAngle = 15.0f;
	spotlight.coneDirection = glm::vec3(0, 0, -1);

	Light directionalLight;
	directionalLight.position = glm::vec4(1, 0.8, 0.6, 0); //w == 0 indications a directional light
	directionalLight.intensities = glm::vec3(0.4, 0.3, 0.1); //weak yellowish light
	directionalLight.ambientCoefficient = 0.06f;

	gLights.push_back(spotlight);
	gLights.push_back(directionalLight);

	//Create second thread
	//HANDLE hThread = CreateThread(NULL, 0, SubThread, NULL, 0, NULL);
	//hMutex = CreateMutex(NULL, FALSE, NULL);
	//gOption_disp = 0;
	//CloseHandle(hThread);
	
	//while (1) {
		//WaitForSingleObject(hMutex, INFINITE);

		//gThreadData += 1;
		//cout << "Main Thread!" << gThreadData << endl;

		//Sleep(1000);

		//ReleaseMutex(hMutex);
		//cout << "1";
	//}

	int integral_step = abs(gFieldSystem->lcs->integrat_T) / gFieldSystem->lcs->integrat_dt;

    // run while the window is open
	double lastTime = glfwGetTime();
	//double lastTime2 = glfwGetTime();

    while(!glfwWindowShouldClose(gWindow)){
        // process pending events
        glfwPollEvents();

		// update the scene based on the time elapsed since last update
		double thisTime = glfwGetTime(); //Tutorial#3
		double interval = thisTime - lastTime;
		//UpdateUI(thisTime - lastTime2);
		//lastTime2 = thisTime;
		
		//UpdateUI();
		//UpdateUI(interval);

		//if (interval < 5) continue;

		
		if (gFieldSystem->lcs->integrat_T < 0 && gStepNum < integral_step) {
			cout << "Skip step " << gStepNum << "." << endl;
		}
		else {
			Update(interval, integral_step);
			// draw one frame
			//cout << "Rendering..." << endl;
			Render();
		}

		lastTime = thisTime;
		gStepNum++;

		cout << "Done!" << endl << endl;

		if(option_pause){
			cout << "Press ENTER to process next step..." << endl;
			getchar();
		}

		// check for errors
		GLenum error = glGetError(); //Tutorial#3
		if (error != GL_NO_ERROR)
			std::cerr << "OpenGL Error " << error << std::endl;

		//exit program if escape key is pressed
		if (glfwGetKey(gWindow, GLFW_KEY_ESCAPE) || (option_maxstep!= 0 && gStepNum >= option_maxstep)) { //Tutorial#4
			if(gFieldSystem) delete gFieldSystem;
			glfwSetWindowShouldClose(gWindow, GL_TRUE);
		}

    }

    // clean up and exit
    glfwTerminate();
}



int main(int argc, char *argv[]) {
    try {
        AppMain();
    } catch (const std::exception& e){
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
