#pragma once

static tdogl::Program* LoadShaders(const char* vertFilename, const char* fragFilename) { //Tutorial#5
	std::vector<tdogl::Shader> shaders;
	shaders.push_back(tdogl::Shader::shaderFromFile("..\\OpenGLEngine\\shaders\\" + std::string(vertFilename), GL_VERTEX_SHADER));
	shaders.push_back(tdogl::Shader::shaderFromFile("..\\OpenGLEngine\\shaders\\" + std::string(fragFilename), GL_FRAGMENT_SHADER));
	return new tdogl::Program(shaders);
}

// returns a new tdogl::Texture created from the given filename
static tdogl::Texture* LoadTexture(const char* filename) { //Tutorial#5
	tdogl::Bitmap bmp = tdogl::Bitmap::bitmapFromFile("..\\OpenGLEngine\\resources\\" + std::string(filename));
	bmp.flipVertically();
	return new tdogl::Texture(bmp);
}

static void LoadWoodenCrateAsset(ModelAsset &asset) {
	// set all the elements of gWoodenCrate
	asset.shaders = LoadShaders("vertex-shader-wooden-crate.txt", "fragment-shader-wooden-crate.txt"); //Tutorial#5
	asset.drawType = GL_TRIANGLES;
	asset.drawStart = 0;
	asset.drawCount = 6 * 2 * 3;
	asset.texture = LoadTexture("wooden-crate.jpg");

	asset.shininess = 80.0; //Tutorial#7: More lighting
	asset.specularColor = glm::vec3(1.0f, 1.0f, 1.0f);

	glGenBuffers(1, &asset.vbo);
	glGenVertexArrays(1, &asset.vao);
	// bind the VAO
	glBindVertexArray(asset.vao);
	// bind the VBO
	glBindBuffer(GL_ARRAY_BUFFER, asset.vbo);

	GLfloat vertexData[] = {	//Tutorial#6: Lighting
								//  X     Y     Z       U     V          Normal
								// bottom
		-1.0f,-1.0f,-1.0f,   0.0f, 0.0f,   0.0f, -1.0f, 0.0f,
		1.0f,-1.0f,-1.0f,   1.0f, 0.0f,   0.0f, -1.0f, 0.0f,
		-1.0f,-1.0f, 1.0f,   0.0f, 1.0f,   0.0f, -1.0f, 0.0f,
		1.0f,-1.0f,-1.0f,   1.0f, 0.0f,   0.0f, -1.0f, 0.0f,
		1.0f,-1.0f, 1.0f,   1.0f, 1.0f,   0.0f, -1.0f, 0.0f,
		-1.0f,-1.0f, 1.0f,   0.0f, 1.0f,   0.0f, -1.0f, 0.0f,

		// top
		-1.0f, 1.0f,-1.0f,   0.0f, 0.0f,   0.0f, 1.0f, 0.0f,
		-1.0f, 1.0f, 1.0f,   0.0f, 1.0f,   0.0f, 1.0f, 0.0f,
		1.0f, 1.0f,-1.0f,   1.0f, 0.0f,   0.0f, 1.0f, 0.0f,
		1.0f, 1.0f,-1.0f,   1.0f, 0.0f,   0.0f, 1.0f, 0.0f,
		-1.0f, 1.0f, 1.0f,   0.0f, 1.0f,   0.0f, 1.0f, 0.0f,
		1.0f, 1.0f, 1.0f,   1.0f, 1.0f,   0.0f, 1.0f, 0.0f,

		// front
		-1.0f,-1.0f, 1.0f,   1.0f, 0.0f,   0.0f, 0.0f, 1.0f,
		1.0f,-1.0f, 1.0f,   0.0f, 0.0f,   0.0f, 0.0f, 1.0f,
		-1.0f, 1.0f, 1.0f,   1.0f, 1.0f,   0.0f, 0.0f, 1.0f,
		1.0f,-1.0f, 1.0f,   0.0f, 0.0f,   0.0f, 0.0f, 1.0f,
		1.0f, 1.0f, 1.0f,   0.0f, 1.0f,   0.0f, 0.0f, 1.0f,
		-1.0f, 1.0f, 1.0f,   1.0f, 1.0f,   0.0f, 0.0f, 1.0f,

		// back
		-1.0f,-1.0f,-1.0f,   0.0f, 0.0f,   0.0f, 0.0f, -1.0f,
		-1.0f, 1.0f,-1.0f,   0.0f, 1.0f,   0.0f, 0.0f, -1.0f,
		1.0f,-1.0f,-1.0f,   1.0f, 0.0f,   0.0f, 0.0f, -1.0f,
		1.0f,-1.0f,-1.0f,   1.0f, 0.0f,   0.0f, 0.0f, -1.0f,
		-1.0f, 1.0f,-1.0f,   0.0f, 1.0f,   0.0f, 0.0f, -1.0f,
		1.0f, 1.0f,-1.0f,   1.0f, 1.0f,   0.0f, 0.0f, -1.0f,

		// left
		-1.0f,-1.0f, 1.0f,   0.0f, 1.0f,   -1.0f, 0.0f, 0.0f,
		-1.0f, 1.0f,-1.0f,   1.0f, 0.0f,   -1.0f, 0.0f, 0.0f,
		-1.0f,-1.0f,-1.0f,   0.0f, 0.0f,   -1.0f, 0.0f, 0.0f,
		-1.0f,-1.0f, 1.0f,   0.0f, 1.0f,   -1.0f, 0.0f, 0.0f,
		-1.0f, 1.0f, 1.0f,   1.0f, 1.0f,   -1.0f, 0.0f, 0.0f,
		-1.0f, 1.0f,-1.0f,   1.0f, 0.0f,   -1.0f, 0.0f, 0.0f,

		// right
		1.0f,-1.0f, 1.0f,   1.0f, 1.0f,   1.0f, 0.0f, 0.0f,
		1.0f,-1.0f,-1.0f,   1.0f, 0.0f,   1.0f, 0.0f, 0.0f,
		1.0f, 1.0f,-1.0f,   0.0f, 0.0f,   1.0f, 0.0f, 0.0f,
		1.0f,-1.0f, 1.0f,   1.0f, 1.0f,   1.0f, 0.0f, 0.0f,
		1.0f, 1.0f,-1.0f,   0.0f, 0.0f,   1.0f, 0.0f, 0.0f,
		1.0f, 1.0f, 1.0f,   0.0f, 1.0f,   1.0f, 0.0f, 0.0f
	};
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexData), vertexData, GL_STATIC_DRAW);
	// connect the xyz to the "vert" attribute of the vertex shader
	glEnableVertexAttribArray(asset.shaders->attrib("vert")); //Tutorial#6: Lighting
	glVertexAttribPointer(asset.shaders->attrib("vert"), 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), NULL);

	// connect the uv coords to the "vertTexCoord" attribute of the vertex shader
	glEnableVertexAttribArray(asset.shaders->attrib("vertTexCoord"));
	glVertexAttribPointer(asset.shaders->attrib("vertTexCoord"), 2, GL_FLOAT, GL_TRUE, 8 * sizeof(GLfloat), (const GLvoid*)(3 * sizeof(GLfloat)));

	// connect the normal to the "vertNormal" attribute of the vertex shader
	glEnableVertexAttribArray(asset.shaders->attrib("vertNormal"));
	glVertexAttribPointer(asset.shaders->attrib("vertNormal"), 3, GL_FLOAT, GL_TRUE, 8 * sizeof(GLfloat), (const GLvoid*)(5 * sizeof(GLfloat)));

	// unbind the VBO and VAO
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

static void LoadTriangleAsset(ModelAsset &asset) {
	asset.shaders = LoadShaders("vertex-shader-triangle.txt", "fragment-shader-triangle.txt"); //Tutorial#5
	asset.drawType = GL_TRIANGLES;
	asset.drawStart = 0;
	asset.drawCount = 3;
	asset.texture = LoadTexture("hazard.png");

	asset.shininess = 80.0; //Tutorial#7: More lighting
	asset.specularColor = glm::vec3(1.0f, 1.0f, 1.0f);

	glGenBuffers(1, &asset.vbo);
	glGenVertexArrays(1, &asset.vao);
	// bind the VAO
	glBindVertexArray(asset.vao);
	// bind the VBO
	glBindBuffer(GL_ARRAY_BUFFER, asset.vbo);

	GLfloat vertexData[] = {	//Tutorial#6: Lighting
		//  X     Y     Z       U     V          Normal
		0.0f, 0.8f, 0.0f,   0.5f, 1.0f,   0.0f, 0.0f, 1.0f,
		-0.8f,-0.8f,0.0f,   0.0f, 0.0f,   0.0f, 0.0f, 1.0f,
		0.8f,-0.8f, 0.0f,   1.0f, 0.0f,   0.0f, 0.0f, 1.0f,
	};
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexData), vertexData, GL_STATIC_DRAW);
	// connect the xyz to the "vert" attribute of the vertex shader
	glEnableVertexAttribArray(asset.shaders->attrib("vert")); //Tutorial#6: Lighting
	glVertexAttribPointer(asset.shaders->attrib("vert"), 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), NULL);

	// connect the uv coords to the "vertTexCoord" attribute of the vertex shader
	glEnableVertexAttribArray(asset.shaders->attrib("vertTexCoord"));
	glVertexAttribPointer(asset.shaders->attrib("vertTexCoord"), 2, GL_FLOAT, GL_TRUE, 8 * sizeof(GLfloat), (const GLvoid*)(3 * sizeof(GLfloat)));

	// connect the normal to the "vertNormal" attribute of the vertex shader
	glEnableVertexAttribArray(asset.shaders->attrib("vertNormal"));
	glVertexAttribPointer(asset.shaders->attrib("vertNormal"), 3, GL_FLOAT, GL_TRUE, 8 * sizeof(GLfloat), (const GLvoid*)(5 * sizeof(GLfloat)));

	// unbind the VBO and VAO
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

static void LoadFilledTriangleAsset(ModelAsset &asset) {
	asset.shaders = LoadShaders("vertex-shader-filledtriangle.txt", "fragment-shader-filledtriangle.txt"); //Tutorial#5
	asset.drawType = GL_TRIANGLES;
	asset.drawStart = 0;
	asset.drawCount = 3;
	//asset.texture = LoadTexture("hazard.png");
	asset.useTexture = false;

	//asset.shininess = 80.0; //Tutorial#7: More lighting
	//asset.specularColor = glm::vec3(1.0f, 1.0f, 1.0f);
	asset.useLight = false;

	glGenBuffers(1, &asset.vbo);
	glGenVertexArrays(1, &asset.vao);
	// bind the VAO
	glBindVertexArray(asset.vao);
	// bind the VBO
	glBindBuffer(GL_ARRAY_BUFFER, asset.vbo);

	GLfloat vertexData[] = {	//Tutorial#6: Lighting
		//  X     Y     Z       U     V          Normal
		0.0f, 0.8f, 0.0f,   0.5f, 1.0f,   0.0f, 0.0f, 1.0f,
		-0.8f,-0.8f,0.0f,   0.0f, 0.0f,   0.0f, 0.0f, 1.0f,
		0.8f,-0.8f, 0.0f,   1.0f, 0.0f,   0.0f, 0.0f, 1.0f,
	};
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexData), vertexData, GL_STATIC_DRAW);
	// connect the xyz to the "vert" attribute of the vertex shader
	glEnableVertexAttribArray(asset.shaders->attrib("vert")); //Tutorial#6: Lighting
	glVertexAttribPointer(asset.shaders->attrib("vert"), 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), NULL);

	// connect the uv coords to the "vertTexCoord" attribute of the vertex shader
	//glEnableVertexAttribArray(asset.shaders->attrib("vertTexCoord"));
	//glVertexAttribPointer(asset.shaders->attrib("vertTexCoord"), 2, GL_FLOAT, GL_TRUE, 8 * sizeof(GLfloat), (const GLvoid*)(3 * sizeof(GLfloat)));

	// connect the normal to the "vertNormal" attribute of the vertex shader
	//glEnableVertexAttribArray(asset.shaders->attrib("vertNormal"));
	//glVertexAttribPointer(asset.shaders->attrib("vertNormal"), 3, GL_FLOAT, GL_TRUE, 8 * sizeof(GLfloat), (const GLvoid*)(5 * sizeof(GLfloat)));

	// unbind the VBO and VAO
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}


static void LoadFilledRectAsset(ModelAsset &asset) {
	asset.shaders = LoadShaders("vertex-shader-filledrect.txt", "fragment-shader-filledrect.txt"); //Tutorial#5
	asset.drawType = GL_TRIANGLES;
	asset.drawStart = 0;
	asset.drawCount = 6;

	asset.useTexture = false;
	asset.useLight = false;

	glGenBuffers(1, &asset.vbo);
	glGenVertexArrays(1, &asset.vao);
	glBindVertexArray(asset.vao);
	glBindBuffer(GL_ARRAY_BUFFER, asset.vbo);

	//float w = 0.016f / 2;
	//float h = 0.016f / 2;
	float w = 1.0f / 2;
	float h = 1.0f / 2;

	GLfloat vertexData[] = {	//Tutorial#6: Lighting
		//  X     Y     Z      
		-w, h, 0.0f,
		-w,-h, 0.0f,
		w,-h, 0.0f,
		-w, h, 0.0f,
		w,-h, 0.0f,
		w, h, 0.0f,
	};
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexData), vertexData, GL_STATIC_DRAW);

	glEnableVertexAttribArray(asset.shaders->attrib("vert")); //Tutorial#6: Lighting
	glVertexAttribPointer(asset.shaders->attrib("vert"), 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), NULL);


	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

static void LoadLineAsset(ModelAsset &asset) {
	asset.shaders = LoadShaders("vertex-shader-line.txt", "fragment-shader-line.txt"); //Tutorial#5
	asset.drawType = GL_LINES;
	asset.drawStart = 0;
	asset.drawCount = 2;
	//asset.texture = LoadTexture("wooden-crate.jpg");
	asset.useTexture = false;

	//asset.shininess = 80.0; //Tutorial#7: More lighting
	//asset.specularColor = glm::vec3(1.0f, 1.0f, 1.0f);
	asset.useLight = false;

	glGenBuffers(1, &asset.vbo);
	glGenVertexArrays(1, &asset.vao);
	// bind the VAO
	glBindVertexArray(asset.vao);
	// bind the VBO
	glBindBuffer(GL_ARRAY_BUFFER, asset.vbo);

	GLfloat vertexData[] = {	//Tutorial#6: Lighting
		//  X     Y     Z
		0.0f, 0.0f, 0.001f,   
		1.0f, 0.0f, 0.001f,   
	};
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexData), vertexData, GL_STATIC_DRAW);
	// connect the xyz to the "vert" attribute of the vertex shader
	glEnableVertexAttribArray(asset.shaders->attrib("vert")); //Tutorial#6: Lighting
	glVertexAttribPointer(asset.shaders->attrib("vert"), 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), NULL);

	// unbind the VBO and VAO
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

static void LoadCircleAsset(ModelAsset &asset) {
	int size = 360;

	asset.shaders = LoadShaders("vertex-shader-circle.txt", "fragment-shader-circle.txt"); //Tutorial#5
	asset.drawType = GL_LINE_LOOP;
	asset.drawStart = 0;
	asset.drawCount = size;
	//asset.texture = LoadTexture("wooden-crate.jpg");
	asset.useTexture = false;

	//asset.shininess = 80.0; //Tutorial#7: More lighting
	//asset.specularColor = glm::vec3(1.0f, 1.0f, 1.0f);
	asset.useLight = false;

	glGenBuffers(1, &asset.vbo);
	glGenVertexArrays(1, &asset.vao);
	// bind the VAO
	glBindVertexArray(asset.vao);
	// bind the VBO
	glBindBuffer(GL_ARRAY_BUFFER, asset.vbo);

	

	int sizeOfVertexData = size * 3;// +size * 4;
	GLfloat *vertexData = new GLfloat[sizeOfVertexData];

	float radius = 1.0f / 2;
	const float DEG2RAD = PI * 2 / size;
	for (int i = 0; i < size; i++)
	{
		float degInRad = i*DEG2RAD;
		vertexData[i * 3 + 0] = cos(degInRad)*radius;
		vertexData[i * 3 + 1] = sin(degInRad)*radius;
		vertexData[i * 3 + 2] = 0;

		//vertexData[i * 7 + 3] = 1;
		//vertexData[i * 7 + 4] = 1;
		//vertexData[i * 7 + 5] = 1;
		//vertexData[i * 7 + 6] = 1;
	}

	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * sizeOfVertexData, vertexData, GL_STATIC_DRAW);
	// connect the xyz to the "vert" attribute of the vertex shader
	glEnableVertexAttribArray(asset.shaders->attrib("vert")); //Tutorial#6: Lighting
	glVertexAttribPointer(asset.shaders->attrib("vert"), 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), NULL);

	//glEnableVertexAttribArray(asset.shaders->attrib("vertexColor"));
	//glVertexAttribPointer(asset.shaders->attrib("vertexColor"), 4, GL_FLOAT, GL_TRUE, 7 * sizeof(GLfloat), (const GLvoid*)(3 * sizeof(GLfloat)));


	// unbind the VBO and VAO
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	delete vertexData;
}

/*
static void LoadFieldSystemAsset(ModelAsset &asset, CFieldSystem *fieldSystem) { //not used
	asset.shaders = LoadShaders("vertex-shader-line.txt", "fragment-shader-line.txt"); //Tutorial#5
	asset.drawType = GL_LINES;
	asset.drawStart = 0;
	asset.drawCount = 2 * fieldSystem->GetHeight() * fieldSystem->GetWidth();
	asset.useTexture = false;
	asset.useLight = false;

	glGenBuffers(1, &asset.vbo);
	glGenVertexArrays(1, &asset.vao);
	glBindVertexArray(asset.vao);// bind the VAO
	glBindBuffer(GL_ARRAY_BUFFER, asset.vbo);// bind the VBO


	int sizeOfVertexData = 3 * 2 * fieldSystem->GetHeight() * fieldSystem->GetWidth();
	GLfloat *vertexData = new GLfloat[sizeOfVertexData];
	int idx = 0;
	for (int i = 0; i < fieldSystem->GetHeight(); i++) {
		for (int j = 0; j < fieldSystem->GetWidth(); j++) {
			vertexData[idx * 6 + 0] = fieldSystem->grid->location[i][j].x;
			vertexData[idx * 6 + 1] = fieldSystem->grid->location[i][j].y;
			vertexData[idx * 6 + 2] = fieldSystem->grid->location[i][j].z;

			vertexData[idx * 6 + 3] = fieldSystem->grid->location[i][j].x + fieldSystem->nodes->velocity[i][j].x;
			vertexData[idx * 6 + 4] = fieldSystem->grid->location[i][j].y + fieldSystem->nodes->velocity[i][j].y;
			vertexData[idx * 6 + 5] = fieldSystem->grid->location[i][j].z + fieldSystem->nodes->velocity[i][j].z;

			idx++;
		}
	}

	glBufferData(GL_ARRAY_BUFFER, sizeOfVertexData * sizeof(GLfloat), vertexData, GL_STATIC_DRAW);
	glEnableVertexAttribArray(asset.shaders->attrib("vert")); 
	glVertexAttribPointer(asset.shaders->attrib("vert"), 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), NULL);

	// unbind the VBO and VAO
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	delete vertexData;
}*/