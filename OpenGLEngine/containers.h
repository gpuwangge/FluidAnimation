#pragma once

/*
Represents a textured geometry asset
Contains everything necessary to draw arbitrary geometry with a single texture:
- shaders
- a texture
- a VBO
- a VAO
- the parameters to glDrawArrays (drawType, drawStart, drawCount)
*/
struct ModelAsset {
	tdogl::Program* shaders; //shader is tdogl::Program
	tdogl::Texture* texture;
	GLuint vbo;
	GLuint vao;
	GLenum drawType;
	GLint drawStart;
	GLint drawCount;

	GLfloat shininess; //Tutorial#7: More lighting
	glm::vec3 specularColor; 

	bool useTexture;
	bool useLight;

	ModelAsset() :
		shaders(NULL),
		texture(NULL),
		vbo(0),
		vao(0),
		drawType(GL_TRIANGLES),
		//drawType(GL_LINE),
		drawStart(0),
		drawCount(0),
		shininess(0.0f),
		specularColor(1.0f, 1.0f, 1.0f),
		useTexture(true),
		useLight(true)
	{}
};



struct ModelInstance {
	ModelAsset* m_asset;
	glm::mat4 m_transform;
	glm::vec4 m_color;

	ModelInstance() :
		m_asset(NULL),
		m_transform(),
		m_color()
	{}

	ModelInstance(ModelAsset* asset, glm::mat4 transform, glm::vec4 color) :
		m_asset(asset),
		m_transform(transform),
		m_color(color)
	{}
};


/*
Represents a point light
*/
struct Light {
	glm::vec4 position; //Tutorial#8: Even more lighting
	glm::vec3 intensities; //a.k.a. the color of the light

	float attenuation; //Tutorial#7: More lighting
	float ambientCoefficient; 

	float coneAngle; //Tutorial#8: Even more lighting
	glm::vec3 coneDirection;
};
