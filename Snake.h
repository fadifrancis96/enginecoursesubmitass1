#pragma once
#include <igl\opengl\glfw\Viewer.h>

class Snake
{
public:
	Snake(int length, igl::opengl::glfw::Viewer* viewer);
	~Snake();

	void init(const std::string& config);



private:
	igl::opengl::glfw::Viewer* viewer;
};

