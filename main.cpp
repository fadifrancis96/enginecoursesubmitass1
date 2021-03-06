
#include "igl/opengl/glfw/renderer.h"
#include "tutorial/sandBox/inputManager.h"
#include "sandBox.h"

int main(int argc, char* argv[])
{
	Display* disp = new Display(1200, 800, "Wellcome");
	Renderer renderer;

	SandBox viewer;

	igl::opengl::glfw::imgui::ImGuiMenu* menu = new igl::opengl::glfw::imgui::ImGuiMenu();


	// ASSIGNMENT 1 - TASK 4 

	viewer.Init("configuration.txt");
	printf("done with init \n");

	Init(*disp, menu);

	renderer.init(&viewer, 2, menu);

	disp->SetRenderer(&renderer);
	disp->launch_rendering(true);
	delete menu;
	delete disp;
}