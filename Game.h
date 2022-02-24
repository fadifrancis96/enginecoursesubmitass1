#pragma once
#include <tutorial\sandBox\Snake.h>



class Game
{
public:
	Game(igl::opengl::glfw::Viewer* viewer);
	~Game();
    void init();
	void loop();
	int high_score;
	int score;
	int level;
	int max_level ;
	int enemies;
	int current_enemies;
	bool running = true;
	bool isOn;
	bool gotOnce = false;
	std::string user_name = "";
	std::string old_user_name = "";

private:
	
	Snake* snake;
	igl::opengl::glfw::Viewer* viewer;
	bool level1=false;
	bool level2=true;
	
};

