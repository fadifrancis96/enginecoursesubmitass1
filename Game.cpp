#include "Game.h"
#include <iostream>
#include <cstdlib>
#include "igl/dqs.h"

Game::Game(igl::opengl::glfw::Viewer* viewer) {

	this->viewer = viewer;
	this->snake = new Snake(17, viewer);

	



}
void Game::init() {
	viewer->right = false;
	viewer->left = false;
	viewer->up = false;
	viewer->down = false;
	viewer->rotDir = false;
	viewer->snakeEye = 0;
	viewer->level = 1;
	viewer->score = 0;
	viewer->finishLevel = false;
	viewer->levelWindow = false;
	viewer->skinning = false;
    viewer->Num_Of_Joints = 16;
	viewer->skelton.resize(viewer->Num_Of_Joints + 1);
	viewer->chain.resize(viewer->Num_Of_Joints + 1);
	viewer->parentsJoints.resize(viewer->Num_Of_Joints + 1);
	viewer->scale = 1;

	//Initialize vT, vQ
	viewer->vT.resize(17);
	viewer->vQ.resize(17);

	this->level1 = true;	
	this->level2 = false;
	viewer->level1 = true;
	viewer->level2 = false;
	this->enemies = 6;
	this->current_enemies = 4;
	this->level = 0;
	this->score = 0;
	this->isOn = false;
	this->max_level = 5;
	snake->init("configuration.txt");
	this->level1 = true;
	this->level2 = false;

    viewer->picked = true;
}







Game::~Game() {
	std::cout << "ending game..." << std::endl;
}



void Game::loop() {

	if (viewer->isActive)
	{
		

		//Option 4 - Only Collision Detection
	
			if (viewer->left)
				viewer->data().MyTranslate(Eigen::Vector3d(-0.03, 0, 0), true);
			if (viewer->right)
				viewer->data().MyTranslate(Eigen::Vector3d(0.03, 0, 0), true);
			if (viewer->up)
				viewer->data().MyTranslate(Eigen::Vector3d(0, 0.03, 0), true);
			if (viewer->down)
				viewer->data().MyTranslate(Eigen::Vector3d(0, -0.03, 0), true);



			
				if (viewer->check_snake_collision())
				{
					printf("im in level 1\n");
					viewer->score += 10;
					if (viewer->score == 50) {
						viewer->finishLevel = true;
						viewer->left = false;
						viewer->right = false;
						viewer->up = false;
						viewer->down = false;
						viewer->isActive = !viewer->isActive;
						
					}
					// now what happens when i move up a notch ? and there is another LEVEL ? 
					   // i want to add one more OBJECT FOR THE SNAKE 


					int y = rand() % 60 - 30;
					double dy = y / 10;
					int x = rand() % 60 - 30;
					double dx = x / 10;
					viewer->data_list[0].MyTranslate(Eigen::Vector3d(dx, dy, 0), true);
					if (viewer->level2) {
						int y = rand() % 60 - 30;
						double dy = y / 10;
						int x = rand() % 60 - 30;
						double dx = x / 10;
						viewer->data_list[2].MyTranslate(Eigen::Vector3d(dx, dy, 0), true);
					}

				
			
			
		}






	}








}





