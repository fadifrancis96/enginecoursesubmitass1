#include "Snake.h"
#include <iostream>
#include <igl\opengl\ViewerCore.h>




Snake::Snake(int length, igl::opengl::glfw::Viewer* viewer) {

	this->viewer = viewer;
}

Snake::~Snake() {
	std::cout << "i'm dying :'(" << std::endl;
}
void Snake::init(const std::string& config) {

	float mesh = 0;
	std::string mesh_path;
	std::string item_name;
	std::ifstream nameFileout;
	viewer->doubleVariable = 0;
	nameFileout.open(config);

	if (!nameFileout.is_open())
	{
		std::cout << "Can't open file " <<"configuration.txt" << std::endl;
	}
	else
	{
	
		while (getline(nameFileout, mesh_path))
		{

			std::cout << "im here trying to open  \n" << item_name << std::endl;

			viewer->load_mesh_from_file(mesh_path);

			//TODO


			double meshtmp = mesh - 0.5;
			double meshtmp1 = meshtmp * 0.5;

			viewer->data().MyTranslate(Eigen::Vector3d(meshtmp1, 0, 0), false);

			igl::AABB<Eigen::MatrixXd, 3>* tree = new igl::AABB<Eigen::MatrixXd, 3>();

			tree->init(viewer->data().V, viewer->data().F);

			viewer->trees.push_back(tree);
			viewer->tree_roots.push_back(tree);
			viewer->last_box.push_back(&tree->m_box);

			viewer->data().original_V = new Eigen::MatrixXd(viewer->data().V);
			viewer->data().original_F = new Eigen::MatrixXi(viewer->data().F);

			viewer->data().show_overlay_depth = false;
			viewer->data().show_lines = false;
			viewer->data().point_size = 1;
			viewer->data().line_width = 1;

			Eigen::RowVector3d red_color(1, 0, 0);

			// this is the function that will draw box around our mesh 
			viewer->draw_m_box(viewer->selected_data_index, tree->m_box, red_color);

			mesh++;

			viewer->parents.push_back(-1);
			viewer->data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(0, 0, 1));
	


		}
		nameFileout.close();
	}


	viewer->data().set_colors(Eigen::RowVector3d(0.9, 0.1, 0.1));
}









