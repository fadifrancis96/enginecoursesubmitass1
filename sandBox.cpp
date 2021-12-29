#include "tutorial/sandBox/sandBox.h"
#include "igl/edge_flaps.h"
#include "igl/collapse_edge.h"
#include "Eigen/dense"
#include <functional>



SandBox::SandBox()
{
	

}

void SandBox::Init(const std::string &config)
{
	bool assignment2 = false ; 
	bool assignment3 = false;

	// choose please what assignment is on  
	if (assignment3) {

		Viewer::init_ik_mesh();
		printf("fadi \n");
		
	}

	else {


		float mesh = 0;
		std::string mesh_path;
		std::string item_name;
		std::ifstream nameFileout;
		doubleVariable = 0;
		nameFileout.open(config);

		if (!nameFileout.is_open())
		{
			std::cout << "Can't open file " << config << std::endl;
		}
		else
		{

			while (getline(nameFileout, mesh_path))
			{

				std::cout << "im here trying to open  \n" << item_name << std::endl;

				this->load_mesh_from_file(mesh_path);

				//TODO


				double meshtmp = mesh - 0.5;
				double meshtmp1 = meshtmp * 0.5;

				Viewer::data().MyTranslate(Eigen::Vector3d(meshtmp1, 0, 0), false);

				igl::AABB<Eigen::MatrixXd, 3>* tree = new igl::AABB<Eigen::MatrixXd, 3>();

				tree->init(Viewer::data().V, Viewer::data().F);

				Viewer::trees.push_back(tree);
				Viewer::tree_roots.push_back(tree);
				Viewer::last_box.push_back(&tree->m_box);

				data().original_V = new Eigen::MatrixXd(this->data().V);
				data().original_F = new Eigen::MatrixXi(this->data().F);

				Viewer::data().show_overlay_depth = false;
				Viewer::data().show_lines = false;
				Viewer::data().point_size = 1;
				Viewer::data().line_width = 1;

				Eigen::RowVector3d red_color(1, 0, 0);

				// this is the function that will draw box around our mesh 
				Viewer::draw_m_box(this->selected_data_index, tree->m_box, red_color);

				mesh++;


				// here should 2 objects already be loaded 

			//	Viewer::collission_init(item_name);
				//getline(infile, mesh_path);
			//	this->load_mesh_from_file(item_name);
				//Viewer::load_mesh_from_configuration(item_name , true );


				//Viewer::init_objs_simpelified();





				parents.push_back(-1);
				data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(0, 0, 1));
				//data().show_overlay_depth = false;
				//data().point_size = 10;
				//data().line_width = 2;
				//data().set_visible(false, 1);


			}
			nameFileout.close();
		}

		MyTranslate(Eigen::Vector3d(0, 0, -1), true);

		data().set_colors(Eigen::RowVector3d(0.9, 0.1, 0.1));
	}
}

SandBox::~SandBox()
{

}

void SandBox::Animate()
{
	if (isActive)
	{
		
		
		
	}
}


