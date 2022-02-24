// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OPENGL_GLFW_VIEWER_H
#define IGL_OPENGL_GLFW_VIEWER_H

#ifndef IGL_OPENGL_4
#define IGL_OPENGL_4
#endif

#include "../../igl_inline.h"
#include "../MeshGL.h"

#include "../ViewerData.h"
#include "ViewerPlugin.h"
#include "../ViewerCore.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vector>
#include <string>
#include <cstdint>
#include <igl\AABB.h>

#define IGL_MOD_SHIFT           0x0001
#define IGL_MOD_CONTROL         0x0002
#define IGL_MOD_ALT             0x0004
#define IGL_MOD_SUPER           0x0008


struct GLFWwindow;


namespace igl
{
namespace opengl
{
namespace glfw
{
  // GLFW-based mesh viewer
  class Viewer : public Movable
  {
  public:
    // UI Enumerations
   // enum class MouseButton {Left, Middle, Right};
   // enum class MouseMode { None, Rotation, Zoom, Pan, Translation} mouse_mode;


      // base for assignment 

      IGL_INLINE void Viewer::assignmentsChoosing(int assignment);


// ASSIGNMENT 1 
      bool simplification_enable ;
      IGL_INLINE bool load_mesh_from_configuration(const std::string config , bool assignment2);
      IGL_INLINE void Viewer::init_objs_simpelified();
      IGL_INLINE void simplify_mesh(int num_of_faces_to_delete);
      IGL_INLINE bool collapse_edge(
          const std::function<void(
              std::vector<Eigen::Matrix<double, 4, 4>>& Qs,
              const int,
              const Eigen::MatrixXd&,
              const Eigen::MatrixXi&,
              const Eigen::MatrixXi&,
              const Eigen::VectorXi&,
              const Eigen::MatrixXi&,
              const Eigen::MatrixXi&,
              double&,
              Eigen::RowVectorXd&)>& cost_and_placement,
          std::vector<Eigen::Matrix<double, 4, 4>>& Qs,
          Eigen::MatrixXd& V,
          Eigen::MatrixXi& F,
          Eigen::MatrixXi& E,
          Eigen::VectorXi& EMAP,
          Eigen::MatrixXi& EF,
          Eigen::MatrixXi& EI,
          std::set<std::pair<double, int> >& Q,
          std::vector<std::set<std::pair<double, int> >::iterator >& Qit,
          Eigen::MatrixXd& C);
      IGL_INLINE bool Viewer::check_if_infinity(double first);


      IGL_INLINE bool Viewer::collapse(
          const int e,
          const Eigen::RowVectorXd& p,
          Eigen::MatrixXd& V,
          Eigen::MatrixXi& F,
          Eigen::MatrixXi& E,
          Eigen::VectorXi& EMAP,
          Eigen::MatrixXi& EF,
          Eigen::MatrixXi& EI,
          std::vector<Eigen::Matrix<double, 4, 4>>& Qs,
          int& a_e1,
          int& a_e2,
          int& a_f1,
          int& a_f2);


   
   
   // ASSIGNMENT 2 
      std::vector< igl::AABB<Eigen::MatrixXd, 3>*> tree_roots;
      std::vector< igl::AABB<Eigen::MatrixXd, 3>*> trees;
      std::vector<Eigen::AlignedBox<double, 3>*> last_box;
      IGL_INLINE void Viewer::draw_m_box(int index, Eigen::AlignedBox<double, 3>& m_box, Eigen::RowVector3d color);
      IGL_INLINE void Viewer::start_collision_render();
      IGL_INLINE bool Viewer::check_collision_occurence(igl::AABB<Eigen::MatrixXd, 3>* first, igl::AABB<Eigen::MatrixXd, 3>* second, Eigen::Vector3d A[], Eigen::Vector3d B[]);
      IGL_INLINE bool Viewer::check_possible_seperate(Eigen::AlignedBox<double, 3>& first, Eigen::AlignedBox<double, 3>& second, Eigen::Vector3d A[], Eigen::Vector3d B[]);
      IGL_INLINE void Viewer::stop_collision(void);

      IGL_INLINE void Viewer::init_data();

    bool render_collison = false;


    //assingment 3 

    bool animateIk = false; 
    //IGL_INLINE void Viewer::start_animate_ik();
    IGL_INLINE void Viewer::init_ik_mesh();
    double doubleVariable;

    IGL_INLINE Eigen::Matrix4f Viewer::get_parent_link_T(int index);
    IGL_INLINE Eigen::Matrix4f Viewer::get_parent_link_Rot(int index);

    int sphere_index;
    int link_number = 4;
    int root_link_index;
    int last_link_index;




    void Viewer::initAxes();
    void Viewer::initLinkAxes();
    void Viewer::updateTipPos();
    void Viewer::updateDestPos();
    void Viewer::CCD();
    void Viewer::Fabrik();
    void Viewer::updateLinksToTips(std::vector<Eigen::Vector4d> newPos);
    std::vector<Eigen::Vector4d> tipPos;
    std::vector<Eigen::Vector4f> tipPos1;
    int linksNum;
    void Viewer::fix_axis_rotation();
    void Viewer::printRotation();
    void Viewer::printTip();
    Eigen::Vector4d destPos;
 

    virtual void Init(const std::string config);
	virtual void Animate() {}
	virtual void WhenTranslate() {}
	virtual Eigen::Vector3d GetCameraPosition() { return Eigen::Vector3d(0, 0, 0); }
	virtual Eigen::Vector3d GetCameraForward() { return Eigen::Vector3d(0, 0, -1); }
	virtual Eigen::Vector3d GetCameraUp() { return Eigen::Vector3d(0, 1, 0); }
    IGL_INLINE Eigen::Matrix4f Viewer::get_parent_link_T();
    IGL_INLINE void Viewer::init_plugins();
    bool picked;
    
     GLFWwindow* window;

    // final project 
     bool Viewer::check_index_collision(int i);
     std::string configFileLevel;
    
     bool level1 ;
     bool level2 ;
     typedef
         std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> >
         RotationList;
     Eigen::MatrixXd V, W, C, U, M;
     Eigen::Vector3d destination_position;
     std::vector<Eigen::Vector3d>chain;
     int scale;
     int Num_Of_Joints;
     std::vector<Eigen::Vector3d>skelton;
     std::vector<Movable> Joints;
     std::vector<Eigen::Vector3d> vT;
     std::vector<int> parentsJoints;
     RotationList vQ;
     bool up;
     bool down;
     bool right;
     bool left;
     bool rotDir;
     bool skinning;

     bool levelWindow;
     bool finishLevel;
     int score;
     int level;

     bool snakeEye;


     void CalcNextPosition();

     bool Is_Collide(Eigen::AlignedBox<double, 3>& A_box, int indexA, Eigen::AlignedBox<double, 3>& B_box, int indexB);

     bool Check_Collision(igl::AABB<Eigen::MatrixXd, 3>& A, int indexA, igl::AABB<Eigen::MatrixXd, 3>& B, int indexB);

     bool Check_Collision();

     IGL_INLINE bool check_snake_collision();

     void Calculate_Weights();

     Eigen::VectorXd creatWiVector(Eigen::Vector4d temp);

     void Initialize_Tree(int index);

     IGL_INLINE void Viewer::init_final_project();
    IGL_INLINE void init_snake_mesh();
    bool running_game = false;
    std::string* data_path;

    bool pause_game = false;
    bool finish_game = false;
    bool end_game = false;
    void* game;
    
    std::string cyl = "/ycylinder.obj";
    std::string sph = "/sphere.obj";
    int frames = 0;
    int left_view;
    int right_view;
    igl::opengl::ViewerCore* second_camera;

  



    std::vector<ViewerPlugin*> plugins;
    IGL_INLINE void setWindow(GLFWwindow* window);
    IGL_INLINE void setRenderer(void* render);
    IGL_INLINE igl::opengl::ViewerCore& core(unsigned core_id = 0);
    IGL_INLINE const igl::opengl::ViewerCore& core(unsigned core_id = 0) const;
    IGL_INLINE void snap_to_canonical_quaternion();
 
    IGL_INLINE void shutdown_plugins();

    void Viewer::updateDestPos1();






	//IGL_INLINE void init_plugins();
    //IGL_INLINE void shutdown_plugins();
    Viewer();
    virtual ~Viewer();
    // Mesh IO
    IGL_INLINE bool load_mesh_from_file(const std::string & mesh_file_name);
    IGL_INLINE bool save_mesh_to_file(const std::string & mesh_file_name);
   
    // Scene IO
    IGL_INLINE bool load_scene();
    IGL_INLINE bool load_scene(std::string fname);
    IGL_INLINE bool save_scene();
    IGL_INLINE bool save_scene(std::string fname);
    // Draw everything
   // IGL_INLINE void draw();
    // OpenGL context resize
   
    // Helper functions

    IGL_INLINE void open_dialog_load_mesh();
    IGL_INLINE void open_dialog_save_mesh();

	IGL_INLINE void draw() {}
    ////////////////////////
    // Multi-mesh methods //
    ////////////////////////

    // Return the current mesh, or the mesh corresponding to a given unique identifier
    //
    // Inputs:
    //   mesh_id  unique identifier associated to the desired mesh (current mesh if -1)
    IGL_INLINE ViewerData& data(int mesh_id = -1);
    IGL_INLINE const ViewerData& data(int mesh_id = -1) const;

    // Append a new "slot" for a mesh (i.e., create empty entries at the end of
    // the data_list and opengl_state_list.
    //
    // Inputs:
    //   visible  If true, the new mesh is set to be visible on all existing viewports
    // Returns the id of the last appended mesh
    //
    // Side Effects:
    //   selected_data_index is set this newly created, last entry (i.e.,
    //   #meshes-1)
    IGL_INLINE int append_mesh(bool visible = true);

    // Erase a mesh (i.e., its corresponding data and state entires in data_list
    // and opengl_state_list)
    //
    // Inputs:
    //   index  index of mesh to erase
    // Returns whether erasure was successful <=> cannot erase last mesh
    //
    // Side Effects:
    //   If selected_data_index is greater than or equal to index then it is
    //   decremented
    // Example:
    //   // Erase all mesh slots except first and clear remaining mesh
    //   viewer.selected_data_index = viewer.data_list.size()-1;
    //   while(viewer.erase_mesh(viewer.selected_data_index)){};
    //   viewer.data().clear();
    //
    IGL_INLINE bool erase_mesh(const size_t index);

    // Retrieve mesh index from its unique identifier
    // Returns 0 if not found
    IGL_INLINE size_t mesh_index(const int id) const;

	Eigen::Matrix4d CalcParentsTrans(int indx);
	inline bool SetAnimation() { return isActive = !isActive; }
public:
    //////////////////////
    // Member variables //
    //////////////////////

    // Alec: I call this data_list instead of just data to avoid confusion with
    // old "data" variable.
    // Stores all the data that should be visualized
    std::vector<ViewerData> data_list;
	
	std::vector<int> parents;

    size_t selected_data_index;
    int next_data_id;
	bool isPicked;
	bool isActive;


    

    // List of registered plugins
//    std::vector<ViewerPlugin*> plugins;

    // Keep track of the global position of the scrollwheel
    float scroll_position;

  public:
    
};

} // end namespace
} // end namespace
} // end namespace

#ifndef IGL_STATIC_LIBRARY
#  include "Viewer.cpp"
#endif

#endif
