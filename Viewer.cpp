// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "Viewer.h"

//#include <chrono>
#include <thread>

#include <Eigen/LU>


#include <cmath>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <cassert>

#include <igl/project.h>
//#include <igl/get_seconds.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/adjacency_list.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/massmatrix.h>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/quat_mult.h>
#include <igl/axis_angle_to_quat.h>
#include <igl/trackball.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/unproject.h>
#include <igl\circulation.h>
#include <igl/serialize.h>
#include<igl\edge_collapse_is_valid.h>
#include <igl\edge_flaps.h>
#include <igl\vertex_triangle_adjacency.h>
#include <igl/collapse_edge.h>



// Internal global variables used for glfw event handling
//static igl::opengl::glfw::Viewer * __viewer;
static double highdpi = 1;
static double scroll_x = 0;
static double scroll_y = 0;



namespace igl
{
namespace opengl
{
namespace glfw
{

  void Viewer::Init(const std::string config)
  {
      
      load_mesh_from_configuration(config);

  }

  IGL_INLINE Viewer::Viewer():
    data_list(1),
    selected_data_index(0),
    next_data_id(1),
	isPicked(false),
	isActive(false)
  {
    data_list.front().id = 0;

  

    // Temporary variables initialization
   // down = false;
  //  hack_never_moved = true;
    scroll_position = 0.0f;

    // Per face
    data().set_face_based(false);

    
#ifndef IGL_VIEWER_VIEWER_QUIET
    const std::string usage(R"(igl::opengl::glfw::Viewer usage:
  [drag]  Rotate scene
  A,a     Toggle animation (tight draw loop)
  F,f     Toggle face based
  I,i     Toggle invert normals
  L,l     Toggle wireframe
  O,o     Toggle orthographic/perspective projection
  T,t     Toggle filled faces
  [,]     Toggle between cameras
  1,2     Toggle between models
  ;       Toggle vertex labels
  :       Toggle face labels
  P    Enable simplefication
  x    Disable simplefication )"
);
    std::cout<<usage<<std::endl;
#endif
  }

  IGL_INLINE Viewer::~Viewer()
  {
  }
  void calculate_cost(
      std::vector<Eigen::Matrix<double, 4, 4>>& Qs,
      const int e,
      const Eigen::MatrixXd& V,
      const Eigen::MatrixXi& /*F*/,
      const Eigen::MatrixXi& E,
      const Eigen::VectorXi& /*EMAP*/,
      const Eigen::MatrixXi& /*EF*/,
      const Eigen::MatrixXi& /*EI*/,
      double& cost,
      Eigen::RowVectorXd& p)
  {
      auto Q1 = Qs[E.row(e)[0]];
      auto Q2 = Qs[E.row(e)[1]];
      Eigen::Matrix4d Q = Q1 + Q2;
      Q.row(3) = Eigen::Vector4d(0, 0, 0, 1);
      bool invertable = false;
      Eigen::Matrix<double, 4, 4> Qi;
      Q.computeInverseWithCheck(Qi, invertable, 0);
      Eigen::Matrix<double, 4, 1> a;
      if (invertable) {
          a = Qi * Eigen::Matrix<double, 4, 1>(0, 0, 0, 1);
      }
      else {
          a = (V.row(E.row(e)[0]) + V.row(E.row(e)[1])) / 2;
          a(3) = 1;
      }
      cost = a.transpose() * (Q1 + Q2) * a;
      p = Eigen::RowVector3d(a(0), a(1), a(2));
  }

  //Assingment 1 task 6
  IGL_INLINE void Viewer::init_objs_simpelified() {
     
      data().Q = new std::set<std::pair<double, int> >(); // priority Q contains the cost of every edge 
      data().EMAP = new Eigen::VectorXi();               // connects faces to edges 
      data().E = new Eigen::MatrixXi();                   // this is the Edges -> edge is represented by <index of source vertex, index of destinations vertex>
      data().EI = new Eigen::MatrixXi();                  //connects edge to vertex index in triangle (0 1 2)   
      data().EF = new Eigen::MatrixXi();                    //connects edges to faces
      data().C = new Eigen::MatrixXd();                 //position of the new vertex after collapsing the corresponding edge
      data().V1 = new Eigen::MatrixXd(data().V);
      data().F1 = new Eigen::MatrixXi(data().F);
      data().Qs = std::vector<Eigen::Matrix<double, 4, 4>>();
      data().Qit = new std::vector<std::set<std::pair<double, int> >::iterator >();
   
      edge_flaps(*data().F1, *data().E, *data().EMAP, *data().EF, *data().EI);

      data().Qit->resize(data().E->rows());
     
      data().C->resize(data().E->rows(), data().V.cols());
   
      Eigen::VectorXd costs(data().E->rows());
  
      data().Q->clear();
     


      auto VF = std::vector<std::vector<int> >(); // VF(v) is a vector of the face ids (index in data().F)
      auto VFi = std::vector<std::vector<int> >();// VFi(v) is a vector of the index of v in the vector that represents f
      vertex_triangle_adjacency(data().V, data().F, VF, VFi);
   

      for (int v = 0; v < data().V.rows(); v++) {

       
          std::vector<int> faces = VF[v];
          Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
          
          for (int face : faces) {
              auto n = data().F_normals.row(face).normalized();
              float d = 0;
              for (int j = 0; j < 3; j++) d += (-n[j] * data().V.row(v)[j]);
              Eigen::RowVector4d p(n[0], n[1], n[2], d);
              p = p.transpose();
              Q += p.transpose() * p;
          }
         
          data().Qs.push_back(Q);
      }
      for (int e = 0; e < data().E->rows(); e++)
      {
     // here i want to calculate all the edges cost
          double cost = 0;
          Eigen::RowVectorXd p(1, 3);

          calculate_cost(data().Qs, e, *data().V1, *data().F1, *data().E, *data().EMAP, *data().EF, *data().EI, cost, p);
          data().C->row(e) = p;
          (*data().Qit)[e] = (data().Q->insert(std::pair<double, int>(cost, e))).first;
      }

     
  }
      

  

  // ASSIGNMENT 1 TASK 7 & TASK 9 

  IGL_INLINE void Viewer::simplify_mesh(int num_of_faces_to_delete) {

      printf("# of faces to delete %d \n", num_of_faces_to_delete);
      bool something_collapsed = false;
      int num_of_collapsed_edges = 0; 




      
      for (int i = 0; i < num_of_faces_to_delete; i++) {
         
         
          if (!collapse_edge(calculate_cost, data().Qs,
              *data().V1,*data().F1,
              *data().E, *data().EMAP,
              *data().EF, *data().EI ,
              *data().Q, *data().Qit,
              * data().C))
          {
              break;

        
 

          }
          
                
          something_collapsed = true;
          data().num_collapsed++;

      }
      if (something_collapsed) {
          printf("succesfully collapsed # %d of edges \n",data().num_collapsed);
          data().clear();
          data().set_mesh(*data().V1, *data().F1);
          data().set_face_based(true);
          data().dirty = 157;



      }

      
				

  }
  IGL_INLINE bool Viewer::collapse_edge(
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
      Eigen::MatrixXd& C)
   {

      printf("in collapse edge func \n");

      using namespace Eigen;

      int e; 
      int e1; 
      int e2; 
      int f1; 
      int f2; 


      if (Q.empty()) {
          printf("Q is empty");
          return false;
      }
      printf("QUEUE SIZE IS %d", Q.size());

      // now we know that the qeueu is not empty and we have some data in it
      // we want to get the first pair , meaning the first edge - with the lowest cost? 

      std::pair<double, int> p = *(Q.begin());

      if (p.first == std::numeric_limits<double>::infinity())
      {
          printf("its infinity - returning \n");
          // min cost edge is infinite cost
          return false;
      }

      // the first edge is not ininity and we want to remove it 
      // remember that the data of first edge is saved in p
      Q.erase(Q.begin());
      e = p.second;
      Qit[e] = Q.end();
      std::vector<int> N = igl::circulation(e, true, EMAP, EF, EI);
      std::vector<int> Nd = igl::circulation(e, false, EMAP, EF, EI);
      N.insert(N.begin(), Nd.begin(), Nd.end());

 

      bool collapsed = this->collapse(e, C.row(e), V, F, E, EMAP, EF, EI, Qs, e1, e2, f1, f2);


      if (collapsed) {

          // Erase the two, other collapsed edges
          Q.erase(Qit[e1]);
          Qit[e1] = Q.end();
          Q.erase(Qit[e2]);
          Qit[e2] = Q.end();
          // update local neighbors
          // loop over original face neighbors
          for (auto n : N)
          {
              if (F(n, 0) != IGL_COLLAPSE_EDGE_NULL ||
                  F(n, 1) != IGL_COLLAPSE_EDGE_NULL ||
                  F(n, 2) != IGL_COLLAPSE_EDGE_NULL)
              {
                  for (int v = 0; v < 3; v++)
                  {
                      // get edge id
                      const int ei = EMAP(v * F.rows() + n);
                      // erase old entry
                      Q.erase(Qit[ei]);
                      // compute cost and potential placement
                      double cost;
                      RowVectorXd place;
                      cost_and_placement(Qs, ei, V, F, E, EMAP, EF, EI, cost, place);
                      // Replace in queue
                      Qit[ei] = Q.insert(std::pair<double, int>(cost, ei)).first;
                      C.row(ei) = place;
                  }
              }
          }
      }
      else
      {
          // reinsert with infinite weight (the provided cost function must **not**
          // have given this un-collapsable edge inf cost already)
          p.first = std::numeric_limits<double>::infinity();
          Qit[e] = Q.insert(p).first;
      }
      return collapsed;
  }

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
      int& a_f2)
  {
      // Assign this to 0 rather than, say, -1 so that deleted elements will get
                 // draw as degenerate elements at vertex 0 (which should always exist and
                 // never get collapsed to anything else since it is the smallest index)
      using namespace Eigen;
      using namespace std;
      const int eflip = E(e, 0) > E(e, 1);
      // source and destination
      const int s = eflip ? E(e, 1) : E(e, 0);
      const int d = eflip ? E(e, 0) : E(e, 1);

      if (!igl::edge_collapse_is_valid(e, F, E, EMAP, EF, EI))
      {
          return false;
      }

      // Important to grab neighbors of d before monkeying with edges
      const std::vector<int> nV2Fd = igl::circulation(e, !eflip, EMAP, EF, EI);

      // The following implementation strongly relies on s<d
      assert(s < d && "s should be less than d");
      Qs[s] = Qs[s] + Qs[d];
      // move source and destination to midpoint
      V.row(s) = p;
      V.row(d) = p;

      // Helper function to replace edge and associate information with NULL
      const auto& kill_edge = [&E, &EI, &EF](const int e)
      {
          E(e, 0) = IGL_COLLAPSE_EDGE_NULL;
          E(e, 1) = IGL_COLLAPSE_EDGE_NULL;
          EF(e, 0) = IGL_COLLAPSE_EDGE_NULL;
          EF(e, 1) = IGL_COLLAPSE_EDGE_NULL;
          EI(e, 0) = IGL_COLLAPSE_EDGE_NULL;
          EI(e, 1) = IGL_COLLAPSE_EDGE_NULL;
      };

      // update edge info
      // for each flap
      const int m = F.rows();
      for (int side = 0; side < 2; side++)
      {
          const int f = EF(e, side);
          const int v = EI(e, side);
          const int sign = (eflip == 0 ? 1 : -1) * (1 - 2 * side);
          // next edge emanating from d
          const int e1 = EMAP(f + m * ((v + sign * 1 + 3) % 3));
          // prev edge pointing to s
          const int e2 = EMAP(f + m * ((v + sign * 2 + 3) % 3));
          assert(E(e1, 0) == d || E(e1, 1) == d);
          assert(E(e2, 0) == s || E(e2, 1) == s);
          // face adjacent to f on e1, also incident on d
          const bool flip1 = EF(e1, 1) == f;
          const int f1 = flip1 ? EF(e1, 0) : EF(e1, 1);
          assert(f1 != f);
          assert(F(f1, 0) == d || F(f1, 1) == d || F(f1, 2) == d);
          // across from which vertex of f1 does e1 appear?
          const int v1 = flip1 ? EI(e1, 0) : EI(e1, 1);
          // Kill e1
          kill_edge(e1);
          // Kill f
          F(f, 0) = IGL_COLLAPSE_EDGE_NULL;
          F(f, 1) = IGL_COLLAPSE_EDGE_NULL;
          F(f, 2) = IGL_COLLAPSE_EDGE_NULL;
          // map f1's edge on e1 to e2
          assert(EMAP(f1 + m * v1) == e1);
          EMAP(f1 + m * v1) = e2;
          // side opposite f2, the face adjacent to f on e2, also incident on s
          const int opp2 = (EF(e2, 0) == f ? 0 : 1);
          assert(EF(e2, opp2) == f);
          EF(e2, opp2) = f1;
          EI(e2, opp2) = v1;
          // remap e2 from d to s
          E(e2, 0) = E(e2, 0) == d ? s : E(e2, 0);
          E(e2, 1) = E(e2, 1) == d ? s : E(e2, 1);
          if (side == 0)
          {
              a_e1 = e1;
              a_f1 = f;
          }
          else
          {
              a_e2 = e1;
              a_f2 = f;
          }
      }

      // finally, reindex faces and edges incident on d. Do this last so asserts
      // make sense.
      //
      // Could actually skip first and last, since those are always the two
      // collpased faces.
      for (auto f : nV2Fd)
      {
          for (int v = 0; v < 3; v++)
          {
              if (F(f, v) == d)
              {
                  const int flip1 = (EF(EMAP(f + m * ((v + 1) % 3)), 0) == f) ? 1 : 0;
                  const int flip2 = (EF(EMAP(f + m * ((v + 2) % 3)), 0) == f) ? 0 : 1;
                  assert(
                      E(EMAP(f + m * ((v + 1) % 3)), flip1) == d ||
                      E(EMAP(f + m * ((v + 1) % 3)), flip1) == s);
                  E(EMAP(f + m * ((v + 1) % 3)), flip1) = s;
                  assert(
                      E(EMAP(f + m * ((v + 2) % 3)), flip2) == d ||
                      E(EMAP(f + m * ((v + 2) % 3)), flip2) == s);
                  E(EMAP(f + m * ((v + 2) % 3)), flip2) = s;
                  F(f, v) = s;
                  break;
              }
          }
      }
      // Finally, "remove" this edge and its information
      kill_edge(e);

      return true;
  }
  IGL_INLINE bool Viewer::check_if_infinity(double first) {
      return first == std::numeric_limits<double>::infinity();

  }





  IGL_INLINE bool Viewer::load_mesh_from_file(
      const std::string & mesh_file_name_string)
  {

    // Create new data slot and set to selected
    if(!(data().F.rows() == 0  && data().V.rows() == 0))
    {
      append_mesh();
    }
    data().clear();

    size_t last_dot = mesh_file_name_string.rfind('.');
    if (last_dot == std::string::npos)
    {
      std::cerr<<"Error: No file extension found in "<<
        mesh_file_name_string<<std::endl;
      return false;
    }

    std::string extension = mesh_file_name_string.substr(last_dot+1);

    if (extension == "off" || extension =="OFF")
    {
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;
      if (!igl::readOFF(mesh_file_name_string, V, F))
        return false;
      data().set_mesh(V,F);
    }
    else if (extension == "obj" || extension =="OBJ")
    {
      Eigen::MatrixXd corner_normals;
      Eigen::MatrixXi fNormIndices;

      Eigen::MatrixXd UV_V;
      Eigen::MatrixXi UV_F;
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;

      if (!(
            igl::readOBJ(
              mesh_file_name_string,
              V, UV_V, corner_normals, F, UV_F, fNormIndices)))
      {
        return false;
      }

      data().set_mesh(V,F);
      if (UV_V.rows() > 0)
      {
          data().set_uv(UV_V, UV_F);
      }

    }
    else
    {
      // unrecognized file type
      printf("Error: %s is not a recognized file type.\n",extension.c_str());
      return false;
    }

    data().compute_normals();
    data().uniform_colors(Eigen::Vector3d(51.0/255.0,43.0/255.0,33.3/255.0),
                   Eigen::Vector3d(255.0/255.0,228.0/255.0,58.0/255.0),
                   Eigen::Vector3d(255.0/255.0,235.0/255.0,80.0/255.0));

    // Alec: why?
    if (data().V_uv.rows() == 0)
    {
      data().grid_texture();
    }
    

    //for (unsigned int i = 0; i<plugins.size(); ++i)
    //  if (plugins[i]->post_load())
    //    return true;

    return true;
  }

  IGL_INLINE bool Viewer::save_mesh_to_file(
      const std::string & mesh_file_name_string)
  {
    // first try to load it with a plugin
    //for (unsigned int i = 0; i<plugins.size(); ++i)
    //  if (plugins[i]->save(mesh_file_name_string))
    //    return true;

    size_t last_dot = mesh_file_name_string.rfind('.');
    if (last_dot == std::string::npos)
    {
      // No file type determined
      std::cerr<<"Error: No file extension found in "<<
        mesh_file_name_string<<std::endl;
      return false;
    }
    std::string extension = mesh_file_name_string.substr(last_dot+1);
    if (extension == "off" || extension =="OFF")
    {
      return igl::writeOFF(
        mesh_file_name_string,data().V,data().F);
    }
    else if (extension == "obj" || extension =="OBJ")
    {
      Eigen::MatrixXd corner_normals;
      Eigen::MatrixXi fNormIndices;

      Eigen::MatrixXd UV_V;
      Eigen::MatrixXi UV_F;

      return igl::writeOBJ(mesh_file_name_string,
          data().V,
          data().F,
          corner_normals, fNormIndices, UV_V, UV_F);
    }
    else
    {
      // unrecognized file type
      printf("Error: %s is not a recognized file type.\n",extension.c_str());
      return false;
    }
    return true;
  }
 
  IGL_INLINE bool Viewer::load_scene()
  {
    std::string fname = igl::file_dialog_open();
    if(fname.length() == 0)
      return false;
    return load_scene(fname);
  }

  IGL_INLINE bool Viewer::load_scene(std::string fname)
  {
   // igl::deserialize(core(),"Core",fname.c_str());
    igl::deserialize(data(),"Data",fname.c_str());
    return true;
  }

  IGL_INLINE bool Viewer::save_scene()
  {
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
      return false;
    return save_scene(fname);
  }

  IGL_INLINE bool Viewer::save_scene(std::string fname)
  {
    //igl::serialize(core(),"Core",fname.c_str(),true);
    igl::serialize(data(),"Data",fname.c_str());

    return true;
  }

  IGL_INLINE void Viewer::open_dialog_load_mesh()
  {
    std::string fname = igl::file_dialog_open();

    if (fname.length() == 0)
      return;
    
    this->load_mesh_from_file(fname.c_str());
  }

  IGL_INLINE void Viewer::open_dialog_save_mesh()
  {
    std::string fname = igl::file_dialog_save();

    if(fname.length() == 0)
      return;

    this->save_mesh_to_file(fname.c_str());
  }

  IGL_INLINE ViewerData& Viewer::data(int mesh_id /*= -1*/)
  {
    assert(!data_list.empty() && "data_list should never be empty");
    int index;
    if (mesh_id == -1)
      index = selected_data_index;
    else
      index = mesh_index(mesh_id);

    assert((index >= 0 && index < data_list.size()) &&
      "selected_data_index or mesh_id should be in bounds");
    return data_list[index];
  }

  IGL_INLINE const ViewerData& Viewer::data(int mesh_id /*= -1*/) const
  {
    assert(!data_list.empty() && "data_list should never be empty");
    int index;
    if (mesh_id == -1)
      index = selected_data_index;
    else
      index = mesh_index(mesh_id);

    assert((index >= 0 && index < data_list.size()) &&
      "selected_data_index or mesh_id should be in bounds");
    return data_list[index];
  }

  IGL_INLINE int Viewer::append_mesh(bool visible /*= true*/)
  {
    assert(data_list.size() >= 1);

    data_list.emplace_back();
    selected_data_index = data_list.size()-1;
    data_list.back().id = next_data_id++;
    //if (visible)
    //    for (int i = 0; i < core_list.size(); i++)
    //        data_list.back().set_visible(true, core_list[i].id);
    //else
    //    data_list.back().is_visible = 0;
    return data_list.back().id;
  }

  IGL_INLINE bool Viewer::erase_mesh(const size_t index)
  {
    assert((index >= 0 && index < data_list.size()) && "index should be in bounds");
    assert(data_list.size() >= 1);
    if(data_list.size() == 1)
    {
      // Cannot remove last mesh
      return false;
    }
    data_list[index].meshgl.free();
    data_list.erase(data_list.begin() + index);
    if(selected_data_index >= index && selected_data_index > 0)
    {
      selected_data_index--;
    }

    return true;
  }

  IGL_INLINE size_t Viewer::mesh_index(const int id) const {
    for (size_t i = 0; i < data_list.size(); ++i)
    {
      if (data_list[i].id == id)
        return i;
    }
    return 0;
  }



  // ASSIGNMENT 1 - TASK 4

  IGL_INLINE bool Viewer::load_mesh_from_configuration(const std::string config) {
    
          std::string mesh_path;
          std::fstream infile;
          infile.open(config);
          if (!infile) {
              std::cout << "Can't open file configuration.txt\n";
              return false;
          }
          else {
              while (getline(infile, mesh_path)) {
                  std::cout << "opening " << mesh_path << std::endl;
                  this->load_mesh_from_file(mesh_path);
                  if (simplification_enable) this->init_objs_simpelified();
              }
              infile.close();
              return true;
          }
      
  }

  Eigen::Matrix4d Viewer::CalcParentsTrans(int indx) 
  {
	  Eigen::Matrix4d prevTrans = Eigen::Matrix4d::Identity();

	  for (int i = indx; parents[i] >= 0; i = parents[i])
	  {
		  //std::cout << "parent matrix:\n" << scn->data_list[scn->parents[i]].MakeTrans() << std::endl;
		  prevTrans = data_list[parents[i]].MakeTransd() * prevTrans;
	  }

	  return prevTrans;
  }

} // end namespace
} // end namespace

}
