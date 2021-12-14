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
      
    //  load_mesh_from_configuration(config);

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

  /* 
  in order to perform a certian contraction we will need to calculate the cost of a contraction 
  to define this cost we attemp to characterize the error at each vertex in the mesh . 

  to do this : we associate a symmatric 4x4 matrix named Q with each vertex .

  and we define the error at a vertix v to be 

  v = = [vx vy vz 1]T  to be the quadratic form tirangle(v) = vTQv. 

  meaning that the error is equal to vTQv. 

  let symbol and name a given constraction (V1,v2) as v~ 

  * now we mush derive a new matrix that we call Q~ which approximates the error at v~ 
  * 
  *  the rule is that Q~ = Q1 + Q2 - WHERE Q1 is the associted matrix for v1 and Q2 is the associated matrix of v2 
  * 
  * 
  * in order to do the constraction we must find a position for v~  - the simple scheme would be to either choose 
  * 1. v1 
  * 2. v2 
  * 3. v1+v2 \2  
  * 
  * 
  * now it depends on which of these have the lowest value at triangelV 
  * 
  * 
  *  meaning which of the following has the lowest value 
  * 
  * 1. v1T X  Q1 X v1
  * 
  * 2. v2T X  Q2 X v2 
  * 
  * 3. ((v1+v2)/2)T  x Q~ x (v1+v2)  
  * 
  * 
  * now it is mentioned how to calculate the v~ 
  * 
  * 
  * 
  * 
  * 
  * 


  
  
  
  
  
  
  */

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
          // this is the v~ 
          a = Qi * Eigen::Matrix<double, 4, 1>(0, 0, 0, 1);
      }
      else {
          // we need to choose the midpoint 
          a = (V.row(E.row(e)[0]) + V.row(E.row(e)[1])) / 2;
          a(3) = 1;
      }


      cost = a.transpose() * (Q1 + Q2) * a;
      // new point s
      p = Eigen::RowVector3d(a(0), a(1), a(2));
  }

  //Assingment 1 task 6
  IGL_INLINE void Viewer::init_objs_simpelified() {
     
      data().Q = new std::set<std::pair<double, int> >(); // priority Q contains the cost of edges <cost , number of edge > 
      data().EMAP = new Eigen::VectorXi();               // connects faces to edges 

      /* 
      E - matrix of edges
      |1 , 3 |
      |2 , 5 |
      |3, 1  |
      
       source , destenation 
      
      
      */

      data().E = new Eigen::MatrixXi();                   // this is the Edges -> edge is represented by <index of source vertex, index of destinations vertex> - matrix []
      
      data().EI = new Eigen::MatrixXi();                  //connects edge to vertex index in triangle (0 1 2) 


      data().EF = new Eigen::MatrixXi();                    //connects edges to faces this is a matrix that have [index of first face , index of second face ]
    
      
      data().C = new Eigen::MatrixXd();                 //position of the new vertex after collapsing the corresponding edge


      data().V1 = new Eigen::MatrixXd(data().V);            // the new back up for v 

      data().F1 = new Eigen::MatrixXi(data().F);            // the new backup for F 

      data().Qs = std::vector<Eigen::Matrix<double, 4, 4>>();
      data().Qit = new std::vector<std::set<std::pair<double, int> >::iterator >(); // saved iterator for every dege so that 
   

      //this gets the EMAP EF E EI ready 
      edge_flaps(*data().F1, *data().E, *data().EMAP, *data().EF, *data().EI);

      data().Qit->resize(data().E->rows());
     
      data().C->resize(data().E->rows(), data().V.cols());
   
      Eigen::VectorXd costs(data().E->rows());
  
      data().Q->clear();
     


      auto VF = std::vector<std::vector<int> >(); 
      auto VFi = std::vector<std::vector<int> >();

      vertex_triangle_adjacency(data().V, data().F, VF, VFi);
   
      // this is for getting the Q ready 
      for (int v = 0; v < data().V.rows(); v++) {

       
          std::vector<int> faces = VF[v];
          Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
          
          for (int face : faces) {
              auto n = data().F_normals.row(face).normalized();
              float d = 0;
              for (int j = 0; j < 3; j++) 
              d += (-n[j] * data().V.row(v)[j]);
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
          double cost_of_Edge = p.first;

          printf("we collapsed the following edges : e -> %d , cost : %f     e1 -> %d  , cost: %f   e2 ->  , cost:    ", e, cost_of_Edge, e1, Qit[e1]->first);
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
                   
                      const int ei = EMAP(v * F.rows() + n);
                  
                      Q.erase(Qit[ei]);
                   
                      double cost;
                      RowVectorXd place;
                      cost_and_placement(Qs, ei, V, F, E, EMAP, EF, EI, cost, place);
                  
                      Qit[ei] = Q.insert(std::pair<double, int>(cost, ei)).first;
                      C.row(ei) = place;
                  }
              }
          }
      }
      else
      {
    
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

  IGL_INLINE bool Viewer::load_mesh_from_configuration(const std::string config , bool assignment2) {
    
      if (assignment2)
          return true;
          //init_obj_collision();
     
     
      else {
          std::string mesh_path;
          std::fstream infile;
          infile.open("configuration.txt");
          if (!infile) {
              std::cout << "Can't open file configuration.txt\n";
              return false;
          }
          else {
              while (getline(infile, mesh_path)) {
                  std::cout << "opening " << mesh_path << std::endl;
                  this->load_mesh_from_file(mesh_path);
                 // if (enable_simplefication) this->init_simplefication_objs();
              }
              infile.close();
              return true;
          }
      }
  }
    


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// ASSIGNMENT 2 //////////////////////////////////////////////////////////////////////




  IGL_INLINE void Viewer::start_collision_render() {
      this->data_list[1].MyTranslate(Eigen::Vector3d(-0.005, 0, 0), true);

      Eigen::Matrix4f t1 = MakeTransScale() * data_list[0].MakeTransScale();
      auto a1 = t1.row(0);
      auto a2 = t1.row(1);
      auto a3 = t1.row(2);
      Eigen::Vector3d A[3] =
      { Eigen::Vector3d(a1(0),a1(1),a1(2)),
          Eigen::Vector3d(a2(0),a2(1),a2(2)),
          Eigen::Vector3d(a3(0),a3(1),a3(2)) };

       t1 = MakeTransScale() * data_list[1].MakeTransScale();
       a1 = t1.row(0);
       a2 = t1.row(1);
       a3 = t1.row(2);
      Eigen::Vector3d B[3] =
      { Eigen::Vector3d(a1(0),a1(1),a1(2)),
          Eigen::Vector3d(a2(0),a2(1),a2(2)),
          Eigen::Vector3d(a3(0),a3(1),a3(2)) };


      if (check_collision_occurence(this->trees.at(0), this->trees.at(1), A, B)) {
          this->stop_collision();
          this->trees.at(0) = this->tree_roots.at(0);
          this->trees.at(1) = this->tree_roots.at(1);
      }
     
  }

  IGL_INLINE bool Viewer::check_collision_occurence(igl::AABB<Eigen::MatrixXd, 3>* first, igl::AABB<Eigen::MatrixXd, 3>* second, Eigen::Vector3d A[], Eigen::Vector3d B[]) {



      bool collision = false;
      igl::AABB<Eigen::MatrixXd, 3>* first_left_child = first->m_left;
      igl::AABB<Eigen::MatrixXd, 3>* first_right_child = first->m_right;
      igl::AABB<Eigen::MatrixXd, 3>* second_left_child = second->m_left;
      igl::AABB<Eigen::MatrixXd, 3>* second_right_child = second->m_right;



      if (first == nullptr || second == nullptr) return collision;

  // if both parents are not leafs- may be at root level and may be deep inside the tree
      if (!(first->is_leaf() || second->is_leaf())) {



                    // if can seprate = true - this means that there is no collision 
          if (!check_possible_seperate(first->m_box, second->m_box, A, B)) {
               // getting here means that the roots / parents couldnt be seperated - there might be a collision we need to check the children 
              // 

              if (!collision && !check_possible_seperate(first_right_child->m_box, second_right_child->m_box, A, B)) {
                  this->last_box.at(0) = &first_right_child->m_box;
                  this->last_box.at(1) = &second_right_child->m_box;

                  collision = check_collision_occurence(first_right_child, second_right_child, A, B);
              }

              if (!collision && !check_possible_seperate(first_right_child->m_box, second_left_child->m_box, A, B)) {
                  this->last_box.at(0) = &first_right_child->m_box;
                  this->last_box.at(1) = &second_left_child->m_box;
                  collision = check_collision_occurence(first_right_child, second_left_child, A, B);
              }
              if (!collision && !check_possible_seperate(first_left_child->m_box, second_right_child->m_box, A, B)) {
                  this->last_box.at(0) = &first_left_child->m_box;
                  this->last_box.at(1) = &second_right_child->m_box;

                  collision = check_collision_occurence(first_left_child, second_right_child, A, B);
              }
            
              if (!collision && !check_possible_seperate(first_left_child->m_box, second_left_child->m_box, A, B)) {
                  this->last_box.at(0) = &first_left_child->m_box;
                  this->last_box.at(1) = &second_left_child->m_box;
                
                  collision = check_collision_occurence(first_left_child, second_left_child, A, B);
              }
            
          }
      }
      else if (!first->is_leaf() && second->is_leaf()) {
          if (!collision && !check_possible_seperate(first_left_child->m_box, second->m_box, A, B)) {
              this->last_box.at(0) = &first_left_child->m_box;

              collision = check_collision_occurence(first_left_child, second, A, B);
          }
          if (!collision && !check_possible_seperate(first_right_child->m_box, second->m_box, A, B)) {
              this->last_box.at(0) = &first_right_child->m_box;

              collision = check_collision_occurence(first_right_child, second, A, B);
          }
      }
      
      else if (first->is_leaf() && !second->is_leaf()) {
          if (!collision && !check_possible_seperate(first->m_box, first_left_child->m_box, A, B)) {
              this->last_box.at(1) = &first_left_child->m_box;
       
              collision = check_collision_occurence(first, first_left_child, A, B);
          }
          if (!collision && !check_possible_seperate(first->m_box, second_right_child->m_box, A, B)) {
              this->last_box.at(1) = &second_right_child->m_box;
           
              collision = check_collision_occurence(first, second_right_child, A, B);
          }
      }
     
      else if (first->is_leaf() && second->is_leaf()) {
          this->last_box.at(0) = &first->m_box;
          this->last_box.at(1) = &second->m_box;
          collision = !check_possible_seperate(first->m_box, second->m_box, A, B);
      }
      if (!collision) {
          this->trees.at(0) = first;
          this->trees.at(1) = second;
      }
      return collision;
  }

  IGL_INLINE bool Viewer::check_possible_seperate(Eigen::AlignedBox<double, 3>& first, Eigen::AlignedBox<double, 3>& second, Eigen::Vector3d A[], Eigen::Vector3d B[]) {

      Eigen::Vector3d C1 = ((this->data_list[0].MakeTransScale()).cast<double>() * Eigen::Vector4d(first.center()(0), first.center()(1), first.center()(2), 1)).block<3, 1>(0, 0);
      Eigen::Vector3d C2 = ((this->data_list[1].MakeTransScale()).cast<double>() * Eigen::Vector4d(second.center()(0), second.center()(1), second.center()(2), 1)).block<3, 1>(0, 0);
      Eigen::Vector3d T = (C2 - C1);

   
      Eigen::Vector3d dims[2] = { first.sizes().cast<double>() / 2  , second.sizes().cast<double>() / 2 };
      for (int dim = 0; dim < 3; dim++) {
          double left_side = (dims[0](dim) + std::abs(dims[1](0) * (A[dim].dot(B[0]))) + std::abs(dims[1](1) * (A[dim].dot(B[1]))) + std::abs(dims[1](2) * (A[dim].dot(B[2]))));
          double right_side = std::abs(T.dot(A[dim]));
          if (right_side > left_side) return true;
      }
      for (int dim = 0; dim < 3; dim++) {
          double left_side = (dims[1](dim) + std::abs(dims[0](0) * (A[0].dot(B[dim]))) + std::abs(dims[0](1) * (A[1].dot(B[dim]))) + std::abs(dims[0](2) * (A[2].dot(B[dim]))));
          double right_side = std::abs(T.dot(B[dim]));
          if (right_side > left_side)return true;
      }

      double left_side = std::abs(dims[0](1) * A[2].dot(B[0])) + std::abs(dims[0](2) * A[1].dot(B[0])) + std::abs(dims[1](1) * A[0].dot(B[2])) + std::abs(dims[1](2) * A[0].dot(B[1]));
      double right_side = std::abs((T.dot(A[2])) * (A[1].dot(B[0])) - (T.dot(A[1])) * (A[2].dot(B[0])));
      if (right_side > left_side) return true;

      left_side = std::abs(dims[0](1) * A[2].dot(B[1])) + std::abs(dims[0](2) * A[1].dot(B[1])) + std::abs(dims[1](0) * A[0].dot(B[2])) + std::abs(dims[1](2) * A[0].dot(B[0]));
      right_side = std::abs((T.dot(A[2])) * (A[1].dot(B[1])) - (T.dot(A[1])) * (A[2].dot(B[1])));
      if (right_side > left_side)return true;

      left_side = std::abs(dims[0](1) * A[2].dot(B[2])) + std::abs(dims[0](2) * A[1].dot(B[2])) + std::abs(dims[1](0) * A[0].dot(B[1])) + std::abs(dims[1](1) * A[0].dot(B[0]));
      right_side = std::abs((T.dot(A[2])) * (A[1].dot(B[2])) - (T.dot(A[1])) * (A[2].dot(B[2])));
      if (right_side > left_side) return true;

      left_side = std::abs(dims[0](0) * A[2].dot(B[0])) + std::abs(dims[0](2) * A[0].dot(B[0])) + std::abs(dims[1](1) * A[1].dot(B[2])) + std::abs(dims[1](2) * A[1].dot(B[1]));
      right_side = std::abs((T.dot(A[0])) * (A[2].dot(B[0])) - (T.dot(A[2])) * (A[0].dot(B[0])));
      if (right_side > left_side) return true;
   
      left_side = std::abs(dims[0](0) * A[2].dot(B[1])) + std::abs(dims[0](2) * A[0].dot(B[1])) + std::abs(dims[1](0) * A[1].dot(B[2])) + std::abs(dims[1](2) * A[1].dot(B[0]));
      right_side = std::abs((T.dot(A[0])) * (A[2].dot(B[1])) - (T.dot(A[2])) * (A[0].dot(B[1])));
      if (right_side > left_side) return true;
   
      left_side = std::abs(dims[0](0) * A[2].dot(B[2])) + std::abs(dims[0](2) * A[0].dot(B[2])) + std::abs(dims[1](0) * A[1].dot(B[1])) + std::abs(dims[1](1) * A[1].dot(B[0]));
      right_side = std::abs((T.dot(A[0])) * (A[2].dot(B[2])) - (T.dot(A[2])) * (A[0].dot(B[2])));
      if (right_side > left_side)return true;

      left_side = std::abs(dims[0](0) * A[1].dot(B[0])) + std::abs(dims[0](1) * A[0].dot(B[0])) + std::abs(dims[1](1) * A[2].dot(B[2])) + std::abs(dims[1](2) * A[2].dot(B[1]));
      right_side = std::abs((T.dot(A[1])) * (A[0].dot(B[0])) - (T.dot(A[0])) * (A[1].dot(B[0])));
      if (right_side > left_side) return true;

      left_side = std::abs(dims[0](0) * A[1].dot(B[1])) + std::abs(dims[0](1) * A[0].dot(B[1])) + std::abs(dims[1](0) * A[2].dot(B[2])) + std::abs(dims[1](2) * A[2].dot(B[0]));
      right_side = std::abs((T.dot(A[1])) * (A[0].dot(B[1])) - (T.dot(A[0])) * (A[1].dot(B[1])));
      if (right_side > left_side)return true;
    
      left_side = std::abs(dims[0](0) * A[1].dot(B[2])) + std::abs(dims[0](1) * A[0].dot(B[2])) + std::abs(dims[1](0) * A[2].dot(B[1])) + std::abs(dims[1](1) * A[2].dot(B[0]));
      right_side = std::abs((T.dot(A[1])) * (A[0].dot(B[2])) - (T.dot(A[0])) * (A[1].dot(B[2])));
      if (right_side > left_side)return true;

      return false;
  }

 

  IGL_INLINE void Viewer::stop_collision(void) {
      Eigen::RowVector3d green_color(1, 1, 0);
      
      this->draw_m_box(0, *this->last_box.at(0), green_color);
      this->draw_m_box(1, *this->last_box.at(1), green_color);
      this->trees.at(0) = this->tree_roots.at(0);
      this->trees.at(1) = this->tree_roots.at(1);

      this->render_collison = false;
    
  }



  IGL_INLINE void Viewer::draw_m_box(int index, Eigen::AlignedBox<double, 3>& m_box, Eigen::RowVector3d color) {
      // Corners of the bounding 
      Eigen::MatrixXd V_box(8, 3);
      V_box << m_box.corner(m_box.BottomLeftCeil).transpose(),
          m_box.corner(m_box.BottomLeftFloor).transpose(),
          m_box.corner(m_box.BottomRightCeil).transpose(),
          m_box.corner(m_box.BottomRightFloor).transpose(),
          m_box.corner(m_box.TopLeftCeil).transpose(),
          m_box.corner(m_box.TopLeftFloor).transpose(),
          m_box.corner(m_box.TopRightCeil).transpose(),
          m_box.corner(m_box.TopRightFloor).transpose();
      // Edges of the bounding box
      Eigen::MatrixXi E_box(12, 2);
      E_box <<
          0, 1,
          1, 3,
          2, 3,
          2, 0,
          4, 5,
          5, 7,
          6, 7,
          6, 4,
          0, 4,
          1, 5,
          2, 6,
          7, 3;
      // Plot the corners of the bounding box as points
      this->data_list[index].add_points(V_box, color);
      // Plot the edges of the bounding box
      for (unsigned i = 0; i < E_box.rows(); ++i)
          this->data_list[index].add_edges(V_box.row(E_box(i, 0)), V_box.row(E_box(i, 1)), color);
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
