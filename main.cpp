/*
 * @file main.cpp
 * @brief This code computes the rotation matrix
 * to use to deform vertices 
 *
 * @author Arya Sinha, Shantanu Dixit
 * @date November 10, 2023
 */

#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/columnize.h>
#include <data_holder.h>
#include <iostream>

using namespace Eigen;
using namespace std;
using namespace igl;

void perturb(MatrixXd & V, MatrixXd & U, data_holder & data);

int main(int argc, char *argv[])
{
    
  MatrixXd V, V_tilde;
	MatrixXi F;
  
  string name = argv[1];
  string mesh = "/home/arya/Downloads/project_mid_eval_2020498_2020118/meshes/" + name;
  readOBJ(mesh, V, F);

  double lambda = stod(argv[2]);
  
  V_tilde = V;
  data_holder perturbed_data(V, F, lambda);

  for (int i=0; i<5; i++)
  {
      std::cout<<i<<"\n";
      perturb(V, V_tilde, perturbed_data);
  }

  opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V_tilde, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}

void perturb(MatrixXd & V, MatrixXd & V_tilde, data_holder & data) {

  //obtain rotation matrix from local step
  MatrixXd R(3,V.rows()*3);
  data.local_step(V, V_tilde, R);
    
  //obtain new vertices from global step (using R)
  
}


