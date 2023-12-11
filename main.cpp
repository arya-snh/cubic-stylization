/*
 * @file main.cpp
 * @brief This code computes the rotation matrix
 * to use to deform vertices 
 *
 * @author Arya Sinha, Shantanu Dixit
 * @date November 10, 2023
 */
#include <fstream>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/columnize.h>
#include <data_holder.h>
#include <iostream>

using namespace Eigen;
using namespace std;
using namespace igl;

double perturb(MatrixXd & V, MatrixXd & U, data_holder & data);

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
  std::ofstream outFile;
  outFile.open("../objective_values_stopping_criteria.csv", std::ios::out | std::ios::app);

  for (int i=0; i<10; i++)
  {
      std::cout<<i<<"\n";
      double objVal = perturb(V, V_tilde, perturbed_data);
      outFile << i << ", " << objVal << "\n"; 
  }

  outFile.close();

  opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V_tilde, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}

double perturb(MatrixXd & V, MatrixXd & V_tilde, data_holder & data) {

  //obtain rotation matrix from local step
  MatrixXd R(3,V.rows()*3);
  data.local_step(V, V_tilde, R);
    
  //obtain new vertices from global step (using R)
  VectorXd Rcol;
  columnize(R, V.rows(), 2, Rcol); //column wise split
  VectorXd Bcol = data.K * Rcol; //

  // For x, y and z coord
  for(int k=0; k<3; k++)
  {
      VectorXd V_tilde_k, B_k, bc_k;
      B_k = Bcol.block(k*V.rows(), 0, V.rows(), 1);
      bc_k = data.boundary_condition.col(k);
      min_quad_with_fixed_solve(data.solver_data, B_k, bc_k, VectorXd(), V_tilde_k);
      V_tilde.col(k) = V_tilde_k;
  }

  return data.objVal;
  
}


