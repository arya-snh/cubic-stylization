#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/columnize.h>
#include <data_holder.h>
#include <iostream>

using namespace Eigen;
using namespace std;

void perturb(MatrixXd & V, MatrixXd & U, data_holder & data);

int main(int argc, char *argv[])
{
  //Command to run: ./binary_file [mesh] [lambda] 
  string MESH_DIRECTORY = "/home/iiitd/Documents/cubic-stylization/meshes/";
    
  MatrixXd V, V_tilde;
	MatrixXi F;
  double lambda;
  
	string meshName = argv[1];
  string mesh = MESH_DIRECTORY + meshName;
  igl::readOBJ(mesh, V, F);

  lambda = stod(argv[2]);
  
  V_tilde = V;
  data_holder perturbed_data(V, F, lambda);

  for (int i=0; i<5; i++)
  {
      std::cout<<i<<"\n";
      perturb(V, V_tilde, perturbed_data);
  }

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V_tilde, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}

void perturb(MatrixXd & V, MatrixXd & V_tilde, data_holder & data) {

  //obtain rotation matrix from local step
  MatrixXd R(3,V.rows()*3);
  data.local_step(V, V_tilde, R);
    
  //obtain new vertices from global step (using R)
  VectorXd Rcol;
  igl::columnize(R, V.rows(), 2, Rcol);
  VectorXd Bcol = data.K * Rcol;
  for(int k=0; k<V.cols(); k++)
  {
      VectorXd V_tilde_c, Bc, bcc;
      Bc = Bcol.block(k*V.rows(), 0, V.rows(), 1);
      bcc = data.boundary_condition.col(k);
      min_quad_with_fixed_solve(data.solver_data, Bc, bcc, VectorXd(), V_tilde_c);
      V_tilde.col(k) = V_tilde_c;
  }
    
}


