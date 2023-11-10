#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/columnize.h>
#include <data_holder.h>
#include <iostream>

using namespace Eigen;
using namespace std;

void single_iteration(const Eigen::MatrixXd & V, Eigen::MatrixXd & U, data_holder & data);

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
      single_iteration(V, V_tilde, perturbed_data);
  }

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V_tilde, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}

void single_iteration(const Eigen::MatrixXd & V, Eigen::MatrixXd & U, data_holder & data) {

  MatrixXd RAll(3,V.rows()*3);
    {
        data.local_step(V, U, RAll, data);
    }

  MatrixXd Upre = U;
    {
        VectorXd Rcol;
        igl::columnize(RAll, V.rows(), 2, Rcol);
        VectorXd Bcol = data.K * Rcol;
        for(int dim=0; dim<V.cols(); dim++)
        {
            VectorXd Uc,Bc,bcc;
            Bc = Bcol.block(dim*V.rows(),0,V.rows(),1);
            bcc = data.boundary_condition.col(dim);
            min_quad_with_fixed_solve(
                data.solver_data,Bc,bcc,VectorXd(),Uc);
            U.col(dim) = Uc;
        }
    }
}


