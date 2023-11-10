#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>

#include <data_holder.h>
#include <iostream>

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
  //Command to run: ./binary_file [mesh] [lambda] 
  string MESH_DIRECTORY = "/home/iiitd/Documents/arya/libigl-example-project-cmd-version/meshes/";
    
  MatrixXd V, V_tilde;
	MatrixXi F;
  double lambda;
  
	string meshName = argv[1];
  string mesh = MESH_DIRECTORY + meshName;
  igl::readOBJ(mesh, V, F);

  lambda = stod(argv[2]);
  
  V_tilde = V;
  data_holder perturbed_data(V, F, lambda);
  
}
