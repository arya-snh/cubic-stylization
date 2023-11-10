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
  data_holder perturb;
	string meshName = argv[1];
  string mesh = MESH_DIRECTORY + meshName;
  perturb.lambda = stod(argv[2]);

  igl::readOBJ(mesh, V, F);
  V_tilde = V;

  {
      perturb.bc.resize(1,3);
      perturb.bc << V.row(F(0,0));

      perturb.b.resize(1);
      perturb.b << F(0,0);
  }

}
