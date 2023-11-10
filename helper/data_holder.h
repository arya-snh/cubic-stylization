#ifndef DATA_HOLDER_H
#define DATA_HOLDER_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <limits>
#include <igl/min_quad_with_fixed.h>
#include <vector>

// needed libigl functions
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/arap_rhs.h>
#include <igl/slice.h>
#include <igl/min_quad_with_fixed.h>

using namespace std;
using namespace Eigen;
using namespace igl;

class data_holder
{
    private:
    //degree of cubeness
    double lambda;

    //Boyd penalty update parameters
	double mu, tao;

    //maximum iterations for ADMM
	double maxi;

    std::vector<Eigen::MatrixXi> hEList;
	std::vector<Eigen::MatrixXd> dVList;
	std::vector<Eigen::VectorXd> W;

	Eigen::SparseMatrix<double> K, cotangent_matrix;
	Eigen::MatrixXd per_vertex_normals, barycentric_area, z_a, u_a;
	Eigen::VectorXd rho_a;

	Eigen::MatrixXd boundary_condition;
	Eigen::VectorXi q;  //fixed vertex

	igl::min_quad_with_fixed_data<double> solver_data;

    public:
    data_holder(Eigen::MatrixXd & V, Eigen::MatrixXi & F, double _lambda);
};

#endif