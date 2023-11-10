#ifndef CUBE_STYLE_DATA_H
#define CUBE_STYLE_DATA_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <limits>
#include <igl/min_quad_with_fixed.h>

class data_holder
{
    public:

	double lambda = 0.0;


    //Boyd penalty update parameters
	double mu = 10;
	double tao = 2;

    //maximum iterations for ADMM
	double maxi = 100;

	std::vector<Eigen::MatrixXi> hEList;
	std::vector<Eigen::MatrixXd> dVList;
	std::vector<Eigen::VectorXd> WVecList;

	Eigen::SparseMatrix<double> K, L;
	Eigen::MatrixXd N, VA, zAll, uAll;
	Eigen::VectorXd rhoAll;

	Eigen::MatrixXd bc; //boundary condition
	Eigen::VectorXi b;  //fixed vertex

	igl::min_quad_with_fixed_data<double> solver_data;
};

#endif