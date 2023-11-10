#include "data_holder.h"

data_holder::data_holder(Eigen::MatrixXd& V, Eigen::MatrixXi& F, double _lambda) {
    
    mu = 10;
    tao = 2;
    lambda = _lambda;
    maxi = 100;

    boundary_condition.resize(1,3);
    boundary_condition << V.row(F(1,1));

    q.resize(1);
    q << F(1,1);

    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
    barycentric_area = M.diagonal();

    igl::per_vertex_normals(V,F, per_vertex_normals);
    igl::cotmatrix(V,F,cotangent_matrix);

    vector<vector<int>> adj_F, nI;
    igl::vertex_triangle_adjacency(V.rows(),F,adj_F, nI);

    igl::arap_rhs(V,F,V.cols(),igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS, K);

    half_edges.resize(V.rows());
    W.resize(V.rows());
    vector<int> adjFacei;

    int v0, v1, v2; //vertices forming a face
    //For all vertices
    for (int i=0; i<V.rows();  i++)
    {   
        adjFacei = adj_F[i];

        half_edges[i].resize(adjFacei.size()*3, 2);
        W[i].resize(adjFacei.size()*3);
        for (int j=0; j<adjFacei.size(); j++)
        {
            //For all faces adjacent to the vertex
            v0 = F(adjFacei[j],0);
            v1 = F(adjFacei[j],1);
            v2 = F(adjFacei[j],2);

            // compute adjacent half-edge indices of a vertex
            half_edges[i](3*j  ,0) = v0;
            half_edges[i](3*j  ,1) = v1;
            half_edges[i](3*j+1,0) = v1;
            half_edges[i](3*j+1,1) = v2;
            half_edges[i](3*j+2,0) = v2;
            half_edges[i](3*j+2,1) = v0;

            W[i](3*j  ) = cotangent_matrix.coeff(v0,v1);
            W[i](3*j+1) = cotangent_matrix.coeff(v1,v2);
            W[i](3*j+2) = cotangent_matrix.coeff(v2,v0);
        }
    }

    igl::min_quad_with_fixed_precompute(cotangent_matrix, q, SparseMatrix<double>(), false, solver_data);

    z_a.resize(3, V.rows()); z_a.setRandom();
    u_a.resize(3, V.rows()); u_a.setRandom();
    rho_a.resize(V.rows()); rho_a.setConstant(1e-4);

}

void data_holder::local_step(const Eigen::MatrixXd & V, Eigen::MatrixXd & U, Eigen::MatrixXd & RAll)
{
    igl::parallel_for(V.rows(), [this, &RAll, &U, &V](const int i)
        {
            VectorXd z = z_a.col(i);
            VectorXd u = u_a.col(i);
            VectorXd n = per_vertex_normals.row(i).transpose();
            double rho = rho_a(i);
            Matrix3d R;
            
            MatrixXi hE = half_edges[i];
            MatrixXd dU(3, hE.rows()); 
            
            MatrixXd U_hE0, U_hE1;
            igl::slice(U,hE.col(0),1,U_hE0);
            igl::slice(U,hE.col(1),1,U_hE1);
            dU = (U_hE1 - U_hE0).transpose();
            
            MatrixXd dV(3, hE.rows()); 
            
            MatrixXd V_hE0, V_hE1;
            igl::slice(V,hE.col(0),1,V_hE0);
            igl::slice(V,hE.col(1),1,V_hE1);
            dV = (V_hE1 - V_hE0).transpose();
            

            Matrix3d M_ = dV * W[i].asDiagonal() * dU.transpose();

            // Scaled Dual ADMM
            for (int k=0; k<maxi; k++)
            {
                // update R
                Matrix3d M = M_ + (rho * n * (z-u).transpose());
                
                JacobiSVD<Matrix3d> svd;
                svd.compute(M, Eigen::ComputeFullU | Eigen::ComputeFullV );
                Matrix3d u__ = svd.matrixU();
                Matrix3d v__ = svd.matrixV();
                R = v__ * u__.transpose();
                if (R.determinant() < 0)
                {
                    u__.col(2) = -u__.col(2);
                    R = v__ * u__.transpose();
                }

                // update z
                VectorXd z_prev = z;

                Eigen::VectorXd x = R*n+u;
                double kk = lambda*barycentric_area(i)/rho;
                VectorXd tmp1 = x.array() - kk;
                VectorXd tmp2 = -x.array() - kk;

                z = tmp1.array().max(0.0) - tmp2.array().max(0.0);

                // update u
                u.noalias() += R*n - z;

                // compute residual
                double r_norm = (z - R*n).norm();
                double s_norm = (-rho * (z - z_prev)).norm();
                
                // from boyd
                // update rho & u
                if (r_norm > mu * s_norm)
                {
                    rho = tao * rho;
                    u = u / tao;
                }
                else if (s_norm > mu * r_norm)
                {
                    rho = rho / tao;
                    u = u * tao;
                }

            }
            z_a.col(i) = z;
            u_a.col(i) = u;
            rho_a(i) = rho;
            RAll.block(0,3*i,3,3) = R; 
        }   
    ,1000);
}