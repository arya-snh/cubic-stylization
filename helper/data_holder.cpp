#include "data_holder.h"

data_holder::data_holder(Eigen::MatrixXd& V, Eigen::MatrixXi& F, double _lambda) {
    
    mu = 10;
    tao = 2;
    lambda = _lambda;
    maxi = 100;

    boundary_condition.resize(1,3);
    boundary_condition << V.row(F(0,0));

    q.resize(1);
    q << F(0,0);

    igl::per_vertex_normals(V,F, per_vertex_normals);
    igl::cotmatrix(V,F,cotangent_matrix);

    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
    barycentric_area = M.diagonal();

    vector<vector<int>> adj_F, nI;
    igl::vertex_triangle_adjacency(V.rows(),F,adj_F, nI);

    igl::arap_rhs(V,F,V.cols(),igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS, K);

    hEList.resize(V.rows());
    W.resize(V.rows());
    dVList.resize(V.rows());
    vector<int> adjFacei;

    int v0, v1, v2; //vertices forming a face
    //For all vertices
    for (int i=0; i<V.rows();  i++)
    {   
        adjFacei = adj_F[i];

        hEList[i].resize(adjFacei.size()*3, 2);
        W[i].resize(adjFacei.size()*3);
        for (int j=0; j<adjFacei.size(); j++)
        {
            //For all faces adjacent to the vertex
            v0 = F(adjFacei[j],0);
            v1 = F(adjFacei[j],1);
            v2 = F(adjFacei[j],2);

            // compute adjacent half-edge indices of a vertex
            hEList[i](3*j  ,0) = v0;
            hEList[i](3*j  ,1) = v1;
            hEList[i](3*j+1,0) = v1;
            hEList[i](3*j+1,1) = v2;
            hEList[i](3*j+2,0) = v2;
            hEList[i](3*j+2,1) = v0;

            W[i](3*j  ) = cotangent_matrix.coeff(v0,v1);
            W[i](3*j+1) = cotangent_matrix.coeff(v1,v2);
            W[i](3*j+2) = cotangent_matrix.coeff(v2,v0);
        }

        dVList[i].resize(3, adjFacei.size()*3);
        MatrixXd V_hE0, V_hE1;
        igl::slice(V,hEList[i].col(0),1,V_hE0);
        igl::slice(V,hEList[i].col(1),1,V_hE1);
        dVList[i] = (V_hE1 - V_hE0).transpose();
    }

    igl::min_quad_with_fixed_precompute(cotangent_matrix, q, SparseMatrix<double>(), false, solver_data);

    z_a.resize(3, V.rows()); z_a.setRandom();
    u_a.resize(3, V.rows()); u_a.setRandom();
    rho_a.resize(V.rows()); rho_a.setConstant(1e-4);

}

void data_holder::local_step(const Eigen::MatrixXd & V, Eigen::MatrixXd & U, Eigen::MatrixXd & RAll)
{
    igl::parallel_for(V.rows(), [this, &RAll, &U](const int ii)
        {
            VectorXd z = z_a.col(ii);
            VectorXd u = u_a.col(ii);
            VectorXd n = per_vertex_normals.row(ii).transpose();
            double rho = rho_a(ii);
            Matrix3d R;

            // get energy parameters
            // Note: dVn = [dV n], dUn = [dU z-u]
            MatrixXi hE = hEList[ii];
            MatrixXd dU(3,hE.rows()); 
            {
                MatrixXd U_hE0, U_hE1;
                igl::slice(U,hE.col(0),1,U_hE0);
                igl::slice(U,hE.col(1),1,U_hE1);
                dU = (U_hE1 - U_hE0).transpose();
            }

            MatrixXd dV = dVList[ii];
            VectorXd WVec = W[ii];
            Matrix3d Spre = dV * WVec.asDiagonal() * dU.transpose();

            // Scaled Dual ADMM
            for (int k=0; k<maxi; k++)
            {
                // R step
                Matrix3d S = Spre + (rho * n * (z-u).transpose());
                
                JacobiSVD<Matrix3d> svd;
                svd.compute(S, Eigen::ComputeFullU | Eigen::ComputeFullV );
                Matrix3d SU = svd.matrixU();
                Matrix3d SV = svd.matrixV();
                R = SV * SU.transpose();
                if (R.determinant() < 0)
                {
                    SU.col(2) = -SU.col(2);
                    R = SV * SU.transpose();
                }

                assert(R.determinant() > 0);

                // z step
                VectorXd zOld = z;

                Eigen::VectorXd x = R*n+u;
                double kk = lambda* barycentric_area(ii)/rho;
                VectorXd tmp1 = x.array() - kk;
                VectorXd posMax = tmp1.array().max(0.0);

                VectorXd tmp2 = -x.array() - kk;
                VectorXd negMax = tmp2.array().max(0.0);

                z = posMax - negMax;

                // u step
                u.noalias() += R*n - z;

                // compute residual
                double r_norm = (z - R*n).norm();
                double s_norm = (-rho * (z - zOld)).norm();
                
                // rho step
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
            z_a.col(ii) = z;
            u_a.col(ii) = u;
            rho_a(ii) = rho;
            RAll.block(0,3*ii,3,3) = R; 
        }   
    ,1000);
}