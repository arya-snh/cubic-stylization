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