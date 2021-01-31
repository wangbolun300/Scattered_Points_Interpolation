#pragma once
#include<igl/read_triangle_mesh.h>
void read_and_visual_mesh(const std::string &filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

void map_vertices_to_square(const Eigen::MatrixXd& V, const Eigen::VectorXi& bnd, Eigen::MatrixXd &bnd_uv);
// remove the vertices not used in contructing the surface
// this function does not consider if there are duplicated vertices (two vertices in 
// exact the same position)

void remove_redundent_mesh_vertices(const Eigen::MatrixXd V, const Eigen::MatrixXi F,
	Eigen::MatrixXd& Vout, Eigen::MatrixXi& Fout);

// parameterization requires the mesh have no unrelated vertices (may cause error) 
// and duplicated vertices(may cause open boundary in unwanted area). this function clean the problem
void generate_clean_mesh_data_for_parametrization(const Eigen::MatrixXd V, const Eigen::MatrixXi F,
	Eigen::MatrixXd& Vout, Eigen::MatrixXi& Fout);