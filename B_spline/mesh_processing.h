#pragma once
#include<igl/read_triangle_mesh.h>
void read_and_visual_mesh(const std::string &filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

// map the vertices loop to a [0,1]x[0,1] square
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

// this is to remove some faces of a mesh accroding to the vertices coordinates.
// axis select from 0, 1, 2. 
void remove_some_faces(const int axis, const double value, const bool remove_larger,
	const Eigen::MatrixXd& V, const Eigen::MatrixXi F, Eigen::MatrixXi& newF);

// an example for parameterization
void mesh_parameterization(
	const std::string &meshfile, Eigen::MatrixXd& V, Eigen::MatrixXd &param, Eigen::MatrixXi &F);