#pragma once
#include<igl/read_triangle_mesh.h>
#include<sparse_interp/Types.hpp>
namespace SIBSplines{
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


// given a 2d point set V, detect its borders loop
void find_border_loop(const Eigen::MatrixXd& V, Eigen::VectorXi& loop);
void constrained_delaunay_triangulation(
	const Eigen::MatrixXd& V, const Eigen::VectorXi& Edge_ids,
	Eigen::MatrixXi& F);

// given parameters and the connectivity, get the U and V parameters, and a map showing the positions of the points
// in U and V
void generate_UV_grid(const Eigen::MatrixXd& param, 
	std::vector<double>& U, std::vector<double>&V, Eigen::MatrixXi& map);


// the returned value is sqrt distance. but sqrtD is un-squared
double point_mesh_distance(const Eigen::MatrixXd ver, const Eigen::MatrixXd& Vmesh,
	const Eigen::MatrixXi& Fmesh, Eigen::VectorXd &sqrD);

}