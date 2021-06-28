#pragma once
#include<igl/read_triangle_mesh.h>
#include<Types.hpp>
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

void mesh_parameter_perturbation(const Eigen::MatrixXd &para_in, const Eigen::MatrixXi &F,
	Eigen::MatrixXd &para_out, int itr);
// given a 2d point set V, detect its borders loop
void find_border_loop(const Eigen::MatrixXd& V, Eigen::VectorXi& loop);
void constrained_delaunay_triangulation(
	const Eigen::MatrixXd& V, const Eigen::VectorXi& Edge_ids,
	Eigen::MatrixXi& F);
void smooth_mesh_vertices(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
	const double h, const std::vector<int> &fixed,
	const int itrs, Eigen::MatrixXd&Vout);
// given parameters and the connectivity, get the U and V parameters, and a map showing the positions of the points
// in U and V
void generate_UV_grid(const Eigen::MatrixXd& param, 
	std::vector<double>& U, std::vector<double>&V, Eigen::MatrixXi& map);

// statement sees function find_one_ring_for_param()
Vector3d linear_interpolation(const int Fid, const double u, const double v,
	const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& param,
	const Eigen::MatrixXi& F);

//  remesh the grid for the smoothing.
void remeshing_based_on_map_grid(const Eigen::MatrixXd& param, const Eigen::MatrixXd& vertices,
	const Eigen::MatrixXi& F,
	const std::vector<double>& U, const std::vector<double>&V, const Eigen::MatrixXi& map,
	Eigen::MatrixXd& paramout, Eigen::MatrixXd& ver_out, Eigen::MatrixXi& Fout);

// the returned value is sqrt distance. but sqrtD is un-squared
double point_mesh_distance(const Eigen::MatrixXd ver, const Eigen::MatrixXd& Vmesh,
	const Eigen::MatrixXi& Fmesh, Eigen::VectorXd &sqrD);