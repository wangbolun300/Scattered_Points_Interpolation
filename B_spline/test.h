#pragma once
#include"curve.h"

void test_curve_knot_fixing();
void visual_curve_fitting(Eigen::MatrixXd& control_pts, Eigen::MatrixXd& control_pts_color,
	Eigen::MatrixXd& curve_pts, Eigen::MatrixXd& curve_pts_color,
	Eigen::MatrixXd& target_pts, Eigen::MatrixXd& target_pts_color);
void test_surface_visual(Eigen::MatrixXd &ver, Eigen::MatrixXi& faces);
void parameter_grid_to_mesh(const Eigen::MatrixXd& uv, Eigen::MatrixXd &ver, Eigen::MatrixXi& edges);
void test_surface_knot_preprocessing(Eigen::MatrixXd &points, Eigen::MatrixXd& knotP, Eigen::MatrixXi& knotE);
void test_knot_fixing(Eigen::MatrixXd &points, Eigen::MatrixXd& knotP, Eigen::MatrixXi& knotE);
void run_ours(const int model, const int nbr_pts, double &per_ours, const std::string path, const std::string tail,
	const double per, const bool enable_local_energy);
void run_Seungyong(const int model, const int nbr_pts, const double tolerance, const std::string path);
void run_piegl(const int model, const int nbr_pts, const double per);
void read_mesh_series(std::string path, std::string namebase, int end);
void run_lofting(const int model, const int nbr_pts, const double per,const std::string path);
void run_mesh_reconstruction(const std::string inpath, const std::string modelname, double &per_ours, const std::string path, const std::string tail,
	const double per, const bool enable_local_energy);
void write_points(const std::string& file, const Eigen::MatrixXd& ver);
void ours_inserting_redundant(const int model, const int nbr_pts, double &per_ours, const std::string path, const std::string tail,
	const double per, const bool enable_insertion);
void run_other_sampling_strategies(const int model, double &per_ours, const std::string path, const std::string tail,
	const double per, const int algorithm);
void initialize_pia_knot_vectors(int degree1, int degree2,
	std::vector<double>& U, std::vector<double>& V, int nu, int nv);
void run_pia(const int model, const int nbr_pts, const int max_itr, const double threadshold, const int cp_nbr_sqrt, const std::string &path);
void run_pia_mesh_reconstruct(const std::string meshfile, const int max_itr, const double threadshold,
	const int cp_nbr_sqrt, const std::string &path);