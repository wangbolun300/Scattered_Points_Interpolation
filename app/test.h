#pragma once
#include"curve.h"


namespace SIBSplines {
	
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

void run_pia(const int model, const int nbr_pts, const int max_itr, const double threadshold, const int cp_nbr_sqrt, const std::string &path);
void run_pia_mesh_reconstruct(const std::string meshfile, const int max_itr, const double threadshold,
	const int cp_nbr_sqrt, const std::string &path);
}