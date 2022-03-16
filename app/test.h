#pragma once
#include<sparse_interp/curve.h>


namespace SIBSplines {
	
void run_ours(const int model, const int nbr_pts, double &per_ours, const std::string path, const std::string tail,
	const double per, const bool enable_local_energy);

// multilevel B-spline [Seungyong Lee, et al., 1997, TVCG]
void run_Seungyong(const int model, const int nbr_pts, const double tolerance, const std::string path);

// averaging method to construct the B-spline knot vectors ([Piegl et al. 1996, The NURBS book])
void run_piegl(const int model, const int nbr_pts, const double per);

// this is an assistant function for [Seungyong Lee, et al., 1997, TVCG]
void read_mesh_series(std::string path, std::string namebase, int end);

// run lofting method constructing the knot vectors ([Wen-Ke Wang et al, 2008, CAD])
void run_lofting(const int model, const int nbr_pts, const double per,const std::string path);

// run our method reconstructing triangle mesh
void run_mesh_reconstruction(const std::string inpath, const std::string modelname, double &per_ours, const std::string path, const std::string tail,
	const double per, const bool enable_local_energy);


void write_points(const std::string& file, const Eigen::MatrixXd& ver);

// show local non-smoothness caused by locally too many rendunt control points
void ours_inserting_redundant(const int model, const int nbr_pts, double &per_ours, const std::string path, const std::string tail,
	const double per, const bool enable_insertion);

// show interpolating non-uniformly sampled points
void run_other_sampling_strategies(const int model, double &per_ours, const std::string path, const std::string tail,
	const double per, const int algorithm);

// IGA (Kineri et al. 2012, CAD )
void run_pia(const int model, const int nbr_pts, const int max_itr, const double threadshold, const int cp_nbr_sqrt, const std::string &path);
void run_pia_mesh_reconstruct(const std::string meshfile, const int max_itr, const double threadshold,
	const int cp_nbr_sqrt, const std::string &path);
}