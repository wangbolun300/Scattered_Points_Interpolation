#pragma once
#include <sparse_interp/surface.h>
namespace SIBSplines
{
	// this file contains the comparison methods:
	// 1. lofting (with knot vector generation method from [Wen-Ke Wang et al, 2008, CAD])
	// 2. multilevel B-splines ([Seungyong Lee, et al., 1997, TVCG])
	// 3. Averaging ([Piegl et al. 1996, The NURBS book])
	// 4. IGA (Kineri et al. 2012, CAD )
	namespace comparison
	{

		// Averaging
		void piegl_method_generate_interpolation_knot_vectors(int degree1, int degree2,
															  std::vector<double> &Uknot, std::vector<double> &Vknot,
															  const Eigen::MatrixXd &param_original,
															  const double per);

		// PIA
		void progressive_iterative_approximation(Bsurface &surface, const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver,
												 const int max_itr, const double threadshold);
		void initialize_pia_knot_vectors(int degree1, int degree2,
										 std::vector<double> &U, std::vector<double> &V, int nu, int nv);

		// lofting
		void lofting_method_generate_interpolation_knot_vectors(const bool start_from_v_direction, int degree1, int degree2,
																std::vector<double> &Uknot, std::vector<double> &Vknot,
																const Eigen::MatrixXd &param_original,
																const double per);

		// multilevel B-splines
		void iteratively_approximate_method(int degree1, int degree2,
											std::vector<double> &Uknot, std::vector<double> &Vknot,
											const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver,
											const double tolerance,
											std::vector<Bsurface> &surfaces, const double per);
	}
}