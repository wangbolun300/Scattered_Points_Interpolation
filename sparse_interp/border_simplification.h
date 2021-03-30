#pragma once
#include <Types.hpp>
#include <vector>

	namespace border {
		// remove the points too close to each other
		std::vector<Vector3d> border_remove_redundant_points(const std::vector<Vector3d>&border);
		
			// get simplified points of the border. d1 is the expected distance of two points, d2 is the maximum search distance
	// constrain 1, the simplified points not too close to each other, constrain 2, the simplified shape is 
	// within an epsilon envelope.
		void get_simp_points(const std::vector<Vector3d>& border, const double & d1, const double d2,
			const double tolerance, std::vector<int>& pid);
		// use expected point number to decide the sample length
		void get_simp_points(const std::vector<Vector3d>& border, const int expect_p_nbr,
			const double tolerance, std::vector<int>& pid);

	}
	
