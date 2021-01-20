#pragma once
#include "Types.hpp"
#include <vector>
namespace workClean {
	namespace border {
		void get_simp_points(const std::vector<Vector3>& border, const int scale, std::vector<int>& pid);
		//get distance of every two neighbouring points,and the total distance 
		void get_border_pp_length(const std::vector<Vector3>&border, std::vector<Scalar>& length, Scalar& dis);

		void get_local_curvature(const std::vector<Vector3>& border, std::vector<Scalar> &cur,
			std::vector<Scalar>& length, Scalar &dis);
	}
	
}