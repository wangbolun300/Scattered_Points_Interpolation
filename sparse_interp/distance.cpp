#include"distance.h"
namespace workClean {


	double two_points_distance(const Vector3& p1, const Vector3 p2) {
		return (p1 - p2).norm();
	}

	// segment s0-s1 
	double point_seg_distance(const Vector3&p, const Vector3 &s0, const Vector3 &s1) {
		if (two_points_distance(s0, s1) < SCALAR_ZERO) return two_points_distance(p, s1);

	}
}