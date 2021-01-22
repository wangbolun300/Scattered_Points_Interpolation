#pragma once
#include<Eigen/Core>
#include"basis.h"
#include"Types.hpp"
Vector3d BsplinePoint(const int degree, const std::vector<double>& U, const double para,
	const std::vector<Vector3d> pts);