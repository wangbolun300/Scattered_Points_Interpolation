#pragma once
#include<basis.h>
Vector3d BSplineSurfacePoint(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V, const double upara,
	const double vpara, const std::vector<std::vector<Vector3d>>& control);
struct Bsurface {
	int degree1;
	int degree2;
	std::vector<double> U;  
	std::vector<double> V;  
	double upara;
	double vpara; 
	std::vector<std::vector<Vector3d>> control_points;
	int nu();
	int nv();
};
Vector3d BSplineSurfacePoint(const Bsurface& surface, const double upara, const double vpara);