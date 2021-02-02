#include<surface.h>
#include<curve.h>
int Bsurface::nu() {
	return U.size() - 2 - degree1;
}
int Bsurface::nv() {
	return V.size() - 2 - degree2;
}
Vector3d BSplineSurfacePoint(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V, const double upara,
	const double vpara, const std::vector<std::vector<Vector3d>>& control) {
	Eigen::Vector3d result = Eigen::Vector3d(0, 0, 0);
	int nu = U.size() - 2 - degree1;// n + 1 = number of control points
	assert(nu + 1 == control.size());
	int nv = V.size() - 2 - degree2;// n + 1 = number of control points
	assert(nv + 1 == control[0].size());

	for (int i = 0; i < nu + 1; i++) {
		double base1 = Nip(i, degree1, upara, U);
		for (int j = 0; j < nv + 1; j++) {
			double base2 = Nip(j, degree2, vpara, V);
			result += base1 * base2 * control[i][j];
		}
		//std::cout << "base " << base << std::endl;

	}
	return result;
}
Vector3d BSplineSurfacePoint(const Bsurface& surface, const double upara, const double vpara) {
	return BSplineSurfacePoint(surface.degree1, surface.degree2, surface.U, surface.V, upara, vpara, surface.control_points);
}
std::array<double, 2> get_the_mean_value(const Eigen::MatrixXd& paras, const std::vector<int> ids) {
	std::array<double, 2> total = { {0,0} };
	for (int i = 0; i < ids.size(); i++) {
		total[0] = total[0]+paras(ids[i],0);
		total[1] = total[1] + paras(ids[i], 1);
	}
	total[0] /= ids.size();
	total[1] /= ids.size();
	return total;
}
void knot_vector_insert_values(const std::vector<double>& U,const std::vector<double>& V,
	const Eigen::MatrixXd& paras,
	std::vector<std::array<int, 2>> need_fix_intervals, std::vector<std::vector<std::vector<int>>> para_ids,
	std::vector<double>& Uout, std::vector<double>& Vout) {
	// for each interval we need to fix, we insert one value
	std::vector<std::array<double, 2>> insert_values;

	for (int i = 0; i < need_fix_intervals.size(); i++) {
		int uinter = need_fix_intervals[i][0];// the u interval id
		int vinter = need_fix_intervals[i][1];
		insert_values.push_back(get_the_mean_value(paras, para_ids[uinter][vinter]));
	}
	// now need to insert insert_values[i] to U
	std::vector<double> Uresult = U;
	std::vector<double> Vresult = V;

	for (int i = 0; i < insert_values.size(); i++) {
		Uresult = knot_vector_insert_one_value(Uresult, insert_values[i][0]);
		Vresult = knot_vector_insert_one_value(Vresult, insert_values[i][1]);
	}
	assert(Uresult.size() == U.size() + insert_values.size());
	Uout = Uresult;
	Vout = Vresult;
}


// since for a block [u_i, u_(i+1)]x[v_j, v_(j+1)] corresponding to at most 
// (degree1+1)x(degree2+1) control points, the target points in this block should no
// more than (degree1+1)x(degree2+1)
void fix_the_grid_not_border(
	const int degree1, const std::vector<double>& Uin,
	const int degree2, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras,
	std::vector<double>& Uout, std::vector<double>& Vout) {

	int nu = Uin.size() - 2 - degree1;
	int nv = Vin.size() - 2 - degree2;
	assert(paras.cols() == 2);// u v parameters
	Uout.clear();
	Vout.clear();
	const int ps1 = Uin.size() - 1;// the number of u intervals
	const int ps2 = Vin.size() - 1;// the number of v intervals
	std::vector<std::vector<std::vector<int>> > para_ids(ps1);
	std::vector<double> uparas, vparas;
	for (int i = 0; i < para_ids.size(); i++) {
		para_ids[i].resize(ps2);
	}
	int tpush = 0;
	for (int i = 0; i < paras.rows(); i++) {
		bool located = false;
		double u = paras(i, 0), v = paras(i, 1);

		// TODO remind that we don't need to worry about u==1&&v==1 since this point will never cause problem

		for (int j = 0; j < ps1; j++) {
			if (u >= Uin[j] && u < Uin[j + 1]) {
				for (int k = 0; k < ps2; k++) {
					if (v >= Vin[k] && v < Vin[k + 1]) {
						para_ids[j][k].push_back(i);
						tpush++;
						located = true;
						break;
					}
				}

			}
			if (located) break;
		}
	}
	
	// now we know in each interval how many points there are. it should no more than (degree1+1)x(degree2+1)
	std::vector<std::array<int, 2>> need_fix_intervals;
	int at_most = (degree1 + 1)*(degree2 + 1);
	for (int i = 0; i < ps1; i++) {
		for (int j = 0; j < ps2; j++) {
			if (para_ids[i][j].size() > at_most) {
				need_fix_intervals.push_back({ {i,j} });
			}
		}
	}

	if (need_fix_intervals.size() == 0) {
		
		Uout = Uin;
		Vout = Vin;
		return;
	}


	// down here, we need to insert values to U and V
	// for each problematic stair, we insert one value; but this may not be enough,
	// so, we recursive this function

	knot_vector_insert_values(Uin, Vin, paras, need_fix_intervals, para_ids, Uout, Vout);

	std::vector<double> Utmp, Vtmp;
	fix_the_grid_not_border(degree1, Uout, degree2, Vout, paras, Utmp, Vtmp);
	Uout = Utmp;
	Vout = Vtmp;
	return;

}
void fix_surface_grid_parameter_too_many(const int degree1, const std::vector<double>& Uin,
	const int degree2, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras,
	std::vector<double>& Uout, std::vector<double>& Vout) {
	
	// deal with u==1 and v==1. using the curve 
	// TODO also need to deal with u==0 and v==0!
	std::vector<double> uparas, vparas;
	for (int i = 0; i < paras.rows(); i++) {
		double u = paras(i, 0);
		double v = paras(i, 1);
		if (u == Uin.back()) {
			uparas.push_back(u);
		}
		if (v == Vin.back()) {
			vparas.push_back(v);
		}
	}
	std::vector<double> Utmp, Vtmp;

	// when interpolating u==1 and v==1, it will be the same as interpolating
	// two curves. so, we use curve processing method
	fix_stairs_row_too_many(degree1, Uin, uparas, Utmp);
	fix_stairs_row_too_many(degree2, Vin, vparas, Vtmp);

	
	fix_the_grid_not_border(degree1, Utmp, degree2, Vtmp, paras, Uout, Vout);
}


