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
	for (int i = 0; i < paras.rows(); i++) {
		bool located = false;
		double u = paras(i, 0), v = paras(i, 1);

		// u==1&&v!=1
		if (u == Uin[Uin.size() - 1] && v != Vin[Vin.size() - 1]) {
			for (int j = 0; j < ps2; j++) {
				if (v >= Vin[j] && v < Vin[j + 1]) {
					para_ids[ps1][j].push_back(i);
					located = true;
					break;
				}
			}
		}

		// u!=1&&v==1
		if (u != Uin[Uin.size() - 1] && v == Vin[Vin.size() - 1]) {
			for (int j = 0; j < ps1; j++) {
				if (u >= Uin[j] && u < Uin[j + 1]) {
					para_ids[j][ps2].push_back(i);
					located = true;
					break;
				}
			}
		}

		// TODO remind that we don't need to worry about u==1&&v==1 since this point will never cause problem
		if (located) continue;




		for (int j = 0; j < ps1; j++) {
			if (u >= Uin[j] && u < Uin[j + 1]) {
				for (int k = 0; k < ps2; k++) {
					if (v >= Vin[k] && v < Vin[k + 1]) {
						para_ids[j][k].push_back(i);
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
			if (para_ids.size() > at_most) {
				need_fix_intervals.push_back({ {i,j} });
			}
		}
	}
	if (need_fix_intervals.size() == 0) {
		// TODO we need to consider about u==1 and v==1
		Uout = Uin;
		Vout = Vin;
		return;
	}


	// down here, we need to insert values to U and V
	// for each problematic stair, we insert one value; but this may not be enough,
	// so, we recursive this function

	Uout = knot_vector_insert_values(Uin, paras, need_fix_intervals, para_ids);

	std::vector<double> Utmp;
	fix_stairs_row_too_many(degree, Uout, paras, Utmp);
	Uout = Utmp;
	return;

}
void fix_surface_grid_parameter_too_many(const int degree1, const std::vector<double>& Uin,
	const int degree2, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras,
	std::vector<double>& Uout, std::vector<double>& Vout) {
	
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
	fix_stairs_row_too_many(degree1, Uin, uparas, Utmp);
	fix_stairs_row_too_many(degree2, Vin, vparas, Vtmp);

	xxxx
		fix_the_grid_not_border()
}


