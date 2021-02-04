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

// build matrix A for surface interpolation (Ax=b)
Eigen::MatrixXd build_matrix_A(const int degree1, const int degree2, 
	const std::vector<double>& U, const std::vector<double>& V,
	const Eigen::MatrixXd& paras) {
	assert(paras.cols() == 2);

	int nu = U.size() - 2 - degree1;// n + 1 = number of u control points
	int nv = V.size() - 2 - degree2;// n + 1 = number of v control points
	int m = paras.size() - 1;// m + 1 = the number of data points
	Eigen::MatrixXd result;
	result.resize(m + 1, (nu + 1)*(nv + 1));
	for (int i = 0; i < result.rows(); i++) {
		double u = paras(i, 0);
		double v = paras(i, 1);
		for (int j = 0; j < result.cols(); j++) {
			// get the indices of N_r(u) and N_q(v)
			int r = j / (nv + 1);
			int q = j - r * (nv + 1);
			double N1 = Nip(r, degree1, u, U);
			double N2 = Nip(q, degree2, v, V);
			result(i, j) = N1 * N2;
		}
	}

	return result;
}
Eigen::MatrixXd build_matrix_A(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V,
	const Eigen::MatrixXd& paras, const std::vector<int> row_id) {
	assert(paras.cols() == 2);

	int nu = U.size() - 2 - degree1;// n + 1 = number of u control points
	int nv = V.size() - 2 - degree2;// n + 1 = number of v control points
	int m = paras.size() - 1;// m + 1 = the number of data points
	Eigen::MatrixXd result;
	result.resize(row_id.size(), (nu + 1)*(nv + 1));
	for (int l = 0; l < row_id.size(); l++) {
		int i = row_id[l];
		double u = paras(i, 0);
		double v = paras(i, 1);
		for (int j = 0; j < result.cols(); j++) {
			// get the indices of N_r(u) and N_q(v)
			int r = j / (nv + 1);
			int q = j - r * (nv + 1);
			double N1 = Nip(r, degree1, u, U);
			double N2 = Nip(q, degree2, v, V);
			result(l, j) = N1 * N2;
		}
	}

	return result;
}

// build vector b for interpolation problem (Ax=b)
Eigen::VectorXd build_Vector_b(const Eigen::MatrixXd& points, const int dimension) {
	Eigen::VectorXd result;
	result.resize(points.rows());
	for (int i = 0; i < result.size(); i++) {
		result(i) = points(i, dimension);
	}
	return result;
}

// build vector b for interpolation problem (Ax=b). Here the returned is from b's start_row row, and has 
// nbr_rows rows.
Eigen::VectorXd build_Vector_b(const Eigen::MatrixXd& points, const int dimension,
	const std::vector<int> row_id) {
	Eigen::VectorXd result;
	result.resize(row_id.size());
	for (int l = 0; l < row_id.size(); l++) {
		int i = row_id[l];
		result(l) = points(i, dimension);
	}
	return result;
}

bool selected_rows_have_solution(const int degree1, const int degree2,
	const std::vector<double>& U, const std::vector<double>& V,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points,
	const std::vector<int> row_id, const int dimension) {

	Eigen::MatrixXd A = build_matrix_A(degree1, degree2, U, V, paras,row_id);
	Eigen::VectorXd b = build_Vector_b(points, dimension, row_id);
	return equation_has_solution(A, b);
}
// UorV shows which interval do we bisect
// TODO now the strategy is always select the left interval as input. maybe we can have a better strategy
void bisect_interval(std::array<double, 2>&Uinterval, std::array<double, 2>&Vinterval,const bool UorV) {
	if (UorV) {
		Uinterval[1] /= 2;
	}
	else {
		Vinterval[1] /= 2;
	}
}

void select_point_id_in_interval(const std::array<double, 2>& Uinterval, const std::array<double, 2>& Vinterval,
	const Eigen::MatrixXd& paras, std::vector<int>&selected) {
	selected.clear();
	for (int i = 0; i < paras.rows(); i++) {
		double u = paras(i, 0);
		double v = paras(i, 1);
		if (u >= Uinterval[0] && u <= Uinterval[1] && v >= Vinterval[0] && v <= Vinterval[1]) {
			selected.push_back(i);
		}
	}
	return;
}

void bisectively_find_solvable_block(const int degree1, const int degree2,
	const std::vector<double>& Uin, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points, const int dimension,
	std::array<double, 2>&Uinterval_out, std::array<double, 2>&Vinterval_out, std::vector<int>& pids) {

	std::vector<int> current_ids(paras.size());
	for (int i = 0; i < current_ids.size(); i++) {
		current_ids[i] = i;// initially all the points are checked
	}

	// the initial block
	std::array<double, 2>Uinterval = { {Uin[0],Uin.back()} };
	std::array<double, 2>Vinterval = { {Vin[0],Vin.back()} };
	
	bool UorV = true;//initially U get bisected
	while (current_ids.size() > 0) {
		if (selected_rows_have_solution(degree1, degree2, Uin, Vin, paras, points, current_ids, dimension)) {
			Uinterval_out = Uinterval;
			Vinterval_out = Vinterval;
			pids = current_ids;
			return;
		}
		
		bisect_interval(Uinterval, Vinterval, UorV);
		UorV = !UorV;

		select_point_id_in_interval(Uinterval, Vinterval, paras, current_ids);
	}
	std::cout << "ERROR OCCURED WHEN TRYING TO LOCATE THE SOLVABLE BLOCK" << std::endl;
	assert(false); //the code cannot go here because there is no solution
	return;
}

void expand_one_point_close_to_interval(
	const std::array<double, 2>& Uinterval,
	const std::array<double, 2>& Vinterval, const Eigen::MatrixXd& paras,
	const std::vector<int> &pid_in, std::vector<int>&pid_out) {
	pid_out = pid_in;

	std::vector<bool> markers(paras.rows(), false);// shows which parameter is already in the block
	for (int i = 0; i < pid_in.size(); i++) {
		markers[pid_in[i]] = true;
	}

	double udis = 1;
	double vdis = 1;
	for (int i = 0; i < paras.rows(); i++) {
		if (markers[i]) continue;

	}

}
void fix_knot_vector_to_interpolate_surface(const int degree1, const int degree2, 
	const std::vector<double>& Uin,const std::vector<double>& Vin,                                                                                           
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points, const int dimension,
	std::vector<double>& Uout, std::vector<double>& Vout
) {
	std::vector<double> Utmp;
	std::vector<double> Vtmp;
	fix_surface_grid_parameter_too_many(degree1, Uin, degree2, Vin, paras, Utmp, Vtmp);

	assert(paras.rows() == points.size());
	std::array<double, 2>Uinterval;
	std::array<double, 2>Vinterval;
	std::vector<int> pids;// point ids of the solvable block

	bisectively_find_solvable_block(degree1, degree2, Uin, Vin, paras, points, dimension, Uinterval, Vinterval, pids);
	
	while (1) {

	}
	
	std::vector<double> expanded_U = init_vec;
	// take init_vec as input, detect if there is any stair whose row is too many (>n+1).
	// if there are such stairs, insert knots to init_vec, the result is expanded_U

	fix_stairs_row_too_many(degree, init_vec, paras, expanded_U);

	assert(points.size() == paras.size());
	int n = expanded_U.size() - 2 - degree;// n + 1 = number of control points
	int m = paras.size() - 1;// m + 1 = the number of data points

	Eigen::MatrixXd A;
	Eigen::VectorXd b;


	A = build_matrix_A(degree, expanded_U, paras);

	b = build_Vector_b(points, dimension);

	bool have_solution = equation_has_solution(A, b);

	int start_row = 0;// initialized row nbr is 0
	int nbr_rows = m + 1;// initialized block has m+1 rows
	int dbg_flag = 0;
	while (!have_solution) {
		if (dbg_flag > 50) exit(0);
		dbg_flag++;

		Eigen::MatrixXd sub_A1, sub_A2;
		Eigen::VectorXd sub_b1, sub_b2;
		int rank_diff1, rank_diff2;
		assert(nbr_rows > 1);

		sub_A1 = build_matrix_A(degree, expanded_U, paras, start_row, nbr_rows - 1);
		sub_b1 = build_Vector_b(points, dimension, start_row, nbr_rows - 1);

		bool have_solution1 = equation_has_solution(sub_A1, sub_b1, rank_diff1);

		sub_A2 = build_matrix_A(degree, expanded_U, paras, start_row + 1, nbr_rows - 1);
		sub_b2 = build_Vector_b(points, dimension, start_row + 1, nbr_rows - 1);

		bool have_solution2 = equation_has_solution(sub_A2, sub_b2, rank_diff2);

		if ((!have_solution1) && (!have_solution2)) {// if no solution, check next level;
			if (rank_diff1 <= rank_diff2) {
				start_row = start_row;
				nbr_rows = nbr_rows - 1;
			}
			else {
				start_row = start_row + 1;
				nbr_rows = nbr_rows - 1;
			}

			continue;
		}
		std::vector<double> tempU;
		if (have_solution1 && (!have_solution2)) {
			// deal with knots, start_row=0, nbr_rows=m+1; calculate have_solution
			// the problematic row is start_row + nbr_rows



			insert_a_knot_to_a_stair(degree, start_row + nbr_rows - 1, expanded_U, paras, tempU, STAIR_FORWARD);

		}
		if ((!have_solution1) && have_solution2) {
			// deal with knots, start_row=0, nbr_rows=m+1; calculate have_solution

			insert_a_knot_to_a_stair(degree, start_row, expanded_U, paras, tempU, STAIR_BACKWARD);

		}
		if (have_solution1 && have_solution2) {
			// deal with knots, start_row=0, nbr_rows=m+1; calculate have_solution
			// if both of them have solution, pick a random one to solve

			insert_a_knot_to_a_stair(degree, start_row + nbr_rows - 1, expanded_U, paras, tempU, STAIR_FORWARD);

		}
		expanded_U = tempU;
		start_row = 0; nbr_rows = m + 1;
		A = build_matrix_A(degree, expanded_U, paras);
		b = build_Vector_b(points, dimension);

		have_solution = equation_has_solution(A, b);

	}
	return expanded_U;

}