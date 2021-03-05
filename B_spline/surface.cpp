#include<surface.h>
#include<curve.h>
#include <iomanip> 
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
	std::cout << "inside fix grid too many" << std::endl;
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
	int m = paras.rows() - 1;// m + 1 = the number of data points
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
	const std::vector<int> &row_id, const int dimension) {

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

// the block contains all the points that intersected with the interval. it means 
// if we add a new point to the block, the new point will be outside the interval
void bisectively_find_solvable_block(const int degree1, const int degree2,
	const std::vector<double>& Uin, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points, const int dimension,
	std::array<double, 2>&Uinterval_out, std::array<double, 2>&Vinterval_out, std::vector<int>& pids) {

	std::vector<int> current_ids(paras.rows());
	for (int i = 0; i < current_ids.size(); i++) {
		current_ids[i] = i;// initially all the points are checked
	}

	// the initial block
	std::array<double, 2>Uinterval = { {Uin[0],Uin.back()} };
	std::array<double, 2>Vinterval = { {Vin[0],Vin.back()} };
	
	bool UorV = true;//initially U get bisected
	while (current_ids.size() > 0) {
		bool solution = selected_rows_have_solution(degree1, degree2, Uin, Vin, paras, points, current_ids, dimension);
		if (solution) {
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

// the expanded point can be on the border of this interval
// or the outside but close to the interval
void expand_one_point_close_to_interval(
	 const Eigen::MatrixXd& paras,
	const std::vector<int> &pid_in, std::vector<int>&pid_out) {
	
	pid_out = pid_in;

	std::vector<bool> markers(paras.rows(), false);// shows which parameter is already in the block
	for (int i = 0; i < pid_in.size(); i++) {
		markers[pid_in[i]] = true;
	}

	double udis = 1;
	double vdis = 1;
	double dis = 1;
	int tmp = -1;
	for (int i = 0; i < paras.rows(); i++) {
		if (markers[i]) continue;
		double u = paras(i, 0);
		double v = paras(0, 1);
		// in this strategy, the closest point will have smallest u or v
		if (u < dis) {
			tmp = i;
			dis = u;
		}
		if (v < dis) {
			tmp = i;
			dis = v;
		}
	}
	if (tmp > -1) {// it means there are points outside. otherwise it will just return the old list
		pid_out.push_back(tmp);
	}
	else {
		std::cout << "ERROR OCCURED WHEN EXPANDING ONE POINT TO THE BLOCK" << std::endl;
	}

}
void double_the_interval_and_expand_points(
	const std::vector<double> &U, const std::vector<double> &V,
	const std::array<double, 2>&Uinterval, const std::array<double, 2>Vinterval,
	const Eigen::MatrixXd& paras,
	const std::vector<int> &pid_in, std::vector<int>&pid_out,
	std::array<double, 2>& Uinterval_out, std::array<double, 2>& Vinterval_out) {
	pid_out = pid_in;
	std::vector<bool> markers(paras.rows(), false);// shows which parameter is already in the block
	for (int i = 0; i < pid_in.size(); i++) {
		markers[pid_in[i]] = true;
	}
	double ul = Uinterval[0];
	double vl = Vinterval[0];
	double ur = Uinterval[1];
	double vr = Vinterval[1];
	if (ur == U.back() && vr == V.back()) {
		for (int i = 0; i < paras.rows(); i++) {
			if (markers[i]) { continue; }
			pid_out.push_back(i);
		}
		Uinterval_out = { {Uinterval[0],ur} };
		Vinterval_out = { {Vinterval[0],vr} };
		return;
	}
	if (ur <= vr) {
		ur = 2 * ur;
	}
	else {
		vr = 2 * vr;
	}
	for (int i = 0; i < paras.rows(); i++) {
		if (markers[i]) { continue; }
		double u = paras(i, 0);
		double v = paras(i, 1);
		if (u >= ul && u <= ur && v >= vl && v <= vr) {
			pid_out.push_back(i);
		}
		
	}
	Uinterval_out = { {Uinterval[0],ur} };
	Vinterval_out = { {Vinterval[0],vr} };
	return;
}

// 1 insert U, 0 insert V.
// TODO maybe have a better way for udiff>0&&vdiff>0
bool insert_U_or_V_direction(const double udiff, const double vdiff) {
	assert(udiff >= 0 || vdiff >= 0);// (udiff<0&&vdiff<0) should not happen
	if (udiff ==0) return 0;
	if (vdiff == 0) return 1;
	if (udiff > 0) return 1;
	if (vdiff > 0) return 0;
	std::cout << "ERROR impossible case when inserting u or v" << std::endl;
	assert(false);
	return false;
}

double get_mean_value(const Eigen::MatrixXd& paras, const std::vector<int> &ids, const bool insert_U) {
	int tuv = insert_U ? 0 : 1;
	double total = 0;
	for (int i = 0; i < ids.size(); i++) {
		int ID = ids[i];
		total += paras(ID, tuv);
	}
	return total / ids.size();
}

// given a id list pids, and pids.back() is the problematic one. return an id list in which the points 
// are in the same U-V interval with the problematic one
std::vector<int> selected_ids_in_this_block(const std::vector<double>& Uin, const std::vector<double>& Vin, 
	const Eigen::MatrixXd& paras, const std::vector<int> &pids_in) {
	int id = pids_in.back();
	double u = paras(id, 0);
	double v = paras(id, 1);
	std::vector<double> knot_vector;
	double uv;
	int tuv;
	std::vector<int> related_id;
	int interval = -1;
	
	std::vector<int> pids = pids_in;
	int counter = 0;
	while (counter < 2) {
		
		if (counter==0) {
			knot_vector = Uin;
			uv = u;
			tuv = 0;
		}
		else {
			knot_vector = Vin;
			uv = v;
			tuv = 1;
		}

		if (uv != knot_vector.back()) {
			for (int i = 0; i < knot_vector.size() - 1; i++) {
				if (uv >= knot_vector[i] && uv < knot_vector[i + 1]) {
					interval = i;
					break;
				}
			}
			for (int i = 0; i < pids.size(); i++) {
				int tmp_id = pids[i];
				double tmp_u = paras(tmp_id, tuv);
				if (tmp_u >= knot_vector[interval] && tmp_u < knot_vector[interval + 1]) {
					related_id.push_back(pids[i]);
				}
			}
		}
		else {
			for (int i = 0; i < pids.size(); i++) {
				int tmp_id = pids[i];
				double tmp_u = paras(tmp_id, tuv);
				if (tmp_u == knot_vector.back()) {
					related_id.push_back(pids[i]);
				}
			}
		}
		pids = related_id;
		related_id.clear();
		counter++;
	}
	return pids;
}

// pids.back() is the problematic one
void gather_points_and_get_mean_value(const std::vector<double>& Uin, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras,
	const bool insert_U, const std::vector<int> &pids,double& mean) {
	std::vector<int> related_id = selected_ids_in_this_block(Uin, Vin, paras, pids);
	// get the mean value of the related points
	mean = get_mean_value(paras, related_id, insert_U);

}

bool gather_insert_and_check_solvable(const int degree1, const int degree2, 
	const std::vector<double>& Uin, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras,
	const bool insert_U, const std::vector<int> &pids, const Eigen::MatrixXd& points,
	const int dimension, std::vector<double>&Uout, std::vector<double>&Vout) {
	double value;
	gather_points_and_get_mean_value(Uin, Vin, paras, insert_U, pids, value);
	if (insert_U) {// insert U
		Uout = knot_vector_insert_one_value(Uin, value);
		Vout = Vin;
	}
	else {// insert V
		Uout = Uin;
		Vout = knot_vector_insert_one_value(Vin, value);
	}
	bool solvable = selected_rows_have_solution(degree1, degree2, Uout, Vout, paras, points, pids, dimension);
	return solvable;
}

// to make the problematic block solvable
// pids.back() is the problematic one
void insert_a_knot_to_problematic_area(const int degree1, const int degree2,
	const std::vector<double>& Uin, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points, const int dimension,
	const std::vector<int> &pids, const std::array<double, 2> &Uinterval, const std::array<double, 2>&Vinterval,
	std::vector<double>& Uout, std::vector<double>& Vout) {
	int id = pids.back();
	double u = paras(id, 0);
	double v = paras(id, 1);
	double udiff = u - Uinterval[1];
	double vdiff = v - Vinterval[1];

	bool insert_which = insert_U_or_V_direction(udiff, vdiff);
	
	double value;
	std::vector<double> Utmp;
	std::vector<double> Vtmp;
	
	// first try 
	bool solvable = gather_insert_and_check_solvable(degree1, degree2, Uin, Vin, paras, insert_which, pids,
		points, dimension, Uout, Vout);
	if (solvable) {
		return;
	}
	else {
		if (insert_which) {
			Utmp = Uout;
		}
		else {
			Vtmp = Vout;
		}
	}

	// second try
	solvable= gather_insert_and_check_solvable(degree1, degree2, Uin, Vin, paras, !insert_which, pids,
		points, dimension, Uout, Vout);
	if (solvable) {
		return;
	}
	else {
		if (!insert_which) {
			Utmp = Uout;
		}
		else {
			Vtmp = Vout;
		}
	}

	// the last try, both U and V get inserted
	solvable= selected_rows_have_solution(degree1, degree2, Utmp, Vtmp, paras, points, pids, dimension);
	if (solvable) {
		Uout = Utmp;
		Vout = Vtmp;
		return;
	}
	else {
		std::cout << "u v diff, " << udiff << " " << vdiff << std::endl;
		std::cout << "Utmp, Vtmp\n";
		print_vector(Utmp);
		print_vector(Vtmp);
		std::cout << "current list length " << pids.size() << std::endl;
		std::cout << "original knot, \n";
		print_vector(Uin);
		print_vector(Vin);
		std::cout << "u and v " << u << " " << v << std::endl;
		std::cout << "ERROR OCCURES HERE, INSERTING KNOT DIDNOT WORK" << std::endl;
		//insert_a_knot_to_problematic_area(degree1, degree2, Utmp, Vtmp, paras, points, dimension, pids, Uinterval, Vinterval,
		//	Uout, Vout);// recursive cannot fix it. i think it is because of numerical problem
		exit(0);
	}
	

}

void DBG_insert_a_knot_to_problematic_area(const int degree1, const int degree2,
	const std::vector<double>& Uin, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points, const int dimension,
	const std::vector<int> &pids, const std::array<double, 2> &Uinterval, const std::array<double, 2>&Vinterval,
	std::vector<double>& Uout, std::vector<double>& Vout) {
	int id = pids.back();
	double u = paras(id, 0);
	double v = paras(id, 1);
	double udiff = u - Uinterval[1];
	double vdiff = v - Vinterval[1];

	bool insert_which = insert_U_or_V_direction(udiff, vdiff);

	double value;
	std::vector<double> Utmp;
	std::vector<double> Vtmp;

	// first try 
	bool solvable = gather_insert_and_check_solvable(degree1, degree2, Uin, Vin, paras, insert_which, pids,
		points, dimension, Uout, Vout);
	if (solvable) {
		return;
	}
	else {
		if (insert_which) {
			Utmp = Uout;
		}
		else {
			Vtmp = Vout;
		}
	}

	// second try
	solvable = gather_insert_and_check_solvable(degree1, degree2, Uin, Vin, paras, !insert_which, pids,
		points, dimension, Uout, Vout);
	if (solvable) {
		return;
	}
	else {
		if (!insert_which) {
			Utmp = Uout;
		}
		else {
			Vtmp = Vout;
		}
	}

	// the last try, both U and V get inserted
	solvable = selected_rows_have_solution(degree1, degree2, Utmp, Vtmp, paras, points, pids, dimension);
	if (solvable) {
		Uout = Utmp;
		Vout = Vtmp;
		return;
	}
	else {
		std::cout << "u v diff, " << udiff << " " << vdiff << std::endl;
		std::cout << "Utmp, Vtmp\n";
		print_vector(Utmp);
		print_vector(Vtmp);
		std::cout << "current list length " << pids.size() << std::endl;
		std::cout << "original knot, \n";
		print_vector(Uin);
		print_vector(Vin);
		std::cout << "u and v " << u << " " << v << std::endl;
		std::cout << "ERROR OCCURES HERE, INSERTING KNOT DIDNOT WORK" << std::endl;
		//insert_a_knot_to_problematic_area(degree1, degree2, Utmp, Vtmp, paras, points, dimension, pids, Uinterval, Vinterval,
		//	Uout, Vout);// recursive cannot fix it. i think it is because of numerical problem
		Uout = Utmp;
		Vout = Vtmp;
		return;
	}


}

// construct the boolean matrix/net to show the key points to be interpolated
// TODO
void construct_boolean_matrix(const std::vector<double>& U, const std::vector<double>& V,
	const Eigen::MatrixXd& paras, Eigen::MatrixXi& bm) {

}
//void fix_knot_vector_to_interpolate_surface_boolean(const int degree1, const int degree2,
//	const std::vector<double>& Uin, const std::vector<double>& Vin,
//	const Eigen::MatrixXd& paras, 
//	std::vector<double>& Uout, std::vector<double>& Vout
//) {
//	std::vector<double> Utmp;
//	std::vector<double> Vtmp;
//	fix_surface_grid_parameter_too_many(degree1, Uin, degree2, Vin, paras, Utmp, Vtmp);
//	std::cout << "surface grid too many get fixed" << std::endl;
//	assert(paras.rows() == points.rows());
//	std::array<double, 2>Uinterval;
//	std::array<double, 2>Vinterval;
//	std::vector<int> pids;// point ids of the solvable block
//	
//    // construct 
//
//
//
//	bisectively_find_solvable_block(degree1, degree2, Utmp, Vtmp, paras, points, dimension, Uinterval, Vinterval, pids);
//	//std::cout << "paras\n" << paras << std::endl;
//	std::cout << "found solvable block" << std::endl;
//	if (pids.size() == points.rows()) {// it means all the points are solvable
//		Uout = Utmp;
//		Vout = Vtmp;
//		std::cout << "return1 pids" << std::endl;
//
//		return;
//	}
//
//	// if goes here, it means the function is not solvable.
//	while (1) {
//		std::vector<int> new_ids;
//
//		// expand a new point to make the block larger
//		// but remind that the intervals are not updated here
//		expand_one_point_close_to_interval(paras, pids, new_ids);
//
//		// and test if the current block is solvable. if is, continue to expand another point;
//		// if isn't, gather the problematic points and insert a knot
//		bool solvable = selected_rows_have_solution(degree1, degree2, Utmp, Vtmp, paras, points, new_ids, dimension);
//		if (solvable) {
//
//			if (new_ids.size() == points.rows()) {// it means all the points are solvable
//				Uout = Utmp;
//				Vout = Utmp;
//				std::cout << "return2 pids" << std::endl;
//
//				return;
//			}
//
//		}
//		else {// if the new inserted point break the solvability, then insert a point
//
//			std::vector<double> Utmpout, Vtmpout;
//			std::vector<int> pid_out;
//			std::array<double, 2> Uinterval_out, Vinterval_out;
//			insert_a_knot_to_problematic_area(degree1, degree2, Utmp, Vtmp, paras, points, dimension, new_ids,
//				Uinterval, Vinterval, Utmpout, Vtmpout);
//
//			// TODO there is a smarter way to select double u or v
//			double_the_interval_and_expand_points(Utmpout, Vtmpout, Uinterval, Vinterval, paras,
//				new_ids, pid_out, Uinterval_out, Vinterval_out);
//
//			solvable = selected_rows_have_solution(degree1, degree2, Utmpout, Vtmpout, paras, points, pid_out, dimension);
//			if (solvable) {// if it is solvable, it means we can directly replace the list with the expanded list
//				if (pid_out.size() == points.rows()) {// it means all the points are solvable
//					Uout = Utmpout;
//					Vout = Vtmpout;
//					std::cout << "return3 pids" << std::endl;
//
//					return;
//				}
//
//				// if it is solvable, update the list and intervals
//				Uinterval = Uinterval_out;
//				Vinterval = Vinterval_out;
//				new_ids = pid_out;
//			}
//			// if the doubled interval is useless, do not update intervals and ids
//			Utmp = Utmpout;
//			Vtmp = Vtmpout;
//		}
//
//		pids = new_ids;
//	}
//
//}
void fix_knot_vector_to_interpolate_surface(const int degree1, const int degree2, 
	const std::vector<double>& Uin,const std::vector<double>& Vin,                                                                                           
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points, const int dimension,
	std::vector<double>& Uout, std::vector<double>& Vout
) {
	std::vector<double> Utmp;
	std::vector<double> Vtmp;
	fix_surface_grid_parameter_too_many(degree1, Uin, degree2, Vin, paras, Utmp, Vtmp);
	std::cout << "surface grid too many get fixed" << std::endl;
	assert(paras.rows() == points.rows());
	std::array<double, 2>Uinterval;
	std::array<double, 2>Vinterval;
	std::vector<int> pids;// point ids of the solvable block

	bisectively_find_solvable_block(degree1, degree2, Utmp, Vtmp, paras, points, dimension, Uinterval, Vinterval, pids);
	//std::cout << "paras\n" << paras << std::endl;
	std::cout << "found solvable block" << std::endl;
	if (pids.size() == points.rows()) {// it means all the points are solvable
		Uout = Utmp;
		Vout = Vtmp;
		std::cout << "return1 pids" << std::endl;
		
		return;
	}

	// if goes here, it means the function is not solvable.
	while (1) {
		std::vector<int> new_ids;

		// expand a new point to make the block larger
		// but remind that the intervals are not updated here
		expand_one_point_close_to_interval(paras, pids, new_ids);

		// and test if the current block is solvable. if is, continue to expand another point;
		// if isn't, gather the problematic points and insert a knot
		bool solvable = 
#ifdef SPARSE_INTERP_WITH_GMP
			selected_rows_have_solution_rational(degree1, degree2, Utmp, Vtmp, paras, points, new_ids, dimension);
#else
			selected_rows_have_solution(degree1, degree2, Utmp, Vtmp, paras, points, new_ids, dimension);
#endif
		if (solvable) {
		
			if (new_ids.size() == points.rows()) {// it means all the points are solvable
				Uout = Utmp;
				Vout = Utmp;
				std::cout << "return2 pids" << std::endl;
				
				return;
			}
			
		}
		else {// if the new inserted point break the solvability, then insert a point
			
			std::vector<double> Utmpout, Vtmpout;
			std::vector<int> pid_out;
			std::array<double, 2> Uinterval_out, Vinterval_out;
			insert_a_knot_to_problematic_area(degree1, degree2, Utmp, Vtmp, paras, points, dimension, new_ids,
				Uinterval, Vinterval, Utmpout, Vtmpout);
			
			// TODO there is a smarter way to select double u or v
			double_the_interval_and_expand_points(Utmpout, Vtmpout,Uinterval,Vinterval,paras,
				new_ids,pid_out, Uinterval_out,Vinterval_out);

			solvable = 
#ifdef SPARSE_INTERP_WITH_GMP
				selected_rows_have_solution_rational(degree1, degree2, Utmpout, Vtmpout, paras, points, pid_out, dimension);
#else
				selected_rows_have_solution(degree1, degree2, Utmpout, Vtmpout, paras, points, pid_out, dimension);
#endif
			if (solvable) {// if it is solvable, it means we can directly replace the list with the expanded list
				if (pid_out.size() == points.rows()) {// it means all the points are solvable
					Uout = Utmpout;
					Vout = Vtmpout;
					std::cout << "return3 pids" << std::endl;
				
					return;
				}

				// if it is solvable, update the list and intervals
				Uinterval = Uinterval_out;
				Vinterval = Vinterval_out;
				new_ids = pid_out;
			}
			// if the doubled interval is useless, do not update intervals and ids
			Utmp = Utmpout;
			Vtmp = Vtmpout;
		}

		pids = new_ids;
	}

}

void fix_knot_vector_to_interpolate_surface(const int degree1, const int degree2,
	const std::vector<double>& Uin, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points,
	std::vector<double>& Uout, std::vector<double>& Vout) {
	std::vector<double> Utmp = Uin, Vtmp = Vin;
	for (int i = 0; i < 3; i++) {
		fix_knot_vector_to_interpolate_surface(degree1, degree2, Utmp, Vtmp, paras, points, i, Uout, Vout);
		Utmp = Uout;
		Vtmp = Vout;
	}
	
}

// theoretically this should interpolate any surface. but by testing we found that this can have numerical problems.
// so, we need boolean operations
void easist_way_to_fix_knot_vector_to_interpolate_surface(const int degree1, const int degree2,
	const std::vector<double>& Uin, const std::vector<double>& Vin,
	const Eigen::MatrixXd& paras, const Eigen::MatrixXd& points, 
	std::vector<double>& Uout, std::vector<double>& Vout) {
	Uout = Uin;
	Vout = Vin;
	//std::cout << "mark1" << std::endl;
	for (int i = 0; i < paras.rows(); i++) {
		double u = paras(i, 0);
		double v = paras(i, 1);
		bool uskip = false, vskip = false;
		for (int j = 0; j < Uout.size(); j++) {
			if (u == Uout[j]) {
				uskip = true;
				break;
			}
		}
		for (int j = 0; j < Vout.size(); j++) {
			if (v == Vout[j]) {
				vskip = true;
				break;
			}
		}
		//std::cout << "i"<<i << std::endl;
		if (!uskip) {
			if (u == 1) {
				std::cout << "error should not happen" << std::endl;
			}

			Uout = knot_vector_insert_one_value(Uout, u);
		}
		if (!vskip) {
			Vout = knot_vector_insert_one_value(Vout, v);
		}
		//std::cout << "i pushed" << i << std::endl;
		
	}
	//std::cout << "mark2" << std::endl;
	std::cout <<std::setprecision(17)<< "vxx " << Vout[Vout.size() - 5]<<" "<< Vout[Vout.size() - 4]<<" "
		<< Vout[Vout.size() - 3]<<" "<< Vout[Vout.size() - 2]<<" "<< Vout[Vout.size() - 1] << std::endl;
	bool eq = (Vout[Vout.size() - 5] == 1);
	std::cout << "eq " << eq << std::endl;
}
