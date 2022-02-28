#include"curve.h"
#include <Eigen/Dense>
#include<iostream>
#include<queue>
int Bcurve::nu() {
	return U.size() - 2 - degree;
}
// return a point position of a given curve
Vector3d BsplinePoint(const int degree, const std::vector<double>& U, const double para,
	const std::vector<Vector3d> &pts) {
	Eigen::Vector3d result = Eigen::Vector3d(0, 0, 0);
	for (int i = 0; i < pts.size(); i++) {
		double base = Nip(i, degree, para, U);
		//std::cout << "base " << base << std::endl;
		result += base * pts[i];
	}
	return result;
}
// return a point position of a given curve
Vector3d BsplinePoint(const int degree, const std::vector<double>& U, const double para,
	const Eigen::MatrixXd& pts) {
	Eigen::Vector3d result = Eigen::Vector3d(0, 0, 0);
	for (int i = 0; i < pts.rows(); i++) {
		double base = Nip(i, degree, para, U);
		//std::cout << "base " << base << std::endl;
		result += base * pts.row(i);
	}
	return result;
}

//TODO this error is not what we want
Eigen::MatrixXd slove_linear_system(const Eigen::MatrixXd& A, const Eigen::MatrixXd &b,
	const bool check_error, double &relative_error) {
	Eigen::MatrixXd x = A.colPivHouseholderQr().solve(b);
	//Eigen::VectorXd x = A.fullPivLu().solve(b);

	if (check_error) {
		relative_error = (A*x - b).norm() / b.norm();
	}
	return x;
}
Eigen::VectorXd slove_linear_system(const Eigen::MatrixXd& A, const Eigen::VectorXd &b,
	const bool check_error, double &relative_error) {
	Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
	//Eigen::VectorXd x = A.fullPivLu().solve(b);

	if (check_error) {
		relative_error = (A*x - b).norm() / b.norm();
	}
	return x;
}
Eigen::MatrixXd build_matrix_N(const int degree, const std::vector<double>& U,
	const std::vector<double>& paras) {
	int n = U.size() - 2 - degree; // there are n+1 control points;
	int m = paras.size() - 1;  // there are m+1 points to fit
	Eigen::MatrixXd result(m - 1, n - 1);
	for (int i = 0; i < m - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			result(i, j) = Nip(j + 1, degree, paras[i + 1], U);
		}
	}
	return result;
}

// use NTN*D = R to solve the n-1 control points D
Eigen::MatrixXd build_matrix_R(const int degree, const std::vector<double>& U,
	const std::vector<double>& paras, const std::vector<Vector3d>& points) {
	int n = U.size() - 2 - degree; // there are n+1 control points;
	int m = paras.size() - 1;  // there are m+1 points to fit
	Eigen::MatrixXd result(n - 1, 3);
	std::vector<Vector3d> ri(m);
	for (int i = 1; i < m; i++) {// from 1 to m-1
		ri[i] =
			points[i] - points[0] * Nip(0, degree, paras[i], U) - points[m] * Nip(n, degree, paras[i], U);
	}

	for (int i = 0; i < n - 1; i++) {
		Vector3d row(0, 0, 0);
		for (int j = 1; j < m; j++) { // j from 1 to m-1
			row += Nip(i + 1, degree, paras[j], U)*ri[j];
		}
		result.row(i) = row;
	}
	return result;
}
// build matrix A for curve interpolation (Ax=b)
Eigen::MatrixXd build_matrix_A(const int degree, const std::vector<double>& U,
	const std::vector<double>& paras) {
	int n = U.size() - 2 - degree;// n + 1 = number of control points
	int m = paras.size() - 1;// m + 1 = the number of data points
	
	Eigen::MatrixXd result(m + 1, n + 1);
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < m + 1; j++) {
			result(j, i) = Nip(i, degree, paras[j], U);
		}
	}
	return result;
}

// build matrix A for curve interpolation (Ax=b), here the returned matrix is the block
// which start from start_row,and have nbr_rows of rows of A
Eigen::MatrixXd build_matrix_A(const int degree, const std::vector<double>& U,
	const std::vector<double>& paras, const int start_row, const int nbr_rows) {
	int n = U.size() - 2 - degree;// n + 1 = number of control points
	int m = paras.size() - 1;// m + 1 = the number of data points
	Eigen::MatrixXd result(nbr_rows, n + 1);
	for (int i = start_row; i < start_row + nbr_rows; i++) {
		for (int j = 0; j < n + 1; j++) {
			result(i - start_row, j) = Nip(j, degree, paras[i], U);
		}
	}

	return result;
}
// build vector b for curve interpolation (Ax=b)
Eigen::VectorXd build_Vector_b(const std::vector<Vector3d>& points, const int dimension) {
	Eigen::VectorXd result;
	result.resize(points.size());
	for (int i = 0; i < result.size(); i++) {
		result(i) = points[i][dimension];
	}
	return result;
}

// build vector b for curve interpolation (Ax=b). Here the returned is from b's start_row row, and has 
// nbr_rows rows.
Eigen::VectorXd build_Vector_b(const std::vector<Vector3d>& points, const int dimension,
	const int start_row, const int nbr_rows) {
	Eigen::VectorXd result;
	result.resize(nbr_rows);
	for (int i = start_row; i < start_row + nbr_rows; i++) {
		result(i - start_row) = points[i][dimension];
	}
	return result;
}
Eigen::MatrixXd solve_curve_control_points(const int degree, const std::vector<double>& U,
	const std::vector<double>& paras, const std::vector<Vector3d>& points) {
	int npoints = points.size();
	int n = U.size() - 2 - degree; // there are n+1 control points;
	Eigen::MatrixXd result(n + 1, 3);
	assert(npoints == paras.size());
	Eigen::MatrixXd N = build_matrix_N(degree, U, paras);
	Eigen::MatrixXd NTN = N.transpose()*N;
	Eigen::MatrixXd R = build_matrix_R(degree, U, paras, points);// NTN*D=R

	bool check_error = true;
	double error;
	Eigen::MatrixXd interior = slove_linear_system(NTN, R, check_error, error);
	std::cout << "error, " << error << std::endl;

	result.row(0) = points[0];
	result.row(n) = points[npoints - 1];
	result.middleRows(1, n - 1) = interior;
	return result;
}
int rank(const Eigen::MatrixXd& matrix) {
	return matrix.completeOrthogonalDecomposition().rank();
}

// check if there is solution(s) for interpolation problem
bool equation_has_solution(const Eigen::MatrixXd& A,
	const Eigen::VectorXd& b) {
	Eigen::MatrixXd  Ab;
	int rankA = rank(A);
	Ab.resize(A.rows(), A.cols() + 1);
	Ab << A, b;
	int rankAb = rank(Ab);
	//std::cout << "Ab,\n" << Ab << std::endl;
	if (rankA == rankAb) {
		return true;
	}
	return false;
}
// rank_diff is the difference between rank(A) and rank(Ab)
bool equation_has_solution(const Eigen::MatrixXd& A,
	const Eigen::VectorXd& b, int& rank_diff) {
	rank_diff = 0;
	Eigen::MatrixXd  Ab;
	int rankA = rank(A);
	Ab.resize(A.rows(), A.cols() + 1);
	Ab << A, b;
	int rankAb = rank(Ab);
	
	if (rankA == rankAb) {
		return true;
	}
	rank_diff = abs(rankAb - rankA);
	return false;
}
// get the mean value of paras[ids[0]],paras[ids[1]], ...
double get_the_mean_value(const std::vector<double>&paras, const std::vector<int> ids) {
	double total = 0;
	for (int i = 0; i < ids.size(); i++) {
		total += paras[ids[i]];
	}
	return total / ids.size();
}
std::vector<double> knot_vector_insert_one_value(const std::vector<double>& U, const double value) {
	std::vector<double> result;
	//result.reserve(U.size() + 1);

	for (int i = 0; i < U.size(); i++) {
		result.push_back(U[i]);
		
		if (i < U.size() - 1) {
			if (value >= U[i] && value < U[i + 1]) {
				assert(value != U[i]);
				if (value == U[i]) {
					continue;
				}
				result.push_back(value);

			}
		}
	}
	
	return result;
}

std::vector<double> knot_vector_insert_values(const std::vector<double>& U, const std::vector<double>& paras,
	std::vector<int> need_fix_intervals, std::vector<std::vector<int>> para_ids) {
	// for each interval we need to fix, we insert one value
	std::vector<double> insert_values;
	int psize = paras.size();
	std::cout << "before first for, need_fix_intervals.size() "<< need_fix_intervals.size() << std::endl;
	for (int i = 0; i < need_fix_intervals.size(); i++) {
		std::cout << i << std::endl;
		std::cout << "this value" << need_fix_intervals[i] << std::endl;
		double insert = need_fix_intervals[i] >= 0 ? get_the_mean_value(paras, para_ids[need_fix_intervals[i]])
			: 0.5*(paras[psize - 1] + paras[psize - 2]);
		insert_values.push_back(insert);
	}
	std::cout << "after first for" << std::endl;
	// now need to insert insert_values[i] to U
	std::vector<double> result = U;
	std::cout << "before second for" << std::endl;
	for (int i = 0; i < insert_values.size(); i++) {
		result = knot_vector_insert_one_value(result, insert_values[i]);
	}
	std::cout << "before second for" << std::endl;
	assert(result.size() == U.size() + insert_values.size());
	return result;
}

// check if there is a stair, the number of rows of the stair > n+1. if happens, fix it
void fix_stairs_row_too_many(const int degree, const std::vector<double>& Uin,
	const std::vector<double>& paras, std::vector<double>& Uout) {
	// Uin has n+degree+2 values, there are n-degree+2 different values, n-degree+1 intervals
	int n = Uin.size() - 2 - degree;
	Uout.clear();
	//std::vector<int> multiplicity(Uin.size() - 1);

	//TODO remind that we don't need to deal with u==1, since it will never cause problem
	std::vector<std::vector<int>> para_ids(Uin.size() - 1);
	for (int i = 0; i < paras.size(); i++) {
		for (int j = 0; j < Uin.size() - 1; j++) {
			if (paras[i] >= Uin[j] && paras[i] < Uin[j + 1]) {
				para_ids[j].push_back(i);
				break;
			}
		}
	}
	std::vector<int> need_fix_intervals;
	// check multiplicity. if larger than n+1, then need to fix
	for (int j = 0; j < Uin.size() - 1; j++) {
		if (para_ids[j].size() > n + 1) {
			need_fix_intervals.push_back(j);
		}
	}
	if (need_fix_intervals.size() == 0) {
		Uout = Uin;
		return;
	}
	// down here, we need to insert values to U
	// for each problematic stair, we insert one value; but this may not be enough,
	// so, we recursive this function
	
	Uout = knot_vector_insert_values(Uin, paras, need_fix_intervals, para_ids);

	std::vector<double> Utmp;
	fix_stairs_row_too_many(degree, Uout, paras, Utmp);
	Uout = Utmp;
	return;

}

// given a problematic row number (this row must be on a stair whose multiplicity is larger than 1, or maybe it 
// is the paras.back())
// and insert one knot to this stair
// since one block (sub_A*x=sub_b) may split one stair into two, construct_method can select from
// STAIR_FORWARD or STAIR_BACKWARD to select certain rows that count into mean value calculation
void insert_a_knot_to_a_stair(const int degree, const int pro_row_id, const std::vector<double>& Uin,
	const std::vector<double>& paras, std::vector<double>& Uout, const int construct_method= STAIR_WHOLE) {

	// Uin has n+degree+2 values, there are n-degree+2 different values, n-degree+1 intervals
	int n = Uin.size() - 2 - degree;

	Uout.clear();

	/*std::vector<int> multiplicity(Uin.size() - 1);*/
	std::vector<std::vector<int>> para_ids(Uin.size() - 1);

	int which_interval = -1; // gives which interval does pro_row_id corresponding to
	bool j_pushed = false;
	for (int i = 0; i < paras.size(); i++) {
		for (int j = 0; j < Uin.size() - 1; j++) {
			if (paras[i] >= Uin[j] && paras[i] < Uin[j + 1]) {
				
				if (i == pro_row_id) {
					which_interval = j;
					if (construct_method == STAIR_BACKWARD) {
						para_ids[j].clear();// if check backward, clear the forward part
					}
				}
				if (construct_method == STAIR_FORWARD) {
					if (j == which_interval) {
						if (j_pushed == false) {
							j_pushed = true;
						}
						else{
							continue;// if which_interval is already set, and check forward, then clear backward part
						}
						
					}
				}
				para_ids[j].push_back(i);
				break;
			}
		}
	}
	std::vector<int> need_fix_intervals;
	need_fix_intervals.push_back(which_interval);


	// down here, we need to insert one value to U
	std::cout << "pro_row_id " << pro_row_id << std::endl;
	std::cout << "before insert" << std::endl;
	Uout = knot_vector_insert_values(Uin, paras, need_fix_intervals, para_ids);
	std::cout << "finish insert" << std::endl;
	return;

}

// this function split the stair with the largest number of rows
// TODO i think this method is not convincing. it may not deal with the real problematic one, but may keep dealing 
// with the irrelevant part. imagine two stairs far a way, both very large, but only one is problematic. this method may
// keep split the large stair which is fine, but delay the processing of the problematic one.

// TODO considering this, do not use this function, although it can also provide a correct result, but not good enough
void insert_a_knot_to_a_stair_largest_rows(const int degree, const int pro_row_id, const int start_row, const int end_row,
	const std::vector<double>& Uin,
	const std::vector<double>& paras, std::vector<double>& Uout, const int construct_method = STAIR_WHOLE) {

	// Uin has n+degree+2 values, there are n-degree+2 different values, n-degree+1 intervals
	int n = Uin.size() - 2 - degree;

	Uout.clear();

	/*std::vector<int> multiplicity(Uin.size() - 1);*/
	std::vector<std::vector<int>> para_ids(Uin.size() - 1);

	int which_interval = -1; // gives which interval does pro_row_id corresponding to
	bool j_pushed = false;
	for (int i = 0; i < paras.size(); i++) {
		if (i<start_row || i>end_row) continue;
		for (int j = 0; j < Uin.size() - 1; j++) {
			if (paras[i] >= Uin[j] && paras[i] < Uin[j + 1]) {

				if (i == pro_row_id) {
					which_interval = j;
					
				}
				
				para_ids[j].push_back(i);
				break;
			}
		}
	}
	int highest_rows_interval = -1;
	int nbr_rows = -1;
	for (int i = 0; i < para_ids.size(); i++) {
		if (para_ids[i].size() > nbr_rows) {
			nbr_rows = para_ids[i].size();
			highest_rows_interval = i;
		}
	}
	assert(nbr_rows > 1);
	std::vector<int> need_fix_intervals;
	if (construct_method == STAIR_HIGHEST) {
		need_fix_intervals.push_back(highest_rows_interval);
	}
	else {
		need_fix_intervals.push_back(which_interval);
	}
	

	
	// down here, we need to insert one value to U

	Uout = knot_vector_insert_values(Uin, paras, need_fix_intervals, para_ids);

	return;

}

// this function uses boolean predicates to check if the equation is row full rank
// if row full rank, Ax = b have solution
bool curve_equation_row_full_rank(const int degree, const std::vector<double>& U, const std::vector<double>& paras,
	const int start_row, const int nbr_rows, int& rankdiff) {
	rankdiff = 0;
	std::vector<std::vector<int>> intervals;
	Eigen::MatrixXd A;
	A = build_matrix_A(degree, U, paras,start_row,nbr_rows);

	int r = 0;
	for (int i = 0; i < nbr_rows; i++) {
		bool found = false;
		for (int j = r; j < A.cols(); j++) {
			if (A(i, j) != 0) {
				r = j + 1;
				found = true;
				break;
			}
		}
		if (!found) {
			rankdiff = nbr_rows - i;
			return false;
		}
	}

	return true;
}

std::vector<double> fix_knot_vector_to_interpolate_curve_boolean(const int degree, const std::vector<double>& init_vec,
	const std::vector<double>& paras) {
	std::vector<double> expanded_U = init_vec;
	// take init_vec as input, detect if there is any stair whose row is too many (>n+1).
	// if there are such stairs, insert knots to init_vec, the result is expanded_U

	fix_stairs_row_too_many(degree, init_vec, paras, expanded_U);

	
	int n = expanded_U.size() - 2 - degree;// n + 1 = number of control points
	int m = paras.size() - 1;// m + 1 = the number of data points

	int start_row = 0;// initialized row nbr is 0
	int nbr_rows = m + 1;// initialized block has m+1 rows
	int rankdiff;
	bool have_solution = curve_equation_row_full_rank(degree, expanded_U, paras, start_row, nbr_rows, rankdiff);


	while (!have_solution) {

		
		
		int rank_diff1, rank_diff2;
		assert(nbr_rows > 1);
	
		bool have_solution1 = curve_equation_row_full_rank(degree, expanded_U, paras, start_row, nbr_rows - 1, rank_diff1);

	
		bool have_solution2 = curve_equation_row_full_rank(degree, expanded_U, paras, start_row + 1, nbr_rows - 1, rank_diff2);
		
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


			std::cout << "branch 1" << std::endl;
			insert_a_knot_to_a_stair(degree, start_row + nbr_rows-1, expanded_U, paras, tempU,STAIR_FORWARD);
			
		}
		if ((!have_solution1) && have_solution2) {
			// deal with knots, start_row=0, nbr_rows=m+1; calculate have_solution
			std::cout << "branch 2" << std::endl;
			insert_a_knot_to_a_stair(degree, start_row, expanded_U, paras, tempU,STAIR_BACKWARD);
			std::cout << "original" << std::endl;
			print_vector(expanded_U);
			std::cout << "fixed" << std::endl;
			print_vector(tempU);
			
		}
		if (have_solution1 && have_solution2) {
			// deal with knots, start_row=0, nbr_rows=m+1; calculate have_solution
			// if both of them have solution, pick a random one to solve
			std::cout << "branch 1" << std::endl;
			insert_a_knot_to_a_stair(degree, start_row + nbr_rows-1, expanded_U, paras, tempU, STAIR_FORWARD);
			
		}
		expanded_U = tempU;
		start_row = 0; nbr_rows = m + 1;


		have_solution = curve_equation_row_full_rank(degree,expanded_U,paras, start_row, nbr_rows,rankdiff);
	}
	return expanded_U;

}
std::vector<double> fix_knot_vector_to_interpolate_curve(const int degree, const std::vector<double>& init_vec,
	const std::vector<double>& paras, const std::vector<Vector3d>& points, const int dimension) {
	std::vector<double> expanded_U = init_vec;
	// take init_vec as input, detect if there is any stair whose row is too many (>n+1).
	// if there are such stairs, insert knots to init_vec, the result is expanded_U

	fix_stairs_row_too_many(degree, init_vec, paras, expanded_U);

	assert(points.size() == paras.size());
	int n = expanded_U.size() - 2 - degree;// n + 1 = number of control points
	int m = paras.size() - 1;// m + 1 = the number of data points
	int prob_id;
	Eigen::MatrixXd A;
	Eigen::VectorXd b;


	A = build_matrix_A(degree, expanded_U, paras);

	b = build_Vector_b(points, dimension);

	bool have_solution = equation_has_solution(A, b);
	//std::cout << "have solution? " << have_solution << std::endl;
	//std::cout << "can interpolate? " << curve_can_be_interpolated(expanded_U, degree, paras) << std::endl;
	if (curve_can_be_interpolated(expanded_U, degree, paras, prob_id)) {
		assert(have_solution);
	}
	if (!have_solution) {
		assert(!curve_can_be_interpolated(expanded_U, degree, paras, prob_id));
	}
	
	
	int start_row = 0;// initialized row nbr is 0
	int nbr_rows = m + 1;// initialized block has m+1 rows

	while (!have_solution) {


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


			std::cout << "branch 1" << std::endl;
			
			insert_a_knot_to_a_stair(degree, start_row + nbr_rows - 1, expanded_U, paras, tempU, STAIR_FORWARD);

		}
		if ((!have_solution1) && have_solution2) {
			// deal with knots, start_row=0, nbr_rows=m+1; calculate have_solution
			std::cout << "branch 2" << std::endl;
			insert_a_knot_to_a_stair(degree, start_row, expanded_U, paras, tempU, STAIR_BACKWARD);
			std::cout << "original" << std::endl;
			print_vector(expanded_U);
			std::cout << "fixed" << std::endl;
			print_vector(tempU);

		}
		if (have_solution1 && have_solution2) {
			// deal with knots, start_row=0, nbr_rows=m+1; calculate have_solution
			// if both of them have solution, pick a random one to solve
			std::cout << "branch 3" << std::endl;
			std::cout << "A\n" << A << std::endl; 
			insert_a_knot_to_a_stair(degree, start_row + nbr_rows - 1, expanded_U, paras, tempU, STAIR_FORWARD);

		}
		expanded_U = tempU;
		start_row = 0; nbr_rows = m + 1;
		A = build_matrix_A(degree, expanded_U, paras);
		b = build_Vector_b(points, dimension);

		have_solution = equation_has_solution(A, b);
		
		if (curve_can_be_interpolated(expanded_U, degree, paras, prob_id)) {
			assert(have_solution);
		}
		if (!have_solution) {
			assert(!curve_can_be_interpolated(expanded_U, degree, paras, prob_id));
		}

	}
	std::cout << "A\n" << A << std::endl;
	return expanded_U;

}
std::vector<double> fix_knot_vector_to_interpolate_curve(const int degree, const std::vector<double>& init_vec,
	const std::vector<double>& paras, const std::vector<Vector3d>& points) {
	std::vector<double> result = init_vec;
	for (int dim = 0; dim <3; dim++) {
		result = fix_knot_vector_to_interpolate_curve(degree, result, paras, points, dim);
	}
	return result;
}

bool curve_can_be_interpolated(const std::vector<double>& U,const int degree, const Eigen::VectorXd & paras, 
	int &prob_id) {
	prob_id = -1;
	std::vector<std::vector<int>> para_ids(U.size() - 1);// there are U.size()-1 intervals
	int nu = U.size() - 2 - degree;// n + 1 = number of control points
	for (int i = 0; i < paras.size(); i++) {
		for (int j = 0; j < U.size() - 1; j++) {
			if (paras[i] >= U[j] && paras[i] < U[j + 1]) {
				para_ids[j].push_back(i);
				break;
			}
		}
	}
	int k = 0;
	// now para_ids contains the information which parameter falls into which interval
	for (int i = 0; i < para_ids.size(); i++) {
		if (para_ids.size() == 0) {
			continue;
		}
		for (int j = 0; j < para_ids[i].size(); j++) {
			// now this parameter is in [U[i], U[i+1]]
			bool itp_this = false; // if this point can be interpolated
			while (!itp_this) {
				if (k > nu) {
					prob_id = para_ids[i][j];
					return false;
				}
				if (k >= i - degree && k <= i) {
					itp_this = true;
				}
				k++;
			}
		}
	}
	return true;

}
bool curve_can_be_interpolated(const std::vector<double>& U, const int degree, const std::vector<double> & paras,
	int &prob_id) {
	Eigen::VectorXd prs;
	prs.resize(paras.size());
	for (int i = 0; i < paras.size(); i++) {
		prs[i] = paras[i];
	}
	return curve_can_be_interpolated(U, degree, prs,prob_id);
}

// this is from [W.K. Wang, 2008, CAD]
bool curve_can_be_interpolated_wkw(const std::vector<double>& U, const int degree, const Eigen::VectorXd & paras,
	int &prob_id) {
	//TODO implement here. maybe no need to keep this function
	return false;
}

std::vector<double> merge_two_knot_vectors(const std::vector<double> &U1, const std::vector<double> &U2, const int degree) {
	std::vector<double>v;
	std::priority_queue<double, std::vector<double>, std::greater<double>> queue;
	
	for (int i = degree + 1; i < U1.size() - 1 - degree; i++) {
		queue.push(U1[i]);
	}
	for (int i = degree + 1; i < U2.size() - 1 - degree; i++) {
		queue.push(U2[i]);
	}
	for (int i = 0; i < degree + 1; i++) {
		v.push_back(U1.front());
	}
	while (!queue.empty()) {
		if (v.back() != queue.top()) {
			v.push_back(queue.top());
		}
		
		queue.pop();
	}
	for (int i = 0; i < degree + 1; i++) {
		v.push_back(U1.back());
	}
	//print_vector(v);
	return v;
}

// fix_nbr is a parameter that decide how many knot do we want to insert in this step. by default it is -1 which means
// the knot will be fully fixed. when fix_nbr>0, if insert more than fix_nbr knots, the function directly returns. 
std::vector<double> fix_knot_vector_to_interpolate_curve_WKW(const int degree, const std::vector<double>& init_vec,
	const std::vector<double>& paras, const double per,  bool &fully_fixed, const int fix_nbr) {
	std::vector<double> result;
	
	int s = paras.size() - 1;// s+1 is the para number
	int m = init_vec.size() - degree - 2;// m+degree+2 is the size 
	int hUsize = //s + degree + 2; // 
		std::max(s + degree + 2, 2 * degree + 2); // the minimal size should be 2*degree+2
	
	if (hUsize == 2 * degree + 2) {// if there are too few parameters to fit, just return the initial one;
		fully_fixed = true;
		return init_vec;
	}
	std::vector<double> hU(hUsize);
	for (int i = 0; i < degree + 1; i++) { // from 0 to degree, the total number is degree+1
		hU[i] = 0;
		result.push_back(0);
	}
	for (int i = hUsize - degree - 1; i < hUsize; i++) {// from size-degree-1 to size-1, degree + 1 elements 
		hU[i] = init_vec.back();
	}
	
	for (int i = degree + 1; i < s + 1; i++) {
		double add = 0;
		for (int j = i - degree; j < i; j++) {
			add += paras[j];
		} 
		hU[i] = add / degree;
	}

	bool breakloop = false;
	int interval = 0;
	int nbr_of_inserted = 0;
	for (int i = degree + 1; i < s + 1; i++) {// from degree+1 to s
		//std::cout << "\ndealing with the ith " << i << std::endl;
		double alpha = (paras[i - 1] - (paras[i - degree] + paras[i - degree - 1]) / 2) / 
			(paras[i - 1] - paras[i - degree - 1]);
		double belta = ((paras[i] + paras[i - 1]) / 2 - paras[i - degree]) /
			(paras[i] - paras[i - degree]);
		double a = (1 - per * alpha)*hU[i] + per * alpha*hU[i - 1];
		double b = (1 - per * belta)*hU[i] + per * belta*hU[i + 1];
		while (init_vec[i + interval] < a) {
			interval += 1;
			if (init_vec[i + interval] >= a) {
				//std::cout << "front interval, " << interval << std::endl;
			}
			if (i + interval > m) { // when i + interval = m+1, it will be 1
				breakloop = true;
				//std::cout << "break here" << std::endl;
				break;
			}
		}
		if (breakloop) {
			break;
		}
		//std::cout << "back interval, " << interval << std::endl;
		if (init_vec[i + interval] <= b) {
			result.push_back(init_vec[i + interval]);
			//std::cout << "directly push value satisfies" << std::endl;
		}
		else {
			result.push_back(hU[i]);
			interval -= 1;
			nbr_of_inserted += 1;
			if (fix_nbr >= 0) {
				if (nbr_of_inserted > fix_nbr) {
					break;
				}
			}
			//std::cout << "push value not satisfies" << std::endl;
		}
	}
	int nowsize = result.size();
	for (int i = nowsize; i < s + 1; i++) {
		nbr_of_inserted += 1;
		result.push_back(hU[i]);
		if (fix_nbr >= 0) {
			if (nbr_of_inserted > fix_nbr) {
				break;
			}
		}
	}
	for (int i = s + 1; i < s + degree + 2; i++) {
		result.push_back(hU[i]);
	}
	//print_vector(hU);
	//print_vector(result);
	
	if (result.size() == hU.size()) {
		fully_fixed = true;
	}
	else {
		
		fully_fixed = false;
	}
	if (fix_nbr < 0) {
		//fully_fixed = true;
		assert(result.size() == hU.size());// by default, the two sizes should be the same. but if not fully fixed, this is not true
	}
	
	//std::cout << "check A validation\n"<< build_matrix_A(degree,result,paras) << std::endl;
	result = merge_two_knot_vectors(result, init_vec, degree);
	
	return result;

}

// TODO add a parameter weight to control the percentage
// per in [0,1]. when per =0, we have larger tolerance for feasible points, which means it converges to tranditional methods
std::vector<int> feasible_control_point_of_given_parameter(const double para, const std::vector<double>&U, 
	const int degree, const double per) {
	int which = -1;
	int nu = U.size() - 2 - degree;// n + 1 = number of control points
	std::vector<int> result;
	if (para != U.back()) {// when para !=1, we check which interval it is in
		for (int i = 0; i < U.size() - 1; i++) {
			if (para >= U[i] && para < U[i + 1]) {
				which = i;
				break;
			}
		}
	}
	else {// if para == 1, the nth control point is feasible 
		result.resize(1);
		result[0] = nu;
		return result;
	}
	// if which >= 0, the para is in [u_which, u_{which+1}).
	// alpha, belta are slightly less than 1
	double alpha_pre = (U[which + degree] - (U[which] + U[which + 1]) / 2) / (U[which + degree] - U[which]);
	double belta_pre = ((U[which + 1] + U[which]) / 2 - U[which - degree + 1]) / (U[which + 1] - U[which - degree + 1]);
	double alpha = (1 - alpha_pre)*(1 - per) + alpha_pre;
	double belta = (1 - belta_pre)*(1 - per) + belta_pre;
	double a = U[which + 1] + alpha * (U[which] - U[which + 1]);
	double b = U[which] + belta * (U[which + 1] - U[which]);
	// U[i]<a<b<U[i+1]
	result.resize(0);
	if (para == U.front()) {
		result.push_back(0);
		return result;
	}
	if (para == U.back()) {
		result.push_back(nu);
		return result;
	}

	std::string compare;
	if (para <= b) { // if para not too close to U[i+1], then N(i-p) is feasible
		result.push_back(which - degree);
	}
	else {
		/*if (para > U[which + 1]) {
			compare = ">";
		}
		if (para < U[which + 1]) {
			compare = "<";
		}
		if (para == U[which + 1]) {
			compare = "=";
		}*/
		//std::cout << para << compare<<" is too close to " << U[which + 1]<<"i-p is "<<which-degree<<" N(i-p) is "<<Nip(which-degree,degree,para,U) << std::endl;
	}
	for (int k = which - degree + 1; k < which; k++) {// from i-p to i-1
		result.push_back(k);
	}
	if (para >= a) { // if para not too close to U[i], then N(i) is feasible
		result.push_back(which);
	}
	else {
		/*if (para > U[which]) {
			compare = ">";
		}
		if (para < U[which]) {
			compare = "<";
		}
		if (para == U[which]) {
			compare = "=";
		}*/
		//std::cout << para <<compare<< " is too close to " << U[which] << "i is " << which  << " N(i) is " << Nip(which, degree, para, U) << std::endl;
	}
	assert(result.size() > 0);
	return result;
}

