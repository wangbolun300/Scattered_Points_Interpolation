#include"curve.h"
#include <Eigen/Dense>
#include<iostream>

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

	Eigen::MatrixXd A;
	Eigen::VectorXd b;


	A = build_matrix_A(degree, expanded_U, paras);

	b = build_Vector_b(points, dimension);

	bool have_solution = equation_has_solution(A, b);

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