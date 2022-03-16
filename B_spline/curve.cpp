#include<sparse_interp/curve.h>
#include <Eigen/Dense>
#include<iostream>
#include<queue>
namespace SIBSplines{
int Bcurve::nu() {
	return U.size() - 2 - degree;
}
// return a point position of a given curve
Vector3d Bcurve::BsplinePoint(const int degree, const std::vector<double>& U, const double para,
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
Vector3d Bcurve::BsplinePoint(const int degree, const std::vector<double>& U, const double para,
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
// TODO maybe merge this with surface solving
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

// least square method, use NTN*D = R to solve the n-1 control points D
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

// least square method solve the control points
Eigen::MatrixXd Bcurve::solve_curve_control_points(const int degree, const std::vector<double>& U,
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



bool Bcurve::curve_can_be_interpolated(const std::vector<double>& U,const int degree, const Eigen::VectorXd & paras, 
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
bool Bcurve::curve_can_be_interpolated(const std::vector<double>& U, const int degree, const std::vector<double> & paras,
	int &prob_id) {
	Eigen::VectorXd prs;
	prs.resize(paras.size());
	for (int i = 0; i < paras.size(); i++) {
		prs[i] = paras[i];
	}
	return curve_can_be_interpolated(U, degree, prs,prob_id);
}

// the following functions are from [W.K. Wang, 2008, CAD]
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
std::vector<double> Bcurve::fix_knot_vector_to_interpolate_curve(const int degree, const std::vector<double>& init_vec,
	const std::vector<double>& paras, const double per,  bool &fully_fixed, const int fix_nbr){
		return fix_knot_vector_to_interpolate_curve_WKW(degree, init_vec,paras,per, fully_fixed, fix_nbr);
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
}
