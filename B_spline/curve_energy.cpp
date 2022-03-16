#include<sparse_interp/energy.h>
#include<sparse_interp/curve.h>
// the ith row, jth column element of the curve energy matrix of function applying least-square to Intergrate(||C'||^2+||C''||^2). 
// which means do partial difference to Pi, and 
// the cofficient of Pj. a and b are the parameters needed for the 1 order and 2 order differential part separately
namespace SIBSplines{
double curve_energy_least_square(const Bcurve& curve, const int i, const int j, const double a, const double b) {

	int degree = curve.degree;
	double result = 0;

	// iterate over all the intervals associate with Pi
	for (int k1 = i; k1 < i + degree + 1; k1++) {
		// now the interval is [U[k1], U[k1+1])
		if (j<k1 - degree || j>k1) {
			continue;
		}
		if (curve.U[k1] == curve.U[k1 + 1]) {
			continue;
		}
		double value1 = a == 0 ? 0 : construct_an_integration(degree, curve.U, 1, 1, i, j, curve.U[k1], curve.U[k1 + 1]);
		
	
		double value2 = b == 0 ? 0 : construct_an_integration(degree, curve.U, 2, 2, i, j, curve.U[k1], curve.U[k1 + 1]);
		//value2 = 0;
		//value1 = 0;
		
		result += a*value1 + b*value2;
		
	}
	result *= 2; // the cofficient before the integrations is 2
	return result;
}


Eigen::MatrixXd energy_part_of_curve_least_square( Bcurve& curve,const double a, const double b) {
	int psize = (curve.nu() + 1);// total number of control points.
	Eigen::MatrixXd result(psize, psize);
	for (int i = 0; i < psize; i++) {
		for (int j = 0; j < psize; j++) {
			result(i, j) = curve_energy_least_square(curve, i, j, a, b);
		}
	}
	return result;
}


Eigen::MatrixXd eqality_part_of_curve_least_square(Bcurve& curve, const std::vector<double> &paras) {
	int rows = paras.size();// how many points are there to be interpolated
	int cols = curve.nu() + 1;// how many control points 
	Eigen::MatrixXd result(rows, cols);
	for (int i = 0; i < rows; i++) {// for every parameter
		for (int j = 0; j < cols; j++) {// for every control point
			result(i, j) = Nip(j, curve.degree, paras[i], curve.U);
		}
	}
	return result;
}
Eigen::MatrixXd lambda_part_of_curve_least_square(Bcurve& curve, const std::vector<double>& paras) {
	int rows = curve.nu() + 1;
	int cols= paras.size();
	Eigen::MatrixXd result(rows, cols);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			result(i, j) = -Nip(i, curve.degree, paras[j], curve.U);
		}
	}
	return result;
}

// assemble the matrix of least square minimizing the curve energy and using lambda multiplier method
// to solve the equality constrains. a and b are the parameters for strengthing and bending parts of the energy.
Eigen::MatrixXd curve_least_square_lambda_multiplier_left_part(Bcurve& curve, const std::vector<double>& paras, 
	const double a,const double b) {
	int rows = curve.nu() + 1;
	int cols = paras.size();
	int size = rows + cols;
	Eigen::MatrixXd result(size, size);
	Eigen::MatrixXd rd = Eigen::MatrixXd::Zero(cols, cols);// right down corner part
	Eigen::MatrixXd lu = energy_part_of_curve_least_square(curve, a, b);
	Eigen::MatrixXd ru = lambda_part_of_curve_least_square(curve, paras);
	Eigen::MatrixXd ld = eqality_part_of_curve_least_square(curve, paras);
	assert(lu.rows() == rows && lu.cols() == rows && ld.rows() == cols && ld.cols() == rows);
	assert(ru.rows() == rows && ru.cols() == cols);
	result << lu, ru,
		ld, rd;
	return result;
}

// dimension implies if we are solving x, y or z coordinates
Eigen::MatrixXd curve_least_square_lambda_multiplier_right_part(Bcurve& curve, const std::vector<double>& paras,
	const Eigen::MatrixXd & points, const int dimension) {
	int rows = curve.nu() + 1;
	int cols = paras.size();
	int size = rows + cols;
	Eigen::MatrixXd result(size, 1);
	for (int i = 0; i < rows; i++) {
		result(i, 0) = 0;
	}
	int counter = 0;
	for (int i = rows; i < size; i++) {
		result(i, 0) = points(counter, dimension);
		counter++;
	}
	assert(counter == cols);
	return result;
}
void push_p_lambda_vector_to_control_points(const Eigen::MatrixXd &pl, 
	  const int dimension, std::vector<Vector3d>& control_points) {
	assert(pl.cols() == 1);
	for (int i = 0; i < control_points.size(); i++) {
		control_points[i][dimension] = pl(i, 0);
	}
}

// the output is the curve.control_points
// trying to find a curve minimizing the energy, while interpolating the points whose parameters are paras.
void Bcurve::solve_control_points_for_fairing_curve(Bcurve& curve, const std::vector<double>& paras,
	const Eigen::MatrixXd & points, const double a, const double b) {
	assert(paras.size() == points.rows());
	std::vector<Vector3d> cps(curve.nu() + 1);// control points
	Eigen::MatrixXd A = curve_least_square_lambda_multiplier_left_part(curve, paras, a, b);
	std::cout << "print A\n"<<A << std::endl;

	for (int i = 0; i < 3; i++) {
		Eigen::MatrixXd b = curve_least_square_lambda_multiplier_right_part(curve, paras, points, i);
		double err = 0.0;

		// solve the matrix contains the p and lambda
		std::cout << "before solving" << std::endl;
		Eigen::MatrixXd p_lambda = slove_linear_system(A, b, false, err);
		std::cout << "after solving" << std::endl;
		push_p_lambda_vector_to_control_points(p_lambda, i, cps);
	}
	curve.control_points = cps;
}
void Bcurve::solve_control_points_for_fairing_curve(Bcurve& curve, const std::vector<double>& paras,
	const std::vector<Vector3d> & pts, const double a, const double b) {
	Eigen::MatrixXd points = vector_to_matrix_3d(pts);
	solve_control_points_for_fairing_curve(curve, paras, points, a, b);
	return;
}
}