#include<energy.h>
#include<surface.h> 
#include<igl/Timer.h>
#include <Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/SparseLU>
igl::Timer timer;
double time0 = 0, time1 = 0, time2 = 0, time3 = 0;
std::vector<double> polynomial_simplify(const std::vector<double>& poly) {
	std::vector<double> result = poly;
	int size = poly.size();
	for (int i = 0; i < size - 1; i++) {
		if (result[size - i - 1] == 0) {
			result.pop_back();
		}
		else {
			break;
		}
	}
	return result;
}

std::vector<double> polynomial_add(const std::vector<double>& poly1, const std::vector<double>& poly2) {
	int size = std::max(poly1.size(), poly2.size());
	std::vector<double> result(size);
	for (int i = 0; i < size; i++) {
		bool flag1 = i < poly1.size();
		bool flag2 = i < poly2.size();
		if (flag1 && flag2) {
			result[i] = poly1[i] + poly2[i];
		}
		else if (flag1) {
			result[i] = poly1[i];
		}
		else {
			result[i] = poly2[i];
		}
	}
	return polynomial_simplify(result);
}
std::vector<double> polynomial_times(const std::vector<double>& poly1, const std::vector<double>& poly2) {
	int size = poly1.size() + poly2.size() - 1;
	std::vector<double> result(size);
	for (int i = 0; i < size; i++) {// initialize the result
		result[i] = 0;
	}

	for (int i = 0; i < poly1.size(); i++) {
		for (int j = 0; j < poly2.size(); j++) {
			result[i + j] += poly1[i] * poly2[j];
		}
	}
	return polynomial_simplify(result);
}
std::vector<double> polynomial_times(const std::vector<double>& poly1, const double& nbr) {
	std::vector<double> result;
	if (nbr == 0) {
		result.resize(1);
		result[0] = 0;
		return result;
	}
	result = poly1;
	for (int i = 0; i < result.size(); i++) {
		result[i] *= nbr;
	}

	return polynomial_simplify(result);
}
std::vector<double> Ni0_func(const int i, const double u, const std::vector<double> &U) {
	std::vector<double> result(1);
	if (u >= U[i] && u < U[i + 1]) {
		result[0] = 1;
		return result;
	}
	result[0] = 0;
	return result;
}
std::vector<double> handle_division_func(const std::vector<double>& a, const double b) {
	std::vector<double> result;
	if (b == 0) {
		result.resize(1);
		result[0] = 0;
		// if the denominator is 0, then this term is 0
		return result;
	}
	else return polynomial_times(a, 1/b);
}

std::vector<double> Nip_func(const int i, const int p, const double u, const std::vector<double> &U) {
	std::vector<double> result;
	if (p == 0) {
		return Ni0_func(i, u, U);
	}
	if (u == U.back()) {
		result.resize(1);
		if (i == U.size() - 2 - p) {
			result[0] = 1;
			return result;
		}
		else {
			result[0] = 0;
			return result;
		}
	}
	std::vector<double> v;
	v = { {-U[i],1} };// u - U[i]
	std::vector<double> result1 = polynomial_times(handle_division_func(v, U[i + p] - U[i]), Nip_func(i, p - 1, u, U));
	/*std::cout << "**this degree " << p << std::endl;
	std::cout << "v "  << std::endl;
	print_vector(v);
	std::cout << "lower degree " << p-1 << std::endl;
	print_vector(Nip(i, p - 1, u, U));*/

	v = { {U[i + p + 1],-1} };// U[i+p+1] - u 
	std::vector<double> result2 = polynomial_times(handle_division_func(v, U[i + p + 1] - U[i + 1]), Nip_func(i + 1, p - 1, u, U));
	return polynomial_add(result1, result2);
}
double polynomial_value(const std::vector<double>& poly, const double para) {
	double result=0;
	for (int i = 0; i < poly.size(); i++) {
		result += poly[i] * std::pow(para, i);
	}
	return result;
}
std::vector<double> polynomial_integration(const std::vector<double>& poly) {
	std::vector<double> result(poly.size() + 1);
	result[0] = 0;
	for (int i = 1; i < result.size(); i++) {
		result[i] = poly[i - 1] / i;
	}
	return polynomial_simplify(result);
}
double polynomial_integration(const std::vector<double>& poly, const double lower, const double upper) {
	double up = polynomial_value(polynomial_integration(poly), upper);
	double lw = polynomial_value(polynomial_integration(poly), lower);
	return up - lw;
}

// order 1 differential
const std::vector<double> polynomial_differential(const std::vector<double>& func) {
	std::vector<double> result;
	if (func.size() == 1) {
		result.resize(1);
		result[0] = 0;
		return result;
	}
	result.resize(func.size() - 1);
	for (int i = 0; i < result.size(); i++) {
		result[i] = func[i + 1] * (i + 1);
	}
	return result;
	
}
const std::vector<double> polynomial_differential(const std::vector<double>& func, const int order) {
	std::vector<double> result = func;
	if (order == 0) return result;
	if (func.size() == 1 && order > 0) {
		result.resize(1);
		result[0] = 0;
		return result;
	}
	std::vector<double> tmp;
	for (int i = order; i > 0; i--) {
		tmp = polynomial_differential(result);
		result = tmp;
	}
	return result;
}

void PolynomialBasis::init(Bsurface& surface) {
	Uknot = surface.U;
	Vknot = surface.V;
	degree1 = surface.degree1;
	degree2 = surface.degree2;
	nu = surface.nu();
	nv = surface.nv();
	Ubasis = calculate(0);
	Vbasis = calculate(1);
	inited = true;
	return;
}
void PolynomialBasis::clear() {
	Uknot.clear();
	Vknot.clear();
	Ubasis.clear();
	Vbasis.clear();
	inited = false;
}
PolynomialBasis::PolynomialBasis(Bsurface& surface) {
	init(surface);
}
PolynomialBasis::PolynomialBasis() {
}
std::vector<double> PolynomialBasis::poly(const int id, const double value, const bool UVknot) {
	if (!inited) {
		std::cout << "WRONG USAGE OF CLASS PolynomialBasis, YOU SHOULD INITIALIZE IT BY CALLING init()" << std::endl;
	}
	
	std::vector<double> kv;
	int degree;
	
	if (UVknot) {
		kv = Vknot;
		degree = degree2;
	}
	else {
		kv = Uknot;
		degree = degree1;

	}

	int which = -1;
	for (int i = 0; i < kv.size() - 1; i++) {
		if (value >= kv[i] && value < kv[i + 1]) {
			which = i;
			break;
		}
	}
	if (which == -1) {
		std::cout << "ERROR: DON'T USE POLYNOMIAL WHEN VALUE = " << value << std::endl;
		exit(0);
	}
	// the value is in [U[i], U[i+1]), the Nip are from N(i-p) to N(i)
	std::vector<double> result;
	int index = id - (which - degree);
	if (UVknot) {// check v
		result = Vbasis[which][index];
	}
	else {
		result = Ubasis[which][index];
	}
	return result;
}

std::vector<std::vector<std::vector<double>>> PolynomialBasis::calculate(const bool uorv) {
	std::vector<double> kv;
	int degree;
	int n;
	if (uorv) {
		kv = Vknot;
		degree = degree2;
		n = nv;
	}
	else {
		kv = Uknot;
		degree = degree1;
		n = nu;
	}
	std::vector<std::vector<std::vector<double>>> pl;
	pl.resize(n + 1);// n+1 intervals;
	for (int i = 0; i < n+1; i++) {// in interval [U[i], U[i+1])
		pl[i].resize(degree + 1);
		for (int j = 0; j < degree + 1; j++) {
			pl[i][j] = Nip_func(i - degree + j, degree, kv[i], kv);
		}
	}
	return pl;

}

// construct an integration of multiplication of two B-spline basis (intergration of partial(Ni1)*partial(Ni2))
// the integration domain is [u1, u2]
double construct_an_integration(const int degree, const std::vector<double>& U,
	const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2, 
	PolynomialBasis &basis, const bool uv) {
	timer.start();
	std::vector<double> func1 = basis.poly(i1, u1, uv);
	//Nip_func(i1, degree, u1, U);
	std::vector<double> func2 = basis.poly(i2, u1, uv);
	//Nip_func(i2, degree, u1, U);
	/*if (basis.poly(i1, u1, uv) != func1) {
		std::cout << "NOT EQUAL" << std::endl;
		exit(0);
	}
	if (basis.poly(i2, u1, uv) != func2) {
		std::cout << "NOT EQUAL" << std::endl;
		exit(0);
	}
	assert(basis.poly(i1, u1, uv) == func1);
	assert(basis.poly(i2, u1, uv) == func2);*/
	timer.stop();
	time0 += timer.getElapsedTimeInMilliSec();
	//std::cout << "degree, "<<degree << std::endl;
	//std::cout << "func1 and func2" << std::endl;
	//print_vector(func1);
	//print_vector(func2);
	timer.start();
	func1 = polynomial_differential(func1, partial1);
	func2 = polynomial_differential(func2, partial2);
	timer.stop();
	time1 += timer.getElapsedTimeInMilliSec();
	//std::cout << "differencial" << std::endl;
	//print_vector(func1);
	//print_vector(func2);
	timer.start();
	std::vector<double> func = polynomial_times(func1, func2);
	//std::cout << "times" << std::endl;
	//print_vector(func);
	double upper = u2;
	if (u2 == U.back()) {
		upper = U.back() - SCALAR_ZERO;
	}
	double result = polynomial_integration(func, u1, upper);
	timer.stop();
	time2 += timer.getElapsedTimeInMilliSec();

	return result;
}
// construct an integration of multiplication of two B-spline basis (intergration of partial(Ni1)*partial(Ni2))
// the integration domain is [u1, u2]
double construct_an_integration(const int degree, const std::vector<double>& U,
	const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2) {
	timer.start();
	std::vector<double> func1 = Nip_func(i1, degree, u1, U);
	std::vector<double> func2 = Nip_func(i2, degree, u1, U);
	timer.stop();
	time0 += timer.getElapsedTimeInMilliSec();
	//std::cout << "degree, "<<degree << std::endl;
	//std::cout << "func1 and func2" << std::endl;
	//print_vector(func1);
	//print_vector(func2);
	timer.start();
	func1 = polynomial_differential(func1, partial1);
	func2 = polynomial_differential(func2, partial2);
	timer.stop();
	time1 += timer.getElapsedTimeInMilliSec();
	//std::cout << "differencial" << std::endl;
	//print_vector(func1);
	//print_vector(func2);
	timer.start();
	std::vector<double> func = polynomial_times(func1, func2);
	//std::cout << "times" << std::endl;
	//print_vector(func);
	double upper = u2;
	if (u2 == U.back()) {
		upper = U.back() - SCALAR_ZERO;
	}
	double result = polynomial_integration(func, u1, upper);
	timer.stop();
	time2 += timer.getElapsedTimeInMilliSec();

	return result;
}
// do partial difference to Pi, the cofficient of jth element Pj.
double surface_energy_least_square( Bsurface& surface, const int i, const int j, PolynomialBasis& basis) {
	// figure out which Pij corresponding to the ith control point
	int partial_i = i / (surface.nv()+1);
	int partial_j = i - partial_i * (surface.nv()+1);

	// figure out which Pij corresponding to the jth control point
	int coff_i = j / (surface.nv() + 1);
	int coff_j = j - coff_i * (surface.nv() + 1);

	// if do partial Pij, the related other Pi1j1 will be: i1 = i-p~i+p, j1= j-q~j+q
	int degree1 = surface.degree1, degree2 = surface.degree2;
	if (coff_i<partial_i - degree1 || coff_i>partial_i + degree1) return 0;
	if (coff_j<partial_j - degree2 || coff_j>partial_j + degree2) return 0;


	// do partial Pij
	// for each block, (U_k1, U_(k1+1)), (V_k2, V_(k2+1)). the related control points are k1-p,...,k1 and k2-q,..., k2
	double result = 0;
	for (int k1 = partial_i; k1 < partial_i + degree1 + 1; k1++) {
		for (int k2 = partial_j; k2 < partial_j + degree2 + 1; k2++) {
			if (coff_i<k1 - degree1 || coff_i>k1 || coff_j<k2 - degree2 || coff_j>k2) {
				continue;
			}
			if (surface.U[k1] == surface.U[k1 + 1] || surface.V[k2] == surface.V[k2 + 1]) {
				continue;// if the area is 0, then no need to compute
			}
			// Suu part
			double value1 = construct_an_integration(degree1, surface.U, 2, 2, 
				partial_i, coff_i, surface.U[k1], surface.U[k1 + 1],basis,false);
			double value2 = construct_an_integration(degree2, surface.V, 0, 0, 
				partial_j, coff_j, surface.V[k2], surface.V[k2 + 1], basis, true);
			// Suv part
			double value3 = construct_an_integration(degree1, surface.U, 1, 1,
				partial_i, coff_i, surface.U[k1], surface.U[k1 + 1], basis, false);
			double value4 = construct_an_integration(degree2, surface.V, 1, 1,
				partial_j, coff_j, surface.V[k2], surface.V[k2 + 1], basis, true);
			// Svv part
			double value5 = construct_an_integration(degree1, surface.U, 0, 0,
				partial_i, coff_i, surface.U[k1], surface.U[k1 + 1], basis, false);
			double value6 = construct_an_integration(degree2, surface.V, 2, 2,
				partial_j, coff_j, surface.V[k2], surface.V[k2 + 1], basis, true);
			result += 2 * value1*value2 + 4 * value3*value4 + 2 * value5*value6;
		}
	}
	return result;
	
}


// which_part = 0: Suu; which_part = 1, Suv; which_part = 2, Svv.
Eigen::MatrixXd energy_part_of_surface_least_square(Bsurface& surface, PolynomialBasis& basis) {
	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	Eigen::MatrixXd result(psize, psize);
	for (int i = 0; i < psize; i++) {
		//std::cout << "the ith row of matrix" << std::endl;
		for (int j = 0; j < psize; j++) {
			result(i, j) = surface_energy_least_square(surface, i, j,basis);
		}
	}
	std::cout << "energy matrix finish calculation" << std::endl;
	return result;
}

Eigen::MatrixXd eqality_part_of_surface_least_square(Bsurface& surface, const Eigen::MatrixXd& paras) {
	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	Eigen::MatrixXd result(paras.rows(), psize);
	int degree1 = surface.degree1;
	int degree2 = surface.degree2;
	std::vector<double> U = surface.U;
	std::vector<double> V = surface.V;
	for (int i = 0; i < result.rows(); i++) {
		for (int j = 0; j < result.cols(); j++) {
			// figure out the jth control point corresponding to which Pij
			int coff_i = j / (surface.nv() + 1);
			int coff_j = j - coff_i * (surface.nv() + 1);
			double u = paras(i, 0);
			double v = paras(i, 1);
			// the corresponding cofficient should be N_coffi(u) and N_coffj(v)
			double N1 = Nip(coff_i, degree1, u, U);
			double N2 = Nip(coff_j, degree2, v, V);
			result(i, j) = N1 * N2;
		}
	}
	return result;
}

Eigen::MatrixXd lambda_part_of_surface_least_square(Bsurface& surface, const Eigen::MatrixXd& paras) {
	Eigen::MatrixXd A = eqality_part_of_surface_least_square(surface, paras);
	Eigen::MatrixXd result = -A.transpose();
	return result;
}

Eigen::MatrixXd surface_least_square_lambda_multiplier_left_part(Bsurface& surface, const Eigen::MatrixXd& paras,PolynomialBasis& basis) {
	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	int target_size = paras.rows();// nbr of target data points
	int size = psize + target_size;
	Eigen::MatrixXd result(size, size);
	Eigen::MatrixXd rd = Eigen::MatrixXd::Zero(target_size, target_size);// right down corner part
	Eigen::MatrixXd lu = energy_part_of_surface_least_square(surface,basis);
	Eigen::MatrixXd ru = lambda_part_of_surface_least_square(surface, paras);
	Eigen::MatrixXd ld = eqality_part_of_surface_least_square(surface, paras);
	std::cout << "sizes" << std::endl;
	std::cout << "lu, " << lu.rows() << " " << lu.cols() << std::endl;
	std::cout << "ru, " << ru.rows() << " " << ru.cols() << std::endl;
	std::cout << "ld, " << ld.rows() << " " << ld.cols() << std::endl;
	std::cout << "rd, " << rd.rows() << " " << rd.cols() << std::endl;
	result << lu, ru,
		ld, rd;
	return result;
}

Eigen::MatrixXd surface_least_square_lambda_multiplier_right_part(Bsurface& surface, const Eigen::MatrixXd& paras,
	const Eigen::MatrixXd & points, const int dimension) {
	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	int target_size = paras.rows();// nbr of target data points
	int size = psize + target_size;
	Eigen::MatrixXd result(size, 1);
	for (int i = 0; i < psize; i++) {
		result(i, 0) = 0;
	}
	int counter = 0;
	for (int i = psize; i < size; i++) {
		result(i, 0) = points(counter, dimension);
		counter++;
	}
	assert(counter == target_size);
	return result;
}

void push_control_point_list_into_surface(Bsurface& surface, const std::vector<Vector3d>& cps) {
	int id = 0;
	std::vector<std::vector<Vector3d>> control;
	control.resize(surface.nu() + 1);
	int vs = surface.nv() + 1;
	for (int i = 0; i < control.size(); i++) {
		control[i].resize(vs);
		for (int j = 0; j < vs; j++) {
			control[i][j] = cps[id];
			id++;
		}
	}
	surface.control_points = control;
	return;
}
void solve_control_points_for_fairing_surface(Bsurface& surface, const Eigen::MatrixXd& paras,
	const Eigen::MatrixXd & points, PolynomialBasis& basis) {
	bool sparse_solve = true;
	using namespace Eigen;
	typedef SparseMatrix<double> SparseMatrixXd;
	assert(paras.rows() == points.rows());
	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	std::vector<Vector3d> cps(psize);// control points
	Eigen::MatrixXd A;
	Eigen::FullPivLU<DenseBase<MatrixXd>::PlainMatrix> decomp;
	A = surface_least_square_lambda_multiplier_left_part(surface, paras,basis);
	SparseMatrixXd matB;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	if (sparse_solve) {
		matB = A.sparseView();

		solver.compute(matB);
		if (solver.info() != Success) {
			// decomposition failed
			std::cout << "solving failed" << std::endl;
			return;
		}
	}
	else {
		decomp = A.fullPivLu();
	}
	
	
	for (int i = 0; i < 3; i++) {
		Eigen::MatrixXd b = surface_least_square_lambda_multiplier_right_part(surface, paras, points, i);
		double err = 0.0;

		// solve the matrix contains the p and lambda
		std::cout << "before solving" << std::endl;
		Eigen::MatrixXd p_lambda;
		if (sparse_solve) {
			p_lambda = solver.solve(b);
			if (solver.info() != Success) {
				std::cout << "solving failed" << std::endl;
				return;
			}
		}
		else {
			p_lambda = decomp.solve(b);
		}
		
		double relative_error = (A*p_lambda - b).norm() / b.norm(); // norm() is L2 norm
		std::cout << "after solving, error is "<<relative_error << std::endl;
		push_p_lambda_vector_to_control_points(p_lambda, i, cps);
	}
	push_control_point_list_into_surface(surface, cps);
	return;
}
void output_timing() {
	std::cout << "get basis time " << time0 << std::endl;
	std::cout << "get differential time " << time1 << std::endl;
	std::cout << "get integration time " << time2 << std::endl;
}