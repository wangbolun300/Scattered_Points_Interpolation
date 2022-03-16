#include<sparse_interp/energy.h>
#include<sparse_interp/surface.h> 
#include<igl/Timer.h>
#include <Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/SparseLU>
#include<sparse_interp/Types.hpp>
#include<cmath>
namespace SIBSplines{
igl::Timer timer;
double time0 = 0, time1 = 0, time2 = 0, time3 = 0;
std::vector<double> ply_operations::polynomial_simplify(const std::vector<double>& poly) {
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

std::vector<double> ply_operations::polynomial_add(const std::vector<double>& poly1, const std::vector<double>& poly2) {
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
	return ply_operations::polynomial_simplify(result);
}
std::vector<double> ply_operations::polynomial_times(const std::vector<double>& poly1, const std::vector<double>& poly2) {
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
	return ply_operations::polynomial_simplify(result);
}
std::vector<double> ply_operations::polynomial_times(const std::vector<double>& poly1, const double& nbr) {
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

	return ply_operations::polynomial_simplify(result);
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
	else return ply_operations::polynomial_times(a, 1 / b);
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
	std::vector<double> result1 = ply_operations::polynomial_times(handle_division_func(v, U[i + p] - U[i]), Nip_func(i, p - 1, u, U));
	/*std::cout << "**this degree " << p << std::endl;
	std::cout << "v "  << std::endl;
	print_vector(v);
	std::cout << "lower degree " << p-1 << std::endl;
	print_vector(Nip(i, p - 1, u, U));*/

	v = { {U[i + p + 1],-1} };// U[i+p+1] - u 
	std::vector<double> result2 = ply_operations::polynomial_times(handle_division_func(v, U[i + p + 1] - U[i + 1]), Nip_func(i + 1, p - 1, u, U));
	return ply_operations::polynomial_add(result1, result2);
}
double ply_operations::polynomial_value(const std::vector<double>& poly, const double para) {
	double result = 0;
	for (int i = 0; i < poly.size(); i++) {
		result += poly[i] * std::pow(para, i);
	}
	return result;
}
std::vector<double> ply_operations::polynomial_integration(const std::vector<double>& poly) {
	std::vector<double> result(poly.size() + 1);
	result[0] = 0;
	for (int i = 1; i < result.size(); i++) {
		result[i] = poly[i - 1] / i;
	}
	return ply_operations::polynomial_simplify(result);
}
double ply_operations::polynomial_integration(const std::vector<double>& poly, const double lower, const double upper) {
	double up = ply_operations::polynomial_value(ply_operations::polynomial_integration(poly), upper);
	double lw = ply_operations::polynomial_value(ply_operations::polynomial_integration(poly), lower);
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
		assert(which < Vbasis.size());
		assert(index < Vbasis[which].size());
		result = Vbasis[which][index];
	}
	else {
		assert(which < Ubasis.size());
		assert(index < Ubasis[which].size());
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

	for (int i = degree; i < n + 1; i++) {// in interval [U[i], U[i+1])
		pl[i].resize(degree + 1);
		for (int j = 0; j < degree + 1; j++) {
			pl[i][j] = Nip_func(i - degree + j, degree, kv[i], kv);
		}
	}

	return pl;

}
std::vector<std::vector<std::vector<double>>> PartialBasis::do_partial(const
	std::vector<std::vector<std::vector<double>>>&basis) {
	std::vector<std::vector<std::vector<double>>> result(basis.size());
	for (int i = 0; i < basis.size(); i++) {
		result[i].resize(basis[i].size());
		for (int j = 0; j < basis[i].size(); j++) {
			result[i][j] = polynomial_differential(basis[i][j]);
		}
	}
	return result;
}
//PartialBasis::PartialBasis(PolynomialBasis& basis, Bsurface& surface) {
//	Ubasis = basis.Ubasis;
//	Vbasis = basis.Vbasis;
//	Ubasis_1 = do_partial(Ubasis); Vbasis_1 = do_partial(Vbasis);
//	Ubasis_2 = do_partial(Ubasis_1); Vbasis_2 = do_partial(Vbasis_1);
//	Uknot = surface.U;
//	Vknot = surface.V;
//	degree1 = surface.degree1;
//	degree2 = surface.degree2;
//}
PartialBasis::PartialBasis() {

}
PartialBasis::PartialBasis(Bsurface& surface) {
	init(surface);
}
void PartialBasis::init(Bsurface& surface) {
	PolynomialBasis pb(surface);
	Ubasis = pb.Ubasis;
	Vbasis = pb.Vbasis;
	Ubasis_1 = do_partial(Ubasis); Vbasis_1 = do_partial(Vbasis);
	Ubasis_2 = do_partial(Ubasis_1); Vbasis_2 = do_partial(Vbasis_1);
	Uknot = surface.U;
	Vknot = surface.V;
	degree1 = surface.degree1;
	degree2 = surface.degree2;
}
void PartialBasis::init(PolynomialBasis &pb){
	Ubasis = pb.Ubasis;
	Vbasis = pb.Vbasis;
	Ubasis_1 = do_partial(Ubasis); Vbasis_1 = do_partial(Vbasis);
	Ubasis_2 = do_partial(Ubasis_1); Vbasis_2 = do_partial(Vbasis_1);
	Uknot = pb.Uknot;
	Vknot = pb.Vknot;
	degree1 = pb.degree1;
	degree2 = pb.degree2;
}
void PartialBasis::clear() {
	Uknot.clear();
	Vknot.clear();
	Ubasis.clear();
	Vbasis.clear();
	Ubasis_1.clear();
	Vbasis_1.clear();
	Ubasis_2.clear();
	Vbasis_2.clear();
}
std::vector<double> PartialBasis::poly(const int id, const double value, const bool UVknot, int partial) {
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
		if (partial == 0) {
			result = Vbasis[which][index];
		}
		if (partial == 1) {
			result = Vbasis_1[which][index];
		}
		if (partial == 2) {
			result = Vbasis_2[which][index];
		}

	}
	else {
		if (partial == 0) {
			result = Ubasis[which][index];
		}
		if (partial == 1) {
			result = Ubasis_1[which][index];
		}
		if (partial == 2) {
			result = Ubasis_2[which][index];
		}
	}
	return result;
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
	std::vector<double> func = ply_operations::polynomial_times(func1, func2);
	//std::cout << "times" << std::endl;
	//print_vector(func);
	double upper = u2;
	if (u2 == U.back()) {
		upper = U.back() - SCALAR_ZERO;
	}
	double result = ply_operations::polynomial_integration(func, u1, upper);
	timer.stop();
	time2 += timer.getElapsedTimeInMilliSec();

	return result;
}
double construct_an_integration(const int degree, const std::vector<double>& U,
	const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2,
	PartialBasis &basis, const bool uv) {
	timer.start();
	//std::vector<double> func1 = basis.poly(i1, u1, uv);
	//std::vector<double> func2 = basis.poly(i2, u1, uv);	
	timer.stop();
	time0 += timer.getElapsedTimeInMilliSec();

	timer.start();
	std::vector<double> func1 = basis.poly(i1, u1, uv, partial1);
	std::vector<double> func2 = basis.poly(i2, u1, uv, partial2);
	timer.stop();
	time1 += timer.getElapsedTimeInMilliSec();

	timer.start();
	std::vector<double> func = ply_operations::polynomial_times(func1, func2);

	double upper = u2;
	if (u2 == U.back()) {
		upper = U.back() - SCALAR_ZERO;
	}
	double result = ply_operations::polynomial_integration(func, u1, upper);
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
	std::vector<double> func = ply_operations::polynomial_times(func1, func2);
	//std::cout << "times" << std::endl;
	//print_vector(func);
	double upper = u2;
	if (u2 == U.back()) {
		upper = U.back() - SCALAR_ZERO;
	}
	double result = ply_operations::polynomial_integration(func, u1, upper);
	timer.stop();
	time2 += timer.getElapsedTimeInMilliSec();

	return result;
}
// do partial difference to Pi, the cofficient of jth element Pj.
double surface_energy_least_square(Bsurface& surface, const int i, const int j, PartialBasis& basis) {
	// figure out which Pij corresponding to the ith control point
	int partial_i = i / (surface.nv() + 1);
	int partial_j = i - partial_i * (surface.nv() + 1);

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
			//std::cout << " k1, k2 " << k1 << " " << k2 << std::endl;
			if (coff_i<k1 - degree1 || coff_i>k1 || coff_j<k2 - degree2 || coff_j>k2) {
				continue;
			}
			assert(k1 + 1 < surface.U.size());
			assert(k2 + 1 < surface.V.size());
			if (surface.U[k1] == surface.U[k1 + 1] || surface.V[k2] == surface.V[k2 + 1]) {
				continue;// if the area is 0, then no need to compute
			}
			// Suu part
			double value1 = construct_an_integration(degree1, surface.U, 2, 2,
				partial_i, coff_i, surface.U[k1], surface.U[k1 + 1], basis, false);
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

// in interval [U[i], U[i+1])x[V[j], V[j+1])
double discrete_surface_partial_value_squared(const int partial1, const int partial2,
	const int i, const int j, Bsurface& surface,
	PartialBasis& basis, const double u, const double v) {
	int p = surface.degree1;
	int q = surface.degree2;
	Eigen::VectorXd Nl(p + 1);
	Eigen::VectorXd Nr(q + 1);
	for (int k = 0; k < p + 1; k++) {
		Nl[k] = ply_operations::polynomial_value(basis.poly(i - p + k, u, 0, partial1), u);
	}
	for (int k = 0; k < q + 1; k++) {
		Nr[k] = ply_operations::polynomial_value(basis.poly(j - q + k, v, 1, partial2), v);
	}
	Eigen::MatrixXd px(p + 1, q + 1), py(p + 1, q + 1), pz(p + 1, q + 1);
	for (int k1 = 0; k1 < p + 1; k1++) {
		for (int k2 = 0; k2 < q + 1; k2++) {
			px(k1, k2) = surface.control_points[i - p + k1][j - q + k2][0];
			py(k1, k2) = surface.control_points[i - p + k1][j - q + k2][1];
			pz(k1, k2) = surface.control_points[i - p + k1][j - q + k2][2];
		}
	}
	double x = (Nl.transpose()*px*Nr);
	double y = (Nl.transpose()*py*Nr);
	double z = (Nl.transpose()*pz*Nr);
	return x * x + y * y + z * z;
}

// calculate thin-plate-energy in region [Ui, U(i+1)]x[Vj, V(j+1)]
Eigen::MatrixXd Bsurface::surface_energy_calculation(Bsurface& surface, PartialBasis& basis,
	const int discrete, Eigen::MatrixXd &energy_uu, Eigen::MatrixXd &energy_vv, Eigen::MatrixXd& energy_uv) {
	int p = surface.degree1;
	int q = surface.degree2;
	std::vector<double> U = surface.U;
	std::vector<double> V = surface.V;
	int nu = surface.nu(); // U[p]=0, U[nu+1]=1
	int nv = surface.nv();
	int uint = nu + 1 - p; //nbr of u intervals
	int vint = nv + 1 - q;
	int n_sample = discrete + 2;// in each interval, there are discrete+2 sample points
	energy_uu.resize(uint, vint);
	energy_vv.resize(uint, vint);
	energy_uv.resize(uint, vint);
	for (int i = 0; i < uint; i++) {
		for (int j = 0; j < vint; j++) {
			double u0 = U[i + p];
			double u1 = U[i + p + 1];
			double v0 = V[j + q];
			double v1 = V[j + q + 1];
			double delta_u = (u1 - u0) / (n_sample - 1);
			double delta_v = (v1 - v0) / (n_sample - 1);
			int Ni = i + p;// interval is [U[Ni], U[Ni+1])
			int Nj = j + q;
			Eigen::MatrixXd values_uu(n_sample, n_sample);
			Eigen::MatrixXd values_uv(n_sample, n_sample);
			Eigen::MatrixXd values_vv(n_sample, n_sample);
			//std::cout << "u0, u1, v0, v1 " << u0 << " " << u1 << " " << v0 << " " << v1 << std::endl;
			for (int k1 = 0; k1 < n_sample; k1++) {
				//std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
				for (int k2 = 0; k2 < n_sample; k2++) {
					//int Ni = k1 == n_sample - 1 ? i + p + 1 : i + p;// if select the last point, 
					double uratio = k1 < n_sample - 1 ? 1 : 1 - SCALAR_ZERO;
					double vratio = k2 < n_sample - 1 ? 1 : 1 - SCALAR_ZERO;
					double uvalue = u0 + delta_u * k1*uratio;
					double vvalue = v0 + delta_v * k2*vratio;
					//std::cout << "uv values " << uvalue << " " << vvalue << std::endl;
					// Suu
					values_uu(k1, k2) = discrete_surface_partial_value_squared(2, 0, Ni, Nj, surface, basis, uvalue, vvalue);
					//Svv
					values_vv(k1, k2) = discrete_surface_partial_value_squared(0, 2, Ni, Nj, surface, basis, uvalue, vvalue);
					//Suv
					values_uv(k1, k2) = discrete_surface_partial_value_squared(1, 1, Ni, Nj, surface, basis, uvalue, vvalue);

				}
			}
			double uusum = 0;
			double vvsum = 0;
			double uvsum = 0;
			double single_area = delta_u * delta_v;
			for (int k1 = 0; k1 < n_sample - 1; k1++) {
				for (int k2 = 0; k2 < n_sample - 1; k2++) {
					uusum += (values_uu(k1, k2) + values_uu(k1, k2 + 1) + values_uu(k1 + 1, k2) + values_uu(k1 + 1, k2 + 1)) / 4;
					vvsum += (values_vv(k1, k2) + values_vv(k1, k2 + 1) + values_vv(k1 + 1, k2) + values_vv(k1 + 1, k2 + 1)) / 4;
					uvsum += (values_uv(k1, k2) + values_uv(k1, k2 + 1) + values_uv(k1 + 1, k2) + values_uv(k1 + 1, k2 + 1)) / 4;

				}
			}
			energy_uu(i, j) = uusum * single_area;
			energy_vv(i, j) = vvsum * single_area;
			energy_uv(i, j) = uvsum * single_area;
		}
	}
	Eigen::MatrixXd energy(uint, vint);
	energy = energy_uu + 2 * energy_uv + energy_vv;
	return energy;
	std::cout << "energy\n" << energy << std::endl;

}

// [U[which],U[which+1]) is the problematic one
void Bsurface::detect_max_energy_interval(Bsurface& surface, const Eigen::MatrixXd& energy, const Eigen::MatrixXd &energy_uu,
	const Eigen::MatrixXd & energy_vv, bool& uorv, int &which, double& em) {
	int k1, k2;
	em = 0;
	for (int i = 0; i < energy.rows(); i++) {
		for (int j = 0; j < energy.cols(); j++) {
			double evalue = energy(i, j);
			if (evalue > em) {
				em = evalue;
				k1 = i;
				k2 = j;
			}
		}
	}
	std::cout << "max energy " << em << std::endl;
	int p = surface.degree1;
	int q = surface.degree2;
	int u_interval = k1 + p;
	int v_interval = k2 + q;
	if (energy_uu(k1, k2) > energy_vv(k1, k2)) {
		uorv = 1;
		which = v_interval;
	}
	else {
		which = u_interval;
		uorv = 0;
	}
	return;
}


// which_part = 0: Suu; which_part = 1, Suv; which_part = 2, Svv.
Eigen::MatrixXd energy_part_of_surface_least_square(Bsurface& surface, PartialBasis& basis) {
	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	Eigen::MatrixXd result(psize, psize);
	for (int i = 0; i < psize; i++) {
		//std::cout << "the ith row of matrix" << std::endl;
		for (int j = 0; j < psize; j++) {
			result(i, j) = surface_energy_least_square(surface, i, j, basis);
		}
	}
	std::cout << "energy matrix finish calculation" << std::endl;
	return result;
}

double surface_error_least_square(Bsurface& surface, const int i, const int j,
	const Eigen::MatrixXd& paras) {
	// figure out which Pij corresponding to the ith control point
	int partial_i = i / (surface.nv() + 1);
	int partial_j = i - partial_i * (surface.nv() + 1);

	// figure out which Pij corresponding to the jth control point
	int coff_i = j / (surface.nv() + 1);
	int coff_j = j - coff_i * (surface.nv() + 1);

	// if do partial Pij, the related other Pi1j1 will be: i1 = i-p~i+p, j1= j-q~j+q
	int degree1 = surface.degree1, degree2 = surface.degree2;
	/*if (coff_i<partial_i - degree1 || coff_i>partial_i + degree1) return 0;
	if (coff_j<partial_j - degree2 || coff_j>partial_j + degree2) return 0;*/


	// do partial Pij
	// for each block, (U_k1, U_(k1+1)), (V_k2, V_(k2+1)). the related control points are k1-p,...,k1 and k2-q,..., k2
	double result = 0;
	for (int k1 = 0; k1 < paras.rows(); k1++) {
		double u = paras(k1, 0);
		double v = paras(k1, 1);
		double N1 = Nip(partial_i, degree1, u, surface.U);
		double N2 = Nip(partial_j, degree2, v, surface.V);
		double N3 = Nip(coff_i, degree1, u, surface.U);
		double N4 = Nip(coff_j, degree2, v, surface.V);
		if (N1 == 0 || N2 == 0 || N3 == 0 || N4 == 0) {
			continue;
		}
		result += N1 * N2*N3*N4;
	}
	return result;

}
Eigen::MatrixXd error_part_of_surface_least_square(Bsurface& surface, const Eigen::MatrixXd& paras) {
	// figure out which Pij corr) {
	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	Eigen::MatrixXd result(psize, psize);
	for (int i = 0; i < psize; i++) {
		//std::cout << "the ith row of matrix" << std::endl;
		for (int j = 0; j < psize; j++) {
			result(i, j) = surface_error_least_square(surface, i, j, paras);
		}
	}
	std::cout << "error matrix finish calculation" << std::endl;
	return result;
}
Eigen::VectorXd right_part_of_least_square_approximation(Bsurface& surface, const Eigen::MatrixXd& paras,
	const Eigen::MatrixXd& ver, const int dimension) {

	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	Eigen::VectorXd result(psize);
	for (int i = 0; i < psize; i++) {
		double res = 0;
		for (int j = 0; j < paras.rows(); j++) {
			int partial_i = i / (surface.nv() + 1);
			int partial_j = i - partial_i * (surface.nv() + 1);
			int degree1 = surface.degree1, degree2 = surface.degree2;

			double u = paras(j, 0);
			double v = paras(j, 1);
			double N1 = Nip(partial_i, degree1, u, surface.U);
			double N2 = Nip(partial_j, degree2, v, surface.V);
			if (N1 == 0 || N2 == 0) {
				continue;
			}
			res += N1 * N2*ver(j, dimension);
		}
		result(i) = res;
	}
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

Eigen::MatrixXd surface_least_square_lambda_multiplier_left_part(Bsurface& surface,
	const Eigen::MatrixXd& paras, PartialBasis& basis) {
	std::cout << "inside left part" << std::endl;
	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	int target_size = paras.rows();// nbr of target data points
	int size = psize + target_size;
	Eigen::MatrixXd result(size, size);
	Eigen::MatrixXd rd = Eigen::MatrixXd::Zero(target_size, target_size);// right down corner part
	std::cout << "finish rd" << std::endl;
	Eigen::MatrixXd lu = energy_part_of_surface_least_square(surface, basis);
	std::cout << "finish lu" << std::endl;
	Eigen::MatrixXd ld = eqality_part_of_surface_least_square(surface, paras);
	std::cout << "finish ld" << std::endl;
	Eigen::MatrixXd ru = -ld.transpose();
	std::cout << "finish ru" << std::endl;
	//lambda_part_of_surface_least_square(surface, paras);

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
void Bsurface::solve_control_points_for_fairing_surface(Bsurface& surface, const Eigen::MatrixXd& paras,
	const Eigen::MatrixXd & points, PartialBasis& basis) {
	//using namespace Eigen;
	typedef Eigen::SparseMatrix<double> SparseMatrixXd;
	assert(paras.rows() == points.rows());
	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	std::vector<Vector3d> cps(psize);// control points
	Eigen::MatrixXd A;
	Eigen::FullPivLU<Eigen::DenseBase<Eigen::MatrixXd>::PlainMatrix> decomp;
	A = surface_least_square_lambda_multiplier_left_part(surface, paras, basis);
	SparseMatrixXd matB;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

	matB = A.sparseView();
	A.resize(0, 0);
	solver.compute(matB);
	if (solver.info() != Eigen::Success) {
		// decomposition failed
		std::cout << "solving failed" << std::endl;
		return;
	}

	for (int i = 0; i < 3; i++) {
		Eigen::MatrixXd b = surface_least_square_lambda_multiplier_right_part(surface, paras, points, i);
		double err = 0.0;

		// solve the matrix contains the p and lambda
		std::cout << "before solving" << std::endl;
		Eigen::MatrixXd p_lambda;

		p_lambda = solver.solve(b);
		if (solver.info() != Eigen::Success) {
			std::cout << "solving failed" << std::endl;
			return;
		}



		double relative_error = (matB*p_lambda - b).norm() / b.norm(); // norm() is L2 norm
		std::cout << "after solving, error is " << relative_error << std::endl;
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
}

