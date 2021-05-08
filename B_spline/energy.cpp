#include<energy.h>
#include<surface.h> 
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

// construct an integration of multiplication of two B-spline basis (intergration of partial(Ni1)*partial(Ni2))
// the integration domain is [u1, u2]
double construct_an_integration(const int degree, const std::vector<double>& U,
	const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2) {
	
	std::vector<double> func1 = Nip_func(i1, degree, u1, U);
	std::vector<double> func2 = Nip_func(i2, degree, u1, U);
	//std::cout << "degree, "<<degree << std::endl;
	//std::cout << "func1 and func2" << std::endl;
	//print_vector(func1);
	//print_vector(func2);
	func1 = polynomial_differential(func1, partial1);
	func2 = polynomial_differential(func2, partial2);
	//std::cout << "differencial" << std::endl;
	//print_vector(func1);
	//print_vector(func2);
	std::vector<double> func = polynomial_times(func1, func2);
	//std::cout << "times" << std::endl;
	//print_vector(func);
	double upper = u2;
	if (u2 == U.back()) {
		upper = U.back() - SCALAR_ZERO;
	}
	double result = polynomial_integration(func, u1, upper);
	std::cout << "integration range, " <<u1<<","<<upper<<", result, "<<result<< std::endl;

	return result;
}

// do partial difference to Pi, the cofficient of jth element Pj.
double energy_element_value( Bsurface& surface, const int i, const int j) {
	// locate Pij
	int partial_i = i / (surface.nv()+1);
	int partial_j = i - partial_i * (surface.nv()+1);

	// locate the cofficient of pij
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
				partial_i, coff_i, surface.U[k1], surface.U[k1 + 1]);
			double value2 = construct_an_integration(degree2, surface.V, 0, 0, 
				partial_j, coff_j, surface.V[k2], surface.V[k2 + 1]);
			// Suv part
			double value3 = construct_an_integration(degree1, surface.U, 1, 1,
				partial_i, coff_i, surface.U[k1], surface.U[k1 + 1]);
			double value4 = construct_an_integration(degree2, surface.V, 1, 1,
				partial_j, coff_j, surface.V[k2], surface.V[k2 + 1]);
			// Svv part
			double value5 = construct_an_integration(degree1, surface.U, 0, 0,
				partial_i, coff_i, surface.U[k1], surface.U[k1 + 1]);
			double value6 = construct_an_integration(degree2, surface.V, 2, 2,
				partial_j, coff_j, surface.V[k2], surface.V[k2 + 1]);
			result += 2 * value1*value2 + 4 * value3*value4 + 2 * value5*value6;
		}
	}
	return result;
	
}


// which_part = 0: Suu; which_part = 1, Suv; which_part = 2, Svv.
Eigen::MatrixXd energy_part_of_solving_control_points(Bsurface& surface) {
	int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
	Eigen::MatrixXd result(psize, psize);
	for (int i = 0; i < psize; i++) {
		for (int j = 0; j < psize; j++) {
			result(i, j) = energy_element_value(surface, i, j);
		}
	}
	return result;
}