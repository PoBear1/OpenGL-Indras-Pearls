//
//  complex_calc.hpp
//
//  Created by Polar Xiong on 11/12/2025.
//

#ifndef complex_calc_hpp
#define complex_calc_hpp

// the standard
#include <exception>
#include <iostream>
#include <complex>
#include <limits>
#include <array>

// gmp/mpfr headers
#include <mp++/mp++.hpp>

// For simplicity's sake!
typedef mppp::complex complex;
// First typedef the matrix construct so the construct will then become
// {
// 		{a[0][0],a[0][1]},
// 		{a[1][0],a[1][1]}
// }
constexpr int prec_N = 1024;
typedef std::array<std::array<complex,2>,2> matrix;
// uhhhhh so um apparently i uh need to compute vector stuff as well??? unlucky
// we typedef vector the obv way but give it c_vector name so it doesn't clash with std::vector
typedef std::array<complex,2> c_vector;
// no-one ever told me to also need to include vectors UNTIL NOW
c_vector O_v={mppp::real("0", prec_N), mppp::real("0", prec_N)};
mppp::real norm(c_vector a) {
	return norm(a[0])+norm(a[1]);
}

// alright i'm annoyed here's an actual function to print matrices
void print(matrix T) {
	for(auto i:T) {
		for(auto j:i) {
			std::cout<<j<<" ";
		}
		std::cout<<"\n";
	}
	std::cout<<"\n";
}


// right so we do matrix stuff
matrix I={{{mppp::real("1", prec_N), mppp::real("0", prec_N)}, {mppp::real("0", prec_N), mppp::real("1", prec_N)}}};
matrix O={{{mppp::real("0", prec_N), mppp::real("0", prec_N)}, {mppp::real("0", prec_N), mppp::real("0", prec_N)}}};
complex det(matrix a) {
	return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}
complex tr(matrix a) {
	return a[0][0] + a[1][1];
}
matrix inv(matrix a) {
	if(det(a).zero_p()) {
		throw std::runtime_error("Error when calling \"matrix inv(matrix)\": attempting to invert matrix with determinant zero");
	} else {
		complex d = det(a);
		return {{{a[1][1]/d, -a[0][1]/d}, {-a[1][0]/d, a[0][0]/d}}};
	}
}
const matrix operator*(const matrix a,const matrix b) {
	return {{
		{a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]},
		{a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]}
	}};
}
matrix operator/(matrix a,matrix b) {return a * inv(b);}
// no one ever told me i forgor about multiplying by constants!
matrix operator*(complex a,matrix b) {
	return {{
		{a * b[0][0], a * b[0][1]},
		{a * b[1][0], a * b[1][1]}
	}};
}
// and adding and subtracting!!!
matrix operator+(matrix a,matrix b) {
	return {{
		{a[0][0] + b[0][0], b[0][1] + a[0][1]},
		{a[1][0] + b[1][0], a[1][1] + b[1][1]}
	}};
}
matrix operator-(matrix b) {
	return {{
		{-b[0][0], -b[0][1]},
		{-b[1][0], -b[1][1]}
	}};
}
matrix operator-(matrix a,matrix b) {
	return a + (-b);
}
// and applying matrices to vectors!
c_vector operator*(matrix a,c_vector b) {
	return {
		a[0][0] * b[0] + a[0][1] * b[1], a[1][0] * b[0] + a[1][1] * b[1]
	};
}
matrix sl2(matrix a) {
	if(!det(a).zero_p()) {
		complex d = det(a);
		return {{{a[0][0]/sqrt(d), a[0][1]/sqrt(d)}, {a[1][0]/sqrt(d), a[1][1]/sqrt(d)}}};
	} else {throw std::runtime_error("Attempting to normalise degenerate matrix");}
}

// assuming you've sanitised the values by normalising already
std::array<complex,2> eigenvalues(matrix a) {
	return {
		(tr(a) + sqrt(tr(a) * tr(a) - complex(mppp::real{"4", prec_N})))/(complex(mppp::real{"2", prec_N})),
		(tr(a) - sqrt(tr(a) * tr(a) - complex(mppp::real{"4", prec_N})))/(complex(mppp::real{"2", prec_N}))
	};
}

// mobius maps
complex operator*(matrix a,complex z) {
	return proj(a[0][0] * z + a[0][1])/(a[1][0] * z + a[1][1]);
}
std::array<complex,2> fixed_pts(matrix a) {
	if(!a[1][0].zero_p()) {
		matrix b = sl2(a);
		return {
			(b[0][0] - b[1][1] + sqrt(tr(b) * tr(b) - complex(mppp::real{"4", prec_N})))/(complex(mppp::real{"2", prec_N}) * b[1][0]),
			(b[0][0] - b[1][1] - sqrt(tr(b) * tr(b) - complex(mppp::real{"4", prec_N})))/(complex(mppp::real{"2", prec_N}) * b[1][0])
		};
	} else {
		return {proj(complex(mppp::real{"1", prec_N})/complex(mppp::real{"0", prec_N})), proj(a[0][1]/(a[1][1] - a[0][0]))};
	}
}
complex attract_fixed(matrix a) {
	if(eigenvalues(a)[0] == eigenvalues(a)[1]) {
		if(!a[1][0].zero_p()) {
			return (a[0][0] - a[1][1]) / (complex(mppp::real{"2", prec_N}) * a[1][0]);
		} else {
			return proj(complex(mppp::real{"1", prec_N}, mppp::real{"0", prec_N}));
		}
	} else {
		complex large = (norm(eigenvalues(a)[0]) > norm(eigenvalues(a)[1])) ? eigenvalues(a)[0] : eigenvalues(a)[1];
		c_vector large_v = (norm({a[0][0] - large, a[0][1]}) > norm({a[1][0], a[1][1] - large})) ? (c_vector({a[0][0] - large, a[0][1]})) : (c_vector({a[1][0], a[1][1] - large}));
		return proj(-large_v[1]/large_v[0]);
	}
}

matrix T(matrix a) {matrix b{a}; swap(b[1][0],b[0][1]); return b;}

// Geometric construct - attempting to consider lines and circles simultaneously which will be hard.
class cline {
private:
	// ok so what do i want...
	complex a; mppp::real b; bool is_line_v;
	// so for circles, a is centre, b is radius
	// for lines with representation 2Re(az)=b, a and b obv.
	complex shoelace(complex a,complex b,complex c) {return (b*conj(c)-conj(b)*c)+(c*conj(a)-a*conj(c))+(a*conj(b)-b*conj(a));}
	complex circumcentre(complex a,complex b,complex c) {
		complex z=shoelace(a,b,c);
		if(z==complex(mppp::real{"0", prec_N}, mppp::real{"0", prec_N})) {
			return complex(mppp::real{"1", prec_N}, mppp::real{"0", prec_N})/complex(mppp::real{"0", prec_N}, mppp::real{"0", prec_N});
		} else {
			complex z1=(a*b*conj(b)-b*a*conj(a))+(b*c*conj(c)-c*b*conj(b))+(c*a*conj(a)-a*c*conj(c));
			return z1/z;
		}
	}
	mppp::real area(mppp::real A, mppp::real B, mppp::real C) {
		mppp::real x=(B + C - A)/mppp::real{"2", prec_N};
		mppp::real y=(C + A - B)/mppp::real{"2", prec_N};
		mppp::real z=(A + B - C)/mppp::real{"2", prec_N};
		return sqrt((x + y + z) * x * y * z);
	}
	mppp::real circumradius(complex a, complex b, complex c) {
		mppp::real A = abs(b - c), B = abs(c - a),C = abs(a - b);
		return (A * B * C)/(mppp::real{"4", prec_N} * area(A, B, C));
	}
	complex line_direction(complex a, complex b) {return conj(a)-conj(b);}
	mppp::real line_offset(complex a, complex b) {return get_real_imag(a*conj(b)-b*conj(a)).first;}
public:
	cline(complex a_in, mppp::real b_in, bool is_line_v_in=false): a(a_in), b(b_in), is_line_v(is_line_v_in) {}
	bool is_line() {return is_line_v;}
	bool is_circle() {return !is_line_v;}
	complex centre() {return a;} mppp::real radius() {return b;}
	complex dir() {return a;} mppp::real offset() {return b;}
	cline mobius(matrix T) {
		complex z1, z2, z3;
		if(!this->is_line_v) {
			// complex omega=complex(1,1)/complex(abs(complex(1,1)));
			complex omega = complex(mppp::real{"-0.5", prec_N}, mppp::real{"0.5", prec_N} * sqrt(mppp::real{"3", prec_N}));
			z1 = a + complex(b);
			z2 = a + complex(b) * omega;
			z3 = a + complex(b) * omega * omega;
		} else {
			z1 = complex(b)/a;
			z2 = (complex(b)/a)*complex(1,1);
			z3 = (complex(b)/a)*complex(1,-1);
			if(b.zero_p()) {
				z2 = complex(mppp::real{"0", prec_N}, mppp::real{"1", prec_N})/a;
				z3 = complex(mppp::real{"0", prec_N}, mppp::real{"-1", prec_N})/a;
			} else {}
		}
		z1 = T * z1;
		z2 = T * z2;
		z3 = T * z3;
		z1 *= complex(mppp::real{"10000000", prec_N});
		z2 *= complex(mppp::real{"10000000", prec_N});
		z3 *= complex(mppp::real{"10000000", prec_N});
		if(shoelace(z1,z2,z3)!=complex(0)) {
			return cline(circumcentre(z1,z2,z3)/complex(mppp::real{"10000000", prec_N}),
						 circumradius(z1,z2,z3)/mppp::real{"10000000", prec_N});
		} else {
			return cline(line_direction(z1,z3)/complex(mppp::real{"10000000", prec_N}),
						 line_offset(z1,z3)/mppp::real{"10000000", prec_N}, true);
		}
	}
	bool contains(complex z) {
		if(this->is_line_v) {
			return false;
		} else {
			return abs(z-this->centre())<=this->radius();
		}
	}
};

#endif /* complex_calc_hpp */
