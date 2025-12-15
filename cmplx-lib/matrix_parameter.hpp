//
//  matrix_parameter.hpp
//  fractal-hell
//
//  Created by Polar Xiong on 11/12/2025.
//

#ifndef matrix_parameter_hpp
#define matrix_parameter_hpp

#include "complex_calc.hpp"
#include <utility>
#include <tuple>
class matrix_generation {
public:
	static std::pair<matrix, matrix> grandma_recipe(complex ta, complex tb, bool x = true) {
		complex tab = (ta * tb - (2 * x - 1) * sqrt(ta * ta * tb * tb - complex(mppp::real{"4", prec_N}) * (ta * ta + tb * tb))) / complex(mppp::real{"2", prec_N});
		complex z0 = ((tab - complex(mppp::real{"2", prec_N})) * tb)/(tb * tab - ta - ta + complex(mppp::real{"0", prec_N}, mppp::real{"2", prec_N}) * tab);
		matrix b{{
			{(tb - complex(mppp::real{"0", prec_N}, mppp::real{"2", prec_N}))/complex(mppp::real{"2", prec_N}), tb/complex(mppp::real{"2", prec_N})},
			{tb/complex(mppp::real{"2", prec_N}), (tb + complex(mppp::real{"0", prec_N}, mppp::real{"2", prec_N}))/complex(mppp::real{"2", prec_N})}
		}},
		ab{{
			{tab/complex(mppp::real{"2", prec_N}), (tab - complex(mppp::real{"2", prec_N}))/(z0 + z0)},
			{(tab + complex(mppp::real{"2", prec_N})) * z0/complex(mppp::real{"2", prec_N}), tab/complex(mppp::real{"2", prec_N})}
		}};
		return {ab * inv(b),b};
	}
	// static std::pair<matrix,matrix> jorgensen_recipe(complex ta,complex tb) {
	// 	complex tab=(ta*tb+sqrt(ta*ta*tb*tb-complex(4,0)*(ta*ta+tb*tb)))/complex(2,0);
	// 	matrix
	// }
	static std::pair<matrix,matrix> three_trace_recipe(complex ta, complex tb, complex tab) {
		complex tc = ta * ta + tb * tb + tab * tab - ta * tb * tab - complex(mppp::real{"2", prec_N});
		complex Q=sqrt(complex(mppp::real{"2", prec_N}) - tc),R;
		if(abs(tc + complex(mppp::real{"0", prec_N}, mppp::real{"1", prec_N}) * Q * sqrt(tc + complex(mppp::real{"2", prec_N}))) >= mppp::real{"2", prec_N}) {
			R = sqrt(tc + complex(mppp::real{"2", prec_N}));
		} else {
			R = -sqrt(tc + complex(mppp::real{"2", prec_N}));
		}
		complex z0=(tab - complex(mppp::real{"2", prec_N})) * (tb + R)/(tb * tab - ta - ta + complex(mppp::real{"0", prec_N}, mppp::real{"1", prec_N}) * Q * tab);
		matrix a{{
			{ta/complex(mppp::real{"2", prec_N}), (ta * tab - complex(mppp::real{"2", prec_N}) * tb + complex(mppp::real{"0", prec_N}, mppp::real{"2", prec_N}) * Q)/((tab + tab + complex(mppp::real{"4", prec_N})) * z0)},
			{((ta*tab-complex(2,0)*tb-complex(0.0,2.0)*Q)*z0)/((tab+tab-complex(4,0))),ta*complex(0.5,0)}
		}},
		b{{
			{(tb-complex(0,1)*Q)*complex(0.5,0.0),(tb*tab-ta-ta-complex(0,1)*Q*tab)/((tab+tab+complex(4,0))*z0)},
			{(tb*tab-ta-ta+complex(0,1)*Q*tab)*z0/(tab+tab-complex(4,0)),(tb+complex(0.0,1.0)*Q)*complex(0.5,0.0)}
		}};
		return {a,b};
	}
};

#endif /* matrix_parameter_hpp */
