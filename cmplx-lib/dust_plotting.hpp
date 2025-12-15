//
//  dust_plotting.hpp
//  fractal-drawing
//
//  Created by Polar Xiong on 11/12/2025.
//

#ifndef dust_plotting_hpp
#define dust_plotting_hpp

#include "complex_calc.hpp"
#include "matrix_parameter.hpp"
#include <iostream>
#include <vector>

class plotting {
public:
	plotting(int levm, complex ta, complex tb, mppp::real eps, std::vector<std::vector<std::vector<int>>> repeats = {{{0}}, {{1}}, {{2}}, {{3}}}, bool x = true) : levmax(levm), epsilon(eps) {
		word.resize(levmax + 1, I);
		tags.resize(levmax + 1);
		matrix a = matrix_generation::grandma_recipe(ta, tb, x).first, b = matrix_generation::grandma_recipe(ta, tb, x).second;
		gens[0] = a; gens[1] = b; gens[2] = inv(a); gens[3] = inv(b);
		tags[0] = 0;
		word[0] = gens[0];
		fix.resize(4);
		for(int i = 0; i < 4; ++i) {fix[i].resize(repeats[i].size());}
		for(int i = 0; i < 4; ++i) {
			for(int j = 0; j < repeats[i].size(); ++j) {
				matrix M = I;
				for(int k : repeats[i][j]) {
					M = M * gens[k];
				}
				fix[i][j] = attract_fixed(M);
			}
		}
	}
	std::vector<complex> dust_points_list() {
		std::vector<complex> points;
		do {
			while(!termination()) {go_forward();}
			points.push_back(newpoint);
			first = false;
			oldpoint = newpoint;
			do {lev--;} while(lev != -1 && tags[lev + 1] == (tags[lev] + 3)%4);
			tags[lev + 1] = (tags[lev + 1] + 3)%4;
			if(lev == -1) {
				word[0] = gens[tags[0]];
			} else {
				word[lev + 1] = word[lev] * gens[tags[lev + 1]];
			}
			lev++;
		} while(lev != 0 || tags[0] != 0);
		return points;
	}
	std::vector<complex> dust_points_list_buffed() {
		tags[0] = 0;
		word[0] = gens[0];
		std::vector<complex> points;
		do {
			while(!buffed_termination()) {go_forward();}
			for(complex z : fix[tags[lev]]) {
				points.push_back(word[lev] * z);
			}
			first = false;
			oldpoint = word[lev] * fix[tags[lev]].back();
			do {lev--;} while(lev != -1 && tags[lev + 1] == (tags[lev] + 3)%4);
			tags[lev + 1] = (tags[lev + 1] + 3)%4;
			if(lev == -1) {
				word[0] = gens[tags[0]];
			} else {
				word[lev + 1] = word[lev] * gens[tags[lev + 1]];
			}
			lev++;
		} while(lev != 0 || tags[0] != 0);
		return points;
	}
private:
	int levmax, lev = 0;
	matrix gens[4];
	std::vector<std::vector<complex>> fix;
	std::vector<int> tags;
	std::vector<matrix> word;
	std::vector<std::vector<matrix>> repeated;
	mppp::real epsilon;
	complex oldpoint, newpoint;
	bool first = true;
	void go_forward() {
		lev++;
		tags[lev] = (tags[lev - 1] + 1) % 4;
		word[lev] = word[lev - 1] * gens[tags[lev]];
	}
	bool termination() {
		newpoint = word[lev] * fix[tags[lev]][0];
		if(!first) {
			return (lev == levmax || abs(newpoint - oldpoint) < epsilon);
		} else {
			return lev == levmax;
		}
	}
	bool buffed_termination() {
		newpoint = word[lev] * fix[tags[lev]][0];
		bool combined{true};
		for(int i = 1; i < fix[tags[lev]].size(); ++i) {
			combined = combined & (abs((word[lev] * fix[tags[lev]][i]) - (word[lev] * fix[tags[lev]][i - 1])) < epsilon);
		}
		if(!first) {
			return (lev == levmax || (combined && abs(newpoint - oldpoint) < epsilon));
		} else {
			return lev == levmax || combined;
		}
	}
};

#endif /* dust_plotting_hpp */
