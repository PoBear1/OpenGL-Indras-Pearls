// Standard stuff
#define GLM_ENABLE_EXPERIMENTAL
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

// OpenGL stuff
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
using namespace glm;

// GMP headers
#include <gmpxx.h>
#include <gmp.h>

// MPFR headers
#ifndef MPFR_USE_INTMAX_T
#define MPFR_USE_INTMAX_T
#endif
#ifndef MPFR_USE_NO_MACRO
#define MPFR_USE_NO_MACRO
#endif
#include <mpfr.h>
#include "../cmplx-lib/mpreal.h"
using mpfr::mpreal;


// local lib
#include "load_shader.hpp"
#include "../cmplx-lib/complex_calc.hpp"
#include "../cmplx-lib/matrix_parameter.hpp"
#include "../cmplx-lib/dust_plotting.hpp"

glm::vec3 hsv_to_rgb(glm::vec3 hsv) {
	glm::vec3 a(0,0,0);
	float H = hsv.x, S = hsv.y, V = hsv.z;
	float C = V * S, X = C * (1 - fabs(fmodf(H/60, (float)2) - 1)), m = V - C;
	if(0 <= H && H < 60) {
		a = glm::vec3(C, X, 0);
	} else if(60 <= H && H < 120) {
		a = glm::vec3(X, C, 0);
	} else if(120 <= H && H < 180) {
		a = glm::vec3(0, C, X);
	} else if(180 <= H && H < 240) {
		a = glm::vec3(0, X, C);
	} else if(240 <= H && H < 300) {
		a = glm::vec3(X, 0, C);
	} else {
		a = glm::vec3(C, 0, X);
	}
	a = a + glm::vec3(m, m, m);
	return a;
}

int main() {
	printf ("MPFR library: %-12s\nMPFR header:  %s (based on %d.%d.%d)\n",
			mpfr_get_version (), MPFR_VERSION_STRING, MPFR_VERSION_MAJOR,
			MPFR_VERSION_MINOR, MPFR_VERSION_PATCHLEVEL);
	mpreal::set_default_prec(mpfr::digits2bits(50));
	if(!glfwInit()) {
		fprintf(stderr, "Failed to initialize GLFW\n");
		return -1;
	} else {}
	
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	
	// GLFWmonitor* primary = glfwGetPrimaryMonitor();
	// const GLFWvidmode *mode = glfwGetVideoMode(primary);
	GLFWwindow* window = glfwCreateWindow(800, 600, "yes us", NULL, NULL);
	
	if(window == NULL) {
		fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
		glfwTerminate();
		return -1;
	} else {}
	glfwMakeContextCurrent(window);
	
	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		return -1;
	} else {}
	
	GLuint programID = LoadShaders("shaders/vertex_shader.vtxshd", "shaders/fragment_shader.fragshd");
	if(programID == -1) {
		fprintf(stderr, "Shaders failed to compile/load, aborting...\n");
		return -1;
	} else {}
	
	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);
	
	complex ta = complex(mppp::real{"1.64213876", prec_N}, mppp::real{"-0.76658841", prec_N}), tb = complex(mppp::real{"2.00", prec_N}, mppp::real{"0.00", prec_N});
	std::vector<std::vector<std::vector<int>>> repetends = {
		{ {1, 2, 3, 0}, {0, 0, 3, 0, 0, 3, 0}, {0, 3, 0, 0, 0, 3, 0}, {0, 3, 0, 0, 3, 0, 0}, {3, 0, 0, 0, 3, 0, 0}, {3, 0, 0, 3, 0, 0, 0}, {3, 2, 1, 0} },
		{ {2, 3, 0, 1}, {2, 2, 2, 1, 2, 2, 1}, {2, 2, 1, 2, 2, 2, 1}, {1}, {0, 3, 2, 1} },
		{ {3, 0, 1, 2}, {2, 2, 1, 2, 2, 1, 2}, {2, 1, 2, 2, 2, 1, 2}, {2, 1, 2, 2, 1, 2, 2}, {1, 2, 2, 2, 1, 2, 2}, {1, 2, 2, 1, 2, 2, 2}, {1, 0, 3, 2} },
		{ {0, 1, 2, 3}, {0, 0, 0, 3, 0, 0, 3}, {0, 0, 3, 0, 0, 0, 3}, {3}, {2, 1, 0, 3} }
	};
	plotting help(100, ta, tb, mppp::real{"0.005", prec_N}, repetends);
	std::vector<complex> v = help.dust_points_list_buffed();
	int N = (int)v.size();
	
	GLfloat *g_vertex_buffer_data = (GLfloat*)calloc(2 * N, sizeof(GLfloat));
	GLfloat *g_color_buffer_data = (GLfloat*)calloc(3 * N, sizeof(GLfloat));
	glm::vec3 useless;
	float max_hor{0}, max_ver{0};
	for(int i = 0; i < N; ++i) {
		g_vertex_buffer_data[2*i] = (float)(mppp::get_real_imag(v[i]).first);
		g_vertex_buffer_data[2*i + 1] = (float)(mppp::get_real_imag(v[i]).second);
		max_hor = max(max_hor, g_vertex_buffer_data[2*i]);
		max_ver = max(max_ver, g_vertex_buffer_data[2*i + 1]);
		useless = hsv_to_rgb(vec3(((float)i)/N * 360, 1.0f, 1.0f));
		g_color_buffer_data[3*i] = useless.x;
		g_color_buffer_data[3*i + 1] = useless.y;
		g_color_buffer_data[3*i + 2] = useless.z;
	}
	
	std::cout<< "N: " << N << std::endl;
	
	GLuint vertexbuffer;
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, 2 * N * sizeof(GLfloat), g_vertex_buffer_data, GL_STATIC_DRAW);
	
	GLuint colorbuffer;
	glGenBuffers(1, &colorbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
	glBufferData(GL_ARRAY_BUFFER, 3 * N * sizeof(GLfloat), g_color_buffer_data, GL_STATIC_DRAW);
	
	glfwSetWindowTitle(window, "Playground");
	printf("Hell YES I've done it");
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GLFW_TRUE);
	
	int k=1, height, width;
	float reduce_size = 1;
	do {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glClearColor(0.0f, 0.0f, 0.3f, 0.0f);
		k++; k%=N;
		GLint ratio_id = glGetUniformLocation(programID, "ratio");
		GLint reduce_size_id = glGetUniformLocation(programID, "reduce_size");
		
		glfwGetWindowSize(window, &width, &height);
		float ratio = ((float)height)/width;
//		if(ratio > 1) {
//			reduce_size = max(max_hor / ratio, max_ver);
//		} else {
//			reduce_size = max(max_hor, max_ver * ratio);
//		}
		glUseProgram(programID);
		glUniform1f(ratio_id, ratio);
		glUniform1f(reduce_size_id, reduce_size);
		
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
		glDrawArrays(GL_LINE_STRIP, 0, N);
		glDisableVertexAttribArray(0);
		
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		
		glfwSwapBuffers(window);
		glfwPollEvents();
		
	}
	while( glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		  !glfwWindowShouldClose(window) );
	
	glfwDestroyWindow(window);
	glfwTerminate();
	
	return 0;
}

