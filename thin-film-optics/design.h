#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <random>


enum layer_type { High, Low, Sub };
enum process_type { Synthesis, Monitoring };
#define PI 3.1415926535897932384626433832795
enum polarization { S, P };

class design {
private:
	std::vector <double> d;
	double d_s;
	double syn_angle;
	double mon_angle;

	std::pair <double, double> n_wl_descriptor; //first: first_wl, second: wl_step

	std::vector <double> n_H;
	std::vector <double> n_L;
	std::vector <double> n_S;
public:
	design() {
		d.emplace_back(0);
		n_L.emplace_back(0);
		n_H.emplace_back(0);
		n_S.emplace_back(0);
		syn_angle = 0;
		mon_angle = 0;
		n_wl_descriptor.first = 0;
		n_wl_descriptor.second = 0;
	}

	design(std::vector <double> thick, double d_s, double syn_angle, double mon_angle, 
		std::pair <double, double> n_wl_descriptor, 
		std::vector <double> nH, 
		std::vector <double> nL, 
		std::vector <double> nS)
		: d(thick), d_s(d_s), syn_angle(syn_angle), mon_angle(mon_angle), n_wl_descriptor(n_wl_descriptor), n_H(nH), n_L(nL), n_S(nS)
	{
	}


	design(const design& copy) {
		d = copy.d;
		n_L = copy.n_L;
		n_H = copy.n_H;
		n_S = copy.n_S;
		d_s = copy.d_s;
		syn_angle = copy.syn_angle;
		mon_angle = copy.mon_angle;
		n_wl_descriptor = copy.n_wl_descriptor;
	}

	design(const design& copy, std::vector <double> d_new) {
		d = d_new;
		n_L = copy.n_L;
		n_H = copy.n_H;
		n_S = copy.n_S;
		d_s = copy.d_s;
		syn_angle = copy.syn_angle;
		mon_angle = copy.mon_angle;
		n_wl_descriptor = copy.n_wl_descriptor;
	}

	~design() {
		d.~vector();
		n_H.~vector();
		n_L.~vector();
		n_S.~vector();
	}

	double get_angle(int n_layer, double n, process_type p) {
		//�������, ������������ ���� ������� � ��������������� ����
		//��� n_layer = 0 - ���� �������, ��� n_layer = layer + 1 - ���� ������� � ��������
		if (p == Synthesis)
		{
			if (n_layer == -1)
				return syn_angle;
			else
				return asin(1 / n * sin(syn_angle));
		}
		else
		{
			if (n_layer == -1)
				return mon_angle;
			else
				return asin(1 / n * sin(mon_angle));
		}
	}

	double get_n(layer_type type, double wavelength) {
		//�������, ������������ �������� ������������ �����������
		//�� ���� � ����� ����� 
		switch (type) {
		case Low:
			return n_L[(wavelength - n_wl_descriptor.first) / n_wl_descriptor.second]; // 1.45;
		case High:
			return n_H[(wavelength - n_wl_descriptor.first) / n_wl_descriptor.second]; // 2.35;
		case Sub:
			return n_S[(wavelength - n_wl_descriptor.first) / n_wl_descriptor.second]; // 1.52;
		}
	}

	std::vector <double> get_n_vec(layer_type type) {
		switch (type) {
		case Low:
			return n_L;
		case High:
			return n_H;
		case Sub:
			return n_S;
		}
	}

	double get_d(int n) {
		return d[n];
	}

	std::vector <double> get_d_vec() {
		return d;
	}

	double get_ds() {
		return d_s;
	}

	int size() {
		return d.size();
	}

	std::pair <double, double> get_n_wl_descriptor() {
		return n_wl_descriptor;
	}
};



double T1(polarization pol, design design, double wavelength, process_type p);
double R1(polarization pol, design design, double wavelength, process_type p);
double T(polarization pol, design design, double wavelength, process_type p);

double SynthesisFunctional(design _design, std::vector<double> target_s, std::vector<double> target_p, double wl_left, double wl_right);
double MeritFunction(design _design, std::vector<double> target_s, std::vector<double> target_p, double wl_left, double wl_right);
void S_factors(std::string csv_file, std::string csv_file_rand, std::string output_file, int number_of_samples, design _design, 
			   std::vector<double> target_s, std::vector<double> target_p, double wl_left, double wl_right);

Eigen::Matrix <double, Eigen::Dynamic, 1> Gradient(polarization pol, design design, double wavelength, process_type p);
Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Hessian(polarization pol, design _design, double wavelength, process_type p);

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Synthesis_Hessian_Pol(polarization pol, design _design, std::vector <double> target, double left_wl, double right_wl);
Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Synthesis_Hessian(design _design, std::vector <double> target_s, std::vector <double> target_p, double left_wl, double right_wl);

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Control_Hessian(design _design, double left_wl, double right_wl);

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Create_W(design _design, std::vector <double> target_s, std::vector <double> target_p, 
																double left_wl, double right_wl, bool extention);
std::vector <double> read_file(std::string filename);
std::vector <double> read_n_file(std::string filename, std::pair <double, double>& n_wl_descriptor);
std::vector <double> Line2Vec(std::string line);

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Calculate_EV(design _design, std::vector <double> target_s, std::vector <double> target_p, double left_wl, double right_wl,
Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>& values, Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>& vectors);
