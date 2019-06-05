#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include "design.h"
#include <complex>
#include <cmath>

using namespace std::complex_literals;

Eigen::Matrix <double, Eigen::Dynamic, 1> Gradient_num(polarization pol, design _design, double wavelength, process_type p);

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Hessian_num(design _design, std::vector <double> target_s, std::vector <double> target_p, double wl_left, double wl_right);


int main() {
	std::string design_name = "jl";
	std::string design_name_cap = "JL";
	std::string mon_range = "(450-1000)";
	std::string syn_range = "(531-765)";

	std::vector <double> thick, thick2, target, errors, thick_er, nH, nL, nS;
	std::string file_design = "Design Data/" + design_name + "_v2.csv";
	std::string file_target = "Design Data/" + design_name + "_target.csv";

	std::string file_errors_csv = "Design Data/Error_Vecs/" + design_name_cap + "_errors_V2.csv";

	std::string file_rand_errors_csv = "Design Data/Error_Vecs/Rand_errors.csv";

	std::string file_nH = "Design Data/Refractive_indices/" + design_name_cap + "/nH_Ta2O5.csv";
	std::string file_nL = "Design Data/Refractive_indices/" + design_name_cap + "/nL_SiO2.csv";
	std::string file_nS = "Design Data/Refractive_indices/" + design_name_cap + "/nS_FS.csv";

	std::string file_output_sv = "Output Data/" + design_name_cap + "/" + mon_range + "_sv.txt";
	std::string file_output_w  = "Output Data/" + design_name_cap + "/" + mon_range + "_W.txt";
	std::string file_output_V  = "Output Data/" + design_name_cap + "/" + mon_range + "_svd_V.txt";

	std::string file_output_ev = "Output Data/" + design_name_cap + "/" + mon_range + "_eVals.txt";
	std::string file_output_eV = "Output Data/" + design_name_cap + "/" + mon_range + "_eVecs.txt";
	std::string file_output_SH = "Output Data/" + design_name_cap + "/" + syn_range + "_SH.txt";
	std::string file_output_SH_num = "Output Data/" + design_name_cap + "/" + syn_range + "_SH_num.txt";


	//std::string file_output_sv_ext = "Output Data/" + design_name_cap + "_sv_ext.txt";
	//std::string file_output_w_ext = "Output Data/" + design_name_cap + "_W_ext.txt";

	std::string file_output_S_factors = "Output Data/" + design_name_cap + "/" + syn_range + "_S_V2.txt";

	thick = read_file(file_design);
	int m = thick.size();
	target = read_file(file_target);

	std::pair <double, double> n_wl_descriptor;
	n_wl_descriptor.first = 450;
	n_wl_descriptor.second = 0.5;

	nH = read_n_file(file_nH, n_wl_descriptor);
	nL = read_n_file(file_nL, n_wl_descriptor);
	nS = read_n_file(file_nS, n_wl_descriptor);

	// ANGLE SHOULD BE ENTERED IN RADIANS!!!!

	design jl(thick, 0, 0, 0, n_wl_descriptor, nH, nL, nS);

	//design gff(thick, 1000000, 0, 0, n_wl_descriptor, nH, nL, nS);
	/*
	errors = read_file("Design Data/gff_errors.csv");
	int s = thick.size();
	for (int i = 0; i < s; i++)
		thick_er.push_back(thick[i] + errors[i]);
	for (int i = 0; i<s; i++)
		std::cout << thick[i] << "-----" << thick_er[i] << std::endl;

	design gff_73_er(thick_er, 1000000, 0, 0, n_range, nH, nL, nS);
	double MF_d = MeritFunction(gff_73_er, target, target, 1528, 1566);
	double MF_0 = MeritFunction(gff_73, target, target, 1528, 1566);
	std::cout << "MF_d = " << MF_d << std::endl << "MF_0 = " << MF_0 << std::endl;

	double norm2 = 0;
	for (int i = 0; i < s; i++)
		norm2 += errors[i] * errors[i];

	std::cout << "||d||^2 = " << norm2 << std::endl;
	double dF = (MF_d * MF_d - MF_0 * MF_0) / norm2;

	std::cout << "dF = " << dF << std::endl;
	*/
	std::ofstream f;

	//S_factors(file_errors_csv, file_rand_errors_csv, file_output_S_factors, 1000000, jl, target, target, 531, 765);
	//------------------------------------------------------------------------------------------------
/*
	Eigen::VectorXd Ts, Tp, grad;
	Ts.resize(201);
	Tp.resize(201);
	grad.resize(50);

	for (int i = 900; i <= 1100; i++)
	{
		Ts(i - 900) = R1(S, main_design, i, Synthesis);
		Tp(i - 900) = R1(P, main_design, i, Synthesis);
	}
	f.open(file_output_T_s);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << Ts;
	f.close();
	f.open(file_output_T_p);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << Tp;
	f.close();
*/
//double t = T(P, main_design, 400, Synthesis);
//std::cout << "t_mon = " << t << std::endl;

	//f.open("Output Data/jl_T_array.txt");
	//for (int w = 531; w <= 765; w++)
	//{
	//	f << w << "\t" << T(S, jl, w, Monitoring) << std::endl;
	//}
	//double MF = MeritFunction(jl, target, target, 531, 765);

	//double delta = 0.01;
	//std::vector <double> thick_p = thick, thick_m = thick;
	//thick_p[0] += delta;
	//thick_m[0] -= delta;

	//design jlp(thick_p, 0, 0, 0, n_range, nH, nL, nS);
	//design jlm(thick_m, 0, 0, 0, n_range, nH, nL, nS);
	//double F0 = pow(MeritFunction(jl, target, target, 531, 765), 2), Fp = pow(MeritFunction(jlp, target, target, 531, 765), 2), Fm = pow(MeritFunction(jlm, target, target, 531, 765), 2);
	//std::cout << "d2F(1,1): " << (Fp - 2 * F0 + Fm) / (delta*delta) << std::endl;

	//std::cout << "MF:\t" << MF << std::endl;

	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> w = Create_W(jl, target, target, 450, 1000, false);
	std::cout << "W has been computed.." << std::endl;

	Eigen::BDCSVD <Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>> svd(w, Eigen::ComputeFullU | Eigen::ComputeFullV);
	//	std::cout << svd.singularValues() << std::endl;

	//Output the W
	f.open(file_output_w);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << w;
	f.close();
	std::cout << "W has been written in the file" << std::endl;

	//output the W singular values
	f.open(file_output_sv);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << svd.singularValues();
	f.close();
	std::cout << "W's SV have been written in the file" << std::endl;

	//output the V matrix
	f.open(file_output_V);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << svd.matrixV();
	f.close();
	std::cout << "W's V have been written in the file" << std::endl;
	/*

	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> evalues, evectors, hessian;
	evalues.resize(m, 1);
	evectors.resize(m, m);
	hessian.resize(m, m);
	hessian = Calculate_EV(jl, target, target, 531, 765, evalues, evectors);

	//output the Q eigen values
	f.open(file_output_ev);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << evalues;
	f.close();
	std::cout << "Q's EigenValues have been written in the file" << std::endl;

	//output the Q eigen vectors
	f.open(file_output_eV);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << evectors;
	f.close();
	std::cout << "Q's EigenVectors have been written in the file" << std::endl;

	//output the Synthesis Hessian
	f.open(file_output_SH);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << hessian;
	f.close();
	std::cout << "Synthesis hessian have been written in the file" << std::endl;

	/*
	//output the Synthesis Hessian numerically
	f.open(file_output_SH_num);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << Hessian_num(gff_73, target, target, 1528, 1566);
	f.close();
	std::cout << "Synthesis hessian_num have been written in the file" << std::endl;

	*/
	//Calculate extended W and its Singular Values

/*	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> w_ext = Create_W(gff_73, target, target, 1100, 1600, 1);
	std::cout << "Extended W has been computed.." << std::endl;

	Eigen::JacobiSVD <Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>> svd_ext(w_ext);
	//	std::cout << svd_ext.singularValues() << std::endl;

	//Output the extended W
	f.open(file_output_w_ext);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << w_ext;
	f.close();
	std::cout << "Extended W has been written in the file" << std::endl;

	//output the extended W singular values
	f.open(file_output_sv_ext);
	if (!f.is_open())
		std::cout << "File hasn't been opened for some reason!" << std::endl;
	else
		f << svd_ext.singularValues();
	f.close();
	std::cout << "Extended W SV have been written in the file" << std::endl;
	*/
	std::cout << "Done..." << std::endl;

	system("pause");
}

Eigen::Matrix <double, Eigen::Dynamic, 1> Gradient_num(polarization pol, design _design, double wavelength, process_type p) {
	double delta = 15;
	int m = _design.size();
	double n_s = _design.get_n(Sub, wavelength);
	Eigen::Matrix <double, Eigen::Dynamic, 1> grad;
	grad.resize(m);
	std::vector <double> d = _design.get_d_vec();
	for (int i = 0; i < m; i++) {

		std::vector <double> dip = d, dim = d;
		dip[i] += delta;
		dim[i] -= delta;

		design desip(dip, _design.get_ds(), _design.get_angle(-1, n_s, Synthesis), _design.get_angle(-1, n_s, Monitoring), _design.get_n_wl_descriptor(), _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub)),
			desim(dim, _design.get_ds(), _design.get_angle(-1, n_s, Synthesis), _design.get_angle(-1, n_s, Monitoring), _design.get_n_wl_descriptor(), _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub));


		grad(i) = (T(pol, desip, wavelength, p) - T(pol, desim, wavelength, p)) / (2 * delta);
	}
	return grad;
}

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Hessian_num(design _design, std::vector <double> target_s, std::vector <double> target_p, double wl_left, double wl_right) {
	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian;
	int m = _design.size();
	hessian.resize(m, m);

	double delta = 20;

	//double n_s = _design.get_n(Sub, wl);
	std::vector <double> d = _design.get_d_vec();

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			if (i == j) {
				std::vector <double> dp = d, dm = d;
				dp[i] += delta;
				dm[i] -= delta;

				design desp(dp, _design.get_ds(), _design.get_angle(-1, 1, Synthesis), _design.get_angle(-1, 1, Monitoring), _design.get_n_wl_descriptor(), _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub));
				design desm(dm, _design.get_ds(), _design.get_angle(-1, 1, Synthesis), _design.get_angle(-1, 1, Monitoring), _design.get_n_wl_descriptor(), _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub));

				long double fp = MeritFunction(desp, target_s, target_p, wl_left, wl_right);
				long double fm = MeritFunction(desm, target_s, target_p, wl_left, wl_right);
				long double f0 = MeritFunction(_design, target_s, target_p, wl_left, wl_right);

				hessian(i, i) = (fp * fp - 2 * f0 * f0 + fm * fm) / (delta * delta);
			}
			else {
				std::vector <double> dpp = d, dpm = d, dmp = d, dmm = d;
				dpp[i] += delta;
				dpp[j] += delta;

				dpm[i] += delta;
				dpm[j] -= delta;

				dmp[i] -= delta;
				dmp[j] += delta;

				dmm[i] -= delta;
				dmm[j] -= delta;

				design despp(dpp, _design.get_ds(), _design.get_angle(-1, 1, Synthesis), _design.get_angle(-1, 1, Monitoring), _design.get_n_wl_descriptor(), _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub)),
					desmm(dmm, _design.get_ds(), _design.get_angle(-1, 1, Synthesis), _design.get_angle(-1, 1, Monitoring), _design.get_n_wl_descriptor(), _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub));
				design desmp(dmp, _design.get_ds(), _design.get_angle(-1, 1, Synthesis), _design.get_angle(-1, 1, Monitoring), _design.get_n_wl_descriptor(), _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub)),
					despm(dpm, _design.get_ds(), _design.get_angle(-1, 1, Synthesis), _design.get_angle(-1, 1, Monitoring), _design.get_n_wl_descriptor(), _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub));

				long double fpp = MeritFunction(despp, target_s, target_p, wl_left, wl_right);
				long double fpm = MeritFunction(despm, target_s, target_p, wl_left, wl_right);
				long double fmp = MeritFunction(desmp, target_s, target_p, wl_left, wl_right);
				long double fmm = MeritFunction(desmm, target_s, target_p, wl_left, wl_right);

				hessian(i, j) = (fpp * fpp - fmp * fmp + fmm * fmm - fpm * fpm) / (4 * delta * delta);
			}
		}
	}
	return hessian;
}