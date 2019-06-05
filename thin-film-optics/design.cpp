#include "design.h"
using namespace std::complex_literals;

double T1(polarization pol, design design, double wavelength, process_type p) {

	Eigen::Matrix <std::complex <double>, 2, 2> M;
	M << 1, 0, 0, 1;
	double k = (double)2 * PI / (double)wavelength;
	double n_a = 1;
	double n_s = design.get_n(Sub, wavelength);
	double n_h = design.get_n(High, wavelength);
	double n_l = design.get_n(Low, wavelength);
	double q_a, q_s;
	const int m = design.size();

	if (pol == S) {
		q_a = n_a * cos(design.get_angle(-1, n_h, p));
		q_s = n_s * cos(design.get_angle(m, n_s, p));
	}
	else {
		q_a = n_a / cos(design.get_angle(-1, n_h, p));
		q_s = n_s / cos(design.get_angle(m, n_s, p));
	}
	for (int i = m - 1; i >= 0; i--) {

		double n = ((i + 1) % 2 == 1 ? n_h : n_l);
		std::complex <double> phi, q;
		phi = k * n * cos(design.get_angle(i, n, p)) * design.get_d(i);
		//		std::cout << "cos_" << i+1 << " = " << cos(design.get_angle(i, n, p)) << std::endl;
		//		std::cout << "phi_" << i+1 << " = " << phi << std::endl;
		if (pol == S)
			q = n * cos(design.get_angle(i, n, p));
		else
			q = n / cos(design.get_angle(i, n, p));
		Eigen::Matrix <std::complex <double>, 2, 2> M_i;
		M_i(0, 0) = cos(phi);
		M_i(0, 1) = 1i * (sin(phi) / q);
		M_i(1, 0) = 1i * q * sin(phi);
		M_i(1, 1) = cos(phi);
		M = M * M_i;
		//		std::cout << "Phi: " << phi.real() << " " << std::endl;
		//		std::cout << "K: " << k << std::endl;
	}

	//	Eigen::Vector2cd vec;
	//	vec(0) = 1;
	//	vec(1) = q_s;

	//	std::cout << std::endl << M*vec << std::endl;

	std::complex <double> t = (2 * q_a) / (q_a * M(0, 0) + q_s * M(1, 1) + q_a * q_s * M(0, 1) + M(1, 0));
	//	std::cout << "t = " << t << std::endl;
	return q_s / q_a * abs(t) * abs(t);
}

double R1(polarization pol, design design, double wavelength, process_type p) {

	Eigen::Matrix <std::complex <double>, 2, 2> M;
	M << 1, 0, 0, 1;
	double k = 2 * PI / wavelength;
	double n_a = 1;
	double n_s = design.get_n(Sub, wavelength);
	double n_h = design.get_n(High, wavelength);
	double n_l = design.get_n(Low, wavelength);
	double q_a, q_s;
	const int m = design.size();

	if (pol == S) {
		q_a = n_a * cos(design.get_angle(-1, n_h, p));
		q_s = n_s * cos(design.get_angle(m, n_s, p));
	}
	else {
		q_a = n_a / cos(design.get_angle(-1, n_h, p));
		q_s = n_s / cos(design.get_angle(m, n_s, p));
	}
	for (int i = m - 1; i >= 0; i--) {

		double n = ((i + 1) % 2 == 1 ? n_h : n_l);
		std::complex <double> phi, q;
		phi = k * n * cos(design.get_angle(i, n, p)) * design.get_d(i);

		if (pol == S)
			q = n * cos(design.get_angle(i, n, p));
		else
			q = n / cos(design.get_angle(i, n, p));
		Eigen::Matrix <std::complex <double>, 2, 2> M_i;
		M_i(0, 0) = cos(phi);
		M_i(0, 1) = 1i * (sin(phi) / q);
		M_i(1, 0) = 1i * q * sin(phi);
		M_i(1, 1) = cos(phi);
		M = M * M_i;
	}
	std::complex <double> r = (q_a * M(0, 0) - q_s * M(1, 1) + q_a * q_s * M(0, 1) - M(1, 0)) / (q_a * M(0, 0) + q_s * M(1, 1) + q_a * q_s * M(0, 1) + M(1, 0));
	return abs(r) * abs(r);
}


double T(polarization pol, design design, double wavelength, process_type p) {
	double khi = 0;//6.3272/1000000000;
	double exponent = exp(-4 * PI * design.get_ds() * khi / wavelength);
	double n_a = 1, n_s = design.get_n(Sub, wavelength);
	double q_a, q_s;
	const int m = design.size();

	if (pol == S) {
		q_a = n_a * cos(design.get_angle(-1, n_a, p));
		q_s = n_s * cos(design.get_angle(m, n_s, p));
	}
	else {
		q_a = n_a / cos(design.get_angle(-1, n_a, p));
		q_s = n_s / cos(design.get_angle(m, n_s, p));
	}
	double T2 = 4 * q_a * q_s / ((q_a + q_s) * (q_a + q_s));
	double R2 = (q_a - q_s) * (q_a - q_s) / ((q_a + q_s) * (q_a + q_s));
	if (design.get_ds() == 0)
	{
		T2 = 1;
		R2 = 0;
	}
	//	std::cout << T2 << std::endl << R2 << std::endl;
	return T1(pol, design, wavelength, p) * T2 * exponent / (1 - R1(pol, design, wavelength, p) * R2 * exponent * exponent);
}

Eigen::Matrix <double, Eigen::Dynamic, 1> Gradient(polarization pol, design design, double wavelength, process_type p) {
	int m = design.size();
	double n_s = design.get_n(Sub, wavelength);
	double n_a = 1;
	double n_h = design.get_n(High, wavelength);
	double n_l = design.get_n(Low, wavelength);
	double k = 2 * PI / wavelength;
	double q_a, q_s;
	if (pol == S) {
		q_a = n_a * cos(design.get_angle(-1, n_a, p));
		q_s = n_s * cos(design.get_angle(m, n_s, p));
	}
	else {
		q_a = n_a / cos(design.get_angle(-1, n_a, p));
		q_s = n_s / cos(design.get_angle(m, n_s, p));
	}
	double R2 = (q_a - q_s) * (q_a - q_s) / ((q_a + q_s) * (q_a + q_s));
	std::vector <double> phi, q;
	std::vector <std::complex <double>> u, v;
	double n;
	for (int i = 0; i < m; i++) {
		n = ((i + 1) % 2 == 1 ? n_h : n_l);
		phi.push_back(k * n * cos(design.get_angle(i, n, p)) * design.get_d(i));
		if (pol == S)
			q.push_back(n * cos(design.get_angle(i, n, p)));
		else
			q.push_back(n / cos(design.get_angle(i, n, p)));

		if (i == 0) {
			u.push_back(cos(phi[i]) + 1i * q_s * sin(phi[i]) / q[i]);
			v.push_back(q_s * cos(phi[i]) + 1i * q[i] * sin(phi[i]));
		}
		else {
			u.push_back(u[i - 1] * cos(phi[i]) + 1i * v[i - 1] * sin(phi[i]) / q[i]);
			v.push_back(v[i - 1] * cos(phi[i]) + 1i * u[i - 1] * q[i] * sin(phi[i]));
		}
	}

	Eigen::Matrix <std::complex <double>, Eigen::Dynamic, Eigen::Dynamic> dU_dD, dV_dD;
	double g, h;
	dU_dD.resize(m, m);
	dV_dD.resize(m, m);
	for (int j = 0; j < m; j++) {
		n = ((j + 1) % 2 == 1 ? n_h : n_l);
		g = n * cos(design.get_angle(j, n, p)) / q[j];
		h = n * cos(design.get_angle(j, n, p)) * q[j];

		for (int i = 0; i < m; i++) {
			if (i < j) {
				dU_dD(i, j) = 0;
				dV_dD(i, j) = 0;
			}
			else if (i == j) {
				dU_dD(i, j) = 1i * k * v[j] * g;
				dV_dD(i, j) = 1i * k * u[j] * h;
			}
			else {
				dU_dD(i, j) = dU_dD(i - 1, j) * cos(phi[i]) + 1i / q[i] * dV_dD(i - 1, j) * sin(phi[i]);
				dV_dD(i, j) = 1i * q[i] * dU_dD(i - 1, j) * sin(phi[i]) + dV_dD(i - 1, j) * cos(phi[i]);
			}
		}
	}

	std::vector <std::complex <double>> dt_dD, dr_dD;
	for (int j = 0; j < m; j++) {
		dt_dD.push_back(-2 * q_a / ((q_a * u[m - 1] + v[m - 1]) * (q_a * u[m - 1] + v[m - 1])) * (q_a * dU_dD(m - 1, j) + dV_dD(m - 1, j)));
		dr_dD.push_back(2 * q_a / ((q_a * u[m - 1] + v[m - 1]) * (q_a * u[m - 1] + v[m - 1])) * (v[m - 1] * dU_dD(m - 1, j) - u[m - 1] * dV_dD(m - 1, j)));
	}

	std::complex <double> t = 2 * q_a / (q_a * u[m - 1] + v[m - 1]);
	std::complex <double> r = (q_a * u[m - 1] - v[m - 1]) / (q_a * u[m - 1] + v[m - 1]);


	std::vector <double> dT1_dD, dR1_dD;
	for (int j = 0; j < m; j++) {
		dT1_dD.push_back(real(q_s / q_a * (dt_dD[j] * conj(t) + conj(dt_dD[j]) * t)));
		//		std::cout << real(q_s / q_a * (dt_dD[j] * conj(t) + conj(dt_dD[j]) * t)) << " + " << imag(q_s / q_a * (dt_dD[j] * conj(t) + conj(dt_dD[j]) * t)) << " i" << "    abs : = " << abs(q_s / q_a * (dt_dD[j] * conj(t) + conj(dt_dD[j]) * t)) <<std::endl;
		dR1_dD.push_back(real(dr_dD[j] * conj(r) + conj(dr_dD[j]) * r));
		//		std::cout << real(dr_dD[j] * conj(r) + conj(dr_dD[j]) * r) << " + " << imag(dr_dD[j] * conj(r) + conj(dr_dD[j]) * r) << " i" << "    abs : = " << abs(dr_dD[j] * conj(r) + conj(dr_dD[j]) * r) << std::endl;
	}

	double khi = 0;//6.3272 /1000000000;
	double exponent = exp(-4 * PI * design.get_ds() * khi / wavelength);
	if (design.get_ds() == 0)
		R2 = 0;
	double dT_dT1 = T(pol, design, wavelength, p) / T1(pol, design, wavelength, p);
	double dT_dR1 = T(pol, design, wavelength, p) * R2 * exponent * exponent / (1 - R1(pol, design, wavelength, p) * R2 * exponent * exponent);

	Eigen::Matrix <double, Eigen::Dynamic, 1> dT_dD;
	dT_dD.resize(m, 1);
	for (int i = 0; i < m; i++)
		dT_dD[i] = (dT_dT1 * dT1_dD[i] + dT_dR1 * dR1_dD[i]);

	return dT_dD;
}

std::vector <double> read_file(std::string filename) {
	std::ifstream f;
	std::string str;
	std::vector <double> t;
	f.open(filename);
	if (!f.is_open())
		std::cout << "File \'" << filename << "\' hasn't been opened for some reason!" << std::endl;
	else {
		while (!f.eof()) {
			std::getline(f, str);
			if (str == "")
				break;
			t.push_back(std::stod(str));
		}
	}
	return t;
}

std::vector <double> read_n_file(std::string filename, std::pair <double, double>& n_wl_descriptor) {
	std::ifstream f;
	std::string str;
	std::vector <double> t;
	double init, fin, step;
	f.open(filename);
	if (!f.is_open())
		std::cout << "File \'" << filename << "\' hasn't been opened for some reason!" << std::endl;
	else {
		std::getline(f, str);
		init = std::stod(str);
		std::getline(f, str);
		fin = std::stod(str);
		std::getline(f, str);
		step = std::stod(str);

		if (init != n_wl_descriptor.first || step != n_wl_descriptor.second)
			std::cout << "Wave length range of file does NOT correspond with the asked range!" << std::endl;
		while (!f.eof()) {
			std::getline(f, str);
			if (str == "")
				break;
			t.push_back(std::stod(str));
		}
	}
	return t;
}

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Hessian(polarization pol, design _design, double wavelength, process_type p) {
	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian;
	int m = _design.size();
	hessian.resize(m, m);

	double delta = 0.01;
	double n_s = _design.get_n(Sub, wavelength);

	std::vector <double> d = _design.get_d_vec();
	Eigen::Matrix <double, Eigen::Dynamic, 1> gradip, gradim, gradjp, gradjm;
	gradip.resize(m, 1);
	gradim.resize(m, 1);
	gradjp.resize(m, 1);
	gradjm.resize(m, 1);

	for (int i = 0; i < m; i++) {
		std::vector <double> dip = d, dim = d;
		dip[i] += delta;
		dim[i] -= delta;

		design desip(dip, _design.get_ds(), _design.get_angle(-1, n_s, Synthesis), _design.get_angle(-1, n_s, Monitoring), _design.get_n_wl_descriptor(), 
					 _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub)),
			   desim(dim, _design.get_ds(), _design.get_angle(-1, n_s, Synthesis), _design.get_angle(-1, n_s, Monitoring), _design.get_n_wl_descriptor(), 
					 _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub));
		gradip = Gradient(pol, desip, wavelength, p);
		gradim = Gradient(pol, desim, wavelength, p);


		for (int j = 0; j < m; j++) {
			std::vector <double> djp = d, djm = d;
			djp[j] += delta;
			djm[j] -= delta;

			design desjp(djp, _design.get_ds(), _design.get_angle(-1, n_s, Synthesis), _design.get_angle(-1, n_s, Monitoring), _design.get_n_wl_descriptor(), 
						 _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub)),
				   desjm(djm, _design.get_ds(), _design.get_angle(-1, n_s, Synthesis), _design.get_angle(-1, n_s, Monitoring), _design.get_n_wl_descriptor(), 
						 _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub));

			gradjp = Gradient(pol, desjp, wavelength, p);
			gradjm = Gradient(pol, desjm, wavelength, p);

			hessian(i, j) = (gradjp[i] - gradjm[i] + gradip[j] - gradim[j]) / (4 * delta);
		}
	}
	return hessian;
}

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Synthesis_Hessian_Pol(polarization pol, design _design, std::vector <double> target, double left_wl, double right_wl) {
	printf("SYNTHESIS HESSIAN DOES NOT WORK PROPERLY!!\n");
	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian, d2T_dD2;
	Eigen::Matrix <double, Eigen::Dynamic, 1> grad;
	int m = _design.size();
	double percent = 0;
	hessian.resize(m, m);
	d2T_dD2.resize(m, m);
	grad.resize(m, 1);
	hessian.setZero();
	for (int wl = left_wl; wl <= right_wl; wl++) {
		if ((wl >= 534 && wl <= 541) || (wl >= 567 && wl <= 574) || (wl >= 578 && wl <= 588) || (wl >= 630 && wl <= 638) ||
			(wl >= 641 && wl <= 648) || (wl >= 681 && wl <= 688) || (wl >= 692 && wl <= 701)) // Exlusion from Target mesh
			continue;
		grad = Gradient(pol, _design, wl, Synthesis);
		hessian += 2 * grad * grad.transpose();

		d2T_dD2 = Hessian(pol, _design, wl, Synthesis);

		hessian += 2 * (T(pol, _design, wl, Synthesis) - target[wl - left_wl]) * d2T_dD2;

		percent = 100 * (double)(wl - left_wl) / (right_wl - left_wl);
		system("cls");
		char c = (pol == S) ? 'S' : 'P';
		printf("Calculating the synthesis hessian for %c polarization: %.2f%%\n", c, percent);
	}
	return hessian;
}

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Synthesis_Hessian(design _design, std::vector <double> target_s, std::vector <double> target_p, 
																		 double left_wl, double right_wl) {
	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian_s, hessian_p;
	hessian_s = Synthesis_Hessian_Pol(S, _design, target_s, left_wl, right_wl);
	hessian_p = Synthesis_Hessian_Pol(P, _design, target_p, left_wl, right_wl);
	return hessian_s + hessian_p;
}


Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Control_Hessian(design _design, double left_wl, double right_wl) {
	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian;
	Eigen::Matrix <double, Eigen::Dynamic, 1> grad_s, grad_p;

	int m = _design.size();
	double percent = 0;
	hessian.resize(m, m);
	grad_s.resize(m, 1);
	grad_p.resize(m, 1);

	hessian.setZero();
	for (int wl = left_wl; wl <= right_wl; wl++) {
		grad_s = Gradient(S, _design, wl, Monitoring);
		grad_p = Gradient(P, _design, wl, Monitoring);

		hessian += grad_s * grad_s.transpose();
		hessian += grad_p * grad_p.transpose();

		percent = 100 * (double)(wl - left_wl) / (right_wl - left_wl);
		//		printf("Calculating the control hessian : %.2f%%\n", percent);
	}
	return hessian;
}

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Create_W(design _design, std::vector <double> target_s, std::vector <double> target_p, 
																double left_wl, double right_wl, bool extention) {
	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> W;
	int m = _design.size();

	//	W = Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>::Zero(m*(m + 1) / 2, m);
	if (extention)
		W.resize(m * (m + 1) / 2 + m, m);
	else
		W.resize(m * (m + 1) / 2, m);
	//	W.setZero();

	std::vector <double> thicks;
	for (int i = 0; i < m; i++) {
		thicks.push_back(_design.get_d(i));
		design temp(thicks, _design.get_ds(), _design.get_angle(-1, 1, Synthesis), _design.get_angle(-1, 1, Monitoring), 
					_design.get_n_wl_descriptor(), _design.get_n_vec(High), _design.get_n_vec(Low), _design.get_n_vec(Sub));
		Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian = Control_Hessian(temp, left_wl, right_wl);
		//		std::cout << hessian << std::endl << std::endl;

		Eigen::SelfAdjointEigenSolver <Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>> es(i + 1);
		es.compute(hessian);
		Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> eigvals = es.eigenvalues();
		Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> eigvecs = es.eigenvectors();

		eigvecs.transposeInPlace();

		Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> blockW;
		blockW.resizeLike(eigvecs);

		for (int j = 0; j <= i; j++) {
			blockW.row(j) = std::sqrt(abs(eigvals(j))) * eigvecs.row(j);
		}

		//		std::cout << es.eigenvectors().transpose() << std::endl << std::endl;
		//		std::cout << es.eigenvalues().transpose() << std::endl << std::endl;
		//		std::cout << blockW << std::endl << std::endl;

		W.block(i * (i + 1) / 2, 0, i + 1, i + 1) = blockW;
		W.block(i * (i + 1) / 2, i + 1, i + 1, m - (i + 1)).setZero();
	}
	if (extention) {
		Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian = Synthesis_Hessian(_design, target_s, target_p, 1528, 1566);

		Eigen::SelfAdjointEigenSolver <Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>> es(m);
		es.compute(hessian);
		Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> eigvals = es.eigenvalues();
		Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> eigvecs = es.eigenvectors();

		eigvecs.transposeInPlace();

		Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> blockW;
		blockW.resizeLike(eigvecs);

		for (int j = 0; j < m; j++) {
			blockW.row(j) = std::sqrt(abs(eigvals(j))) * eigvecs.row(j);
		}
		W.block((m + 1) * m / 2, 0, m, m) = blockW;
	}
	return W;
}

Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> Calculate_EV(design _design, std::vector <double> target_s, std::vector <double> target_p, double left_wl, double right_wl,
																	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>& values, 
																	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>& vectors) {
	int m = _design.size();

	Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian = Synthesis_Hessian(_design, target_s, target_p, left_wl, right_wl);
	Eigen::SelfAdjointEigenSolver <Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>> es(m);
	es.compute(hessian);

	values = es.eigenvalues();
	vectors = es.eigenvectors();
	return hessian;
}

double MeritFunction(design _design, std::vector<double> target_s, std::vector<double> target_p, double wl_left, double wl_right) {
	double MF = 0;
	double wl_grid = 0;
	std::pair <double, double> n_wl_descriptor = _design.get_n_wl_descriptor();
	double step = n_wl_descriptor.second;
	for (double wl = wl_left; wl <= wl_right; wl += step) {
		//if ((wl >= 534 && wl <= 541) || (wl >= 567 && wl <= 574) || (wl >= 578 && wl <= 588) || (wl >= 630 && wl <= 638) ||
		//	(wl >= 641 && wl <= 648) || (wl >= 681 && wl <= 688) || (wl >= 692 && wl <= 701)) // Exlusion from Target mesh
		//	continue;
		if (target_s[(wl - wl_left) / step] < 0)
			continue;
		MF += pow(T(S, _design, wl, Synthesis) - target_s[(wl - wl_left) / step], 2);
		MF += pow(T(P, _design, wl, Synthesis) - target_p[(wl - wl_left) / step], 2);
		wl_grid += 1;
	}
	return sqrt(MF / wl_grid);
}

double SynthesisFunctional(design _design, std::vector<double> target_s, std::vector<double> target_p, double wl_left, double wl_right) {
	double F = 0;
	std::pair <double, double> n_wl_descriptor = _design.get_n_wl_descriptor();
	double step = n_wl_descriptor.second;
	for (double wl = wl_left; wl <= wl_right; wl += step) {
		//if ((wl >= 534 && wl <= 541) || (wl >= 567 && wl <= 574) || (wl >= 578 && wl <= 588) || (wl >= 630 && wl <= 638) ||
		//	(wl >= 641 && wl <= 648) || (wl >= 681 && wl <= 688) || (wl >= 692 && wl <= 701)) // Exlusion from Target mesh
		//	continue;
		if (target_s[(wl - wl_left) / step] < 0)
			continue;
		//printf("[S] wl = %1.1f:\tT = %1.7f,\tTar = %1.1f\tdT = %1.12f\n", wl, T(S, _design, wl, Synthesis), target_s[(wl - wl_left) / step], 
		//	pow(T(S, _design, wl, Synthesis) - target_s[(wl - wl_left) / step], 2)*10000);
		//printf("[P] wl = %1.1f:\tT = %1.7f,\tTar = %1.1f\tdT = %1.12f\n", wl, T(P, _design, wl, Synthesis), target_p[(wl - wl_left) / step],
		//	pow(T(P, _design, wl, Synthesis) - target_p[(wl - wl_left) / step], 2) * 10000);

		F += pow(T(S, _design, wl, Synthesis) - target_s[(wl - wl_left) / step], 2);
		F += pow(T(P, _design, wl, Synthesis) - target_p[(wl - wl_left) / step], 2);
	}
	return F;
}


std::vector<double> Line2Vec(std::string line) {
	std::vector <double> errors_vec;
	std::string number;
	for (int i = 0; i < line.length(); i++) {
		if (line[i] == ',') {
			errors_vec.push_back(std::stod(number));
			number = "";
		}
		else
			number += line[i];
	}
	errors_vec.push_back(std::stod(number));
	//	std::cout << '\n';
	return errors_vec;
}

std::vector<double> normalize(std::vector<double> thickness) {
	std::vector<double> result;
	double sums = 0;
	for (int i = 0; i < thickness.size(); i++)
		sums += thickness[i] * thickness[i];
	for (int i = 0; i < thickness.size(); i++)
		result.push_back(thickness[i] / sqrt(sums));
	//printf("%f\n", sqrt(sums));
	return result;
}

std::vector<double> add(std::vector<double> one, std::vector<double> two) {
	std::vector<double> result;
	//	printf("%d, %d", one.size(), two.size());
	assert(one.size() == two.size());
	for (int i = 0; i < one.size(); i++)
		result.push_back(one[i] + two[i]);
	return result;
}


void S_factors(std::string csv_file, std::string csv_file_rand, std::string output_file, int number_of_samples, 
	design _design, std::vector<double> target_s, std::vector<double> target_p, double wl_left, double wl_right) {

	std::ifstream input_file_stream;
	std::ofstream output_file_stream;
	output_file_stream.open(output_file);
	if (!output_file_stream.is_open())
		std::cout << "File " << output_file << " hasn't been opened for some reason...\n";
	//calculate dF for vectors;
	Eigen::VectorXd dF;
	dF.resize(number_of_samples);

	std::string error_line;
	std::vector <double> errors_vec;
	std::vector <double> thickness = _design.get_d_vec();;
	std::vector <double> thick_with_err_norm;

	int rand_attempts = 10000;
	double dF_rand_mean = 0;
	double dF_rand = 0;
	input_file_stream.open(csv_file_rand);
	double F = SynthesisFunctional(_design, target_s, target_p, wl_left, wl_right);
	for (int i = 0; i < rand_attempts; i++)
	{
		//printf("%d\n", i);
		std::getline(input_file_stream, error_line);
		errors_vec = Line2Vec(error_line);

		thick_with_err_norm = add(thickness, normalize(errors_vec));
		design err_design(_design, thick_with_err_norm);

		double F_err = SynthesisFunctional(err_design, target_s, target_p, wl_left, wl_right);
		dF_rand = F_err - F;
		if (dF_rand < 0)
			printf("Negative dF (rand)!!!\n");
		dF_rand_mean += dF_rand;

		if (i % 1000 == 0)
			printf("Done (Rand): %.3f%%\n", 100 * double(i) / rand_attempts);
	}
	dF_rand_mean /= rand_attempts;
	input_file_stream.close();
	printf("Mean dF(rand): %f\n", dF_rand_mean);
	Eigen::VectorXd S;
	S.resize(number_of_samples);
	input_file_stream.open(csv_file); 
	if (!input_file_stream.is_open())
		std::cout << "File " << csv_file << " hasn't been opened for some reason...\n";
	for (int i = 0; i < number_of_samples; i++)
	{
		std::getline(input_file_stream, error_line);
		errors_vec = Line2Vec(error_line);
		//printf("\n---%d---\n", i);
		thick_with_err_norm = add(thickness, normalize(errors_vec));
		design err_design(_design, thick_with_err_norm);

		double F_err = SynthesisFunctional(err_design, target_s, target_p, wl_left, wl_right);
		//printf("%0.10f\t%0.10f\n", F_err, F);
		dF[i] = F_err - F;
		//printf("%0.10f\n", dF[i]);
		S[i] = dF_rand_mean / dF[i];

		if (i % 1000 == 0)
 			printf("Done (Data): %.3f%%\n", 100 * double(i) / number_of_samples);
	}
	input_file_stream.close();
	if (!output_file_stream.is_open())
		printf("Don't know where to write!!!\n");
	output_file_stream << S;
	output_file_stream.close();
}