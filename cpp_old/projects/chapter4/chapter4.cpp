// Chapter 4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <valarray>
#include <cmath>

using namespace std;

struct Range
{
	double m_t_step;
	double bethe_range;
};

double stop_power(int Z, double A, double E, double J);
double range(int Z, double A, double E0, double J, double density);
Range range_HD(int Z, double A, double E0, double J, double density);


std::random_device random;
std::mt19937 generator(random());
std::uniform_real_distribution<double> distribution(0, 1);


int main()
{
	int Z, repeat;
	double A, density, E0;

	cout << "Chapter 4: the Plural Scattering Model " << endl;
	cout << "please enter Z(atomic number), A (atomic weight,g/mol), density(g/cm3), and E0 (keV) " << endl;
	cin >> Z >> A >> density >> E0;
	cout << "please enter the number of electrons" << endl;
	cin >> repeat;

	ofstream ofile;
	ofile.open("myfile_input.txt");
	ofile << "Z=" << Z << "\nA=" << A << " \ndensity=" << density << "\nE0=" << E0 << "\n";

	ofstream ofile1;
	ofile1.open("myfile_output.csv");

	double J = (9.76*Z + 58.5 / pow(Z, 0.19)) / 1000;  // keV

	double distance_step = range(Z, A, E0, J, density);
	double m_t_step = distance_step;
	//Range range_bethe = range_HD(Z, A, E0, J, density);
	//double distance_step = range_bethe.bethe_range / 50;
	//double m_t_step = range_bethe.m_t_step;

	//cout << distance_step << endl;
	double p0 = 0.394*pow(Z, 0.4) / E0;

	double x0 = 0, y0 = 0, z0 = 0;
	double cx = 0, cy = 0, cz = -1;

	int backscatter = 0;

	ofile1 << "m" << "," << "n" << "," << "xn" << "," << "yn" << "," << "zn" << "," << "E" << "," << "cx" << "," << "cy" << "," << "cz" << "," << "cosphi" << "," << "sinphi" << "," << "angle" << endl;
	for (int m = 1; m <= repeat; ++m)
	{
		double xn = x0, yn = y0, zn = z0;
		double ca = cx, cb = cy, cc = cz;
		double E = E0;

		double cosphi = 1;
		double sinphi = 0;

		double RND3 = distribution(generator);
		double angle = 2 * 3.1415926*RND3;

		for (int n = 1; n <= 55 & E>0.05; ++n)
		{
			ofile1 << m << "," << n << "," << xn << "," << yn << "," << zn << "," << E << "," << cx << "," << cy << "," << cz << "," << cosphi << "," << sinphi << "," << angle << endl;
			//cout <<m<<"  " <<n << " " << xn << " " << yn << " " << zn << " " << E << "\n";

			double AN = -(ca / cc);
			double AM = 1 / sqrt(1 + AN*AN);
			double V1 = AM*sinphi;
			double V2 = AM*AN*sinphi;
			double V3 = cos(angle);
			double V4 = sin(angle);

			double ca_n = ca*cosphi + V1*V3 + cb*V2*V4;
			double cb_n = cb*cosphi + (V4*(cc*V1 - ca*V2));
			double cc_n = cc*cosphi + V2*V3 - cb*V1*V4;

			ca = ca_n;
			cb = cb_n;
			cc = cc_n;

			xn = xn + distance_step*ca;
			yn = yn + distance_step*cb;
			zn = zn + distance_step*cc;          // cm

			const double RND2 = distribution(generator);
			const double ta = 0.0144*((1 / sqrt(RND2)) - 1)*Z / E / 2 / p0;

			cosphi = (1 - (ta* ta)) / (1 + (ta*ta));
			sinphi = (ta + ta) / (1 + ta*ta);

			RND3 = distribution(generator);
			angle = 2 * 3.1415926*RND3;

			double dES = stop_power(Z, A, E, J);
			double A1 = dES*m_t_step;
			double A2 = m_t_step*stop_power(Z, A, E - A1 / 2, J);
			double A3 = m_t_step*stop_power(Z, A, E - A2 / 2, J);
			double A4 = m_t_step*stop_power(Z, A, E - A3, J);
			E = E - (A1 + 2 * A2 + 2 * A3 + A4) / 6;

			if (zn >= 0)
			{
				backscatter = backscatter + 1;
				break;
			}

		}
	}

	double backscatteredcoefficient = backscatter / static_cast<double>(repeat);
	std::cout << "backscattered coefficient=" << backscatteredcoefficient << "\n";
	ofile << "backscattered coefficient=" << backscatteredcoefficient << endl;

	ofile.close();

	return 0;
}


double stop_power(int Z, double A, double E, double J)
{
	double dES = 78500 * Z / A / E*log(1.166*(E + 0.85*J) / J);      //kev/(g/cm2)
	return dES;
}

double range(int Z, double A, double E0, double J, double density)
{
	double bethe_range = 0;
	double E_step = (E0 - 0.05) / 200;
	double E = E0;
	for (int n = 1; n <= 200; n++)
	{
		double dES = stop_power(Z, A, E, J);
		bethe_range = bethe_range + 1 / dES*E_step / density;
		E = E - E_step;
	}
	bethe_range = bethe_range / 50;       //step:50
	return bethe_range;
}

Range range_HD(int Z, double A, double E0, double J, double density)
{
	const int number_energy_step = 20;
	double fs = 0.0;
	double E_step = (E0 - 0.05) / 200;
	double E = E0;

	// A Simpson's rule integration.
	for (int m = 1; m <= number_energy_step + 1; m++)
	{
		const double energy_keV = (m - 1) * E0 / number_energy_step;
		const double dES = stop_power(Z, A, energy_keV, J);
		const double f = 1.0 / dES;

		double l = 2.0;
		if (m % 2 == 0)
		{
			l = 4.0;
		}
		if (m == 1 || m == number_energy_step + 1)
		{
			l = 1.0;
		}

		fs += l*f;
	}
	

	double bethe_range = fs * E0 / 60.0;
	double m_t_step = bethe_range / 50.0;

	bethe_range *= 10000.0 / density;

	Range range;
	range.m_t_step = m_t_step;
	range.bethe_range = bethe_range;
	return range;
}
