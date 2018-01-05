// Chapter3.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <valarray>
#include <math.h>

using namespace std;

std::random_device random;
std::mt19937 generator(random());
std::uniform_real_distribution<double> distribution(0, 1);

double compute_alpha(double atomic_number, double energy_keV);

int main()
{
	int Z,step;
	double A, density, E0, E, cx, cy,cz, x0, y0 , z0 ,xn,yn,zn,ca,cb,cc,crosssection,alpha,meanfreepath,s, cosphi,RND1,RND2,RND3,angle,AN,AM,V1,V2,V3,V4,dES,J,distance,repeat;

	cout << "Chapter 3: Single Scattering Model " << endl;
	cout << "please enter Z(atomic number), A (atomic weight,g/mol), density(g/cm3), and E0 (keV) " << endl;
	cin >> Z >> A >> density >> E0;
	cout << "please enter step and the number of electrons" << endl;
	cin >> step>>repeat;

	ofstream ofile;
	ofile.open("myfile_input.txt");
	ofile << "Z=" << Z << "\nA=" << A << " \ndensity=" << density << "\nE0=" << E0<<"\n";

	ofstream ofile1;
	ofile1.open("myfile_output.csv");

	x0 = 0;
	y0 = 0;
	z0 = 0;
	double cx0 = 0;
	double cy0 = 0;
	double cz0 = -1;              //tilt=0

	int backscatter = 0;

	for (int m = 1;m <= repeat;++m)
	{
		xn = x0;
		yn = y0;
		zn = z0;
		cx = cx0;
		cy = cy0;
		cz = cz0;
		E = E0;
		dES = 1;

		cosphi = 1;

        //alpha = 3.4 / 1000 * pow(Z, 0.67) / E;           //screening factor
        //crosssection = 5.21*pow(10, -21)*Z*Z / E / E * 4 * 3.1415926 / (alpha*(1 + alpha))*pow((E + 511) / (E + 1024), 2);        //cm2/atom
        //meanfreepath = A / (6.02*pow(10, 23)*density*crosssection);

        //RND1 = distribution(generator);
        //s = -meanfreepath*log(RND1);               //cm

        //xn = xn + s*cx;
        //yn = yn + s*cy;
        //zn = zn + s*cz;          // cm

        //J = (9.76*Z + 58.5 / pow(Z, 0.19)) / 1000;                // keV
        //dES = 78500 * Z / A / E*log(1.166*(E + 0.85*J) / J);     // keV/cm
        //E = E - s*density*dES;

        for (int n = 1; E >= 0.05 && n <= step && dES > 0; ++n) {

			//ofile1 <<m<<","<< n << "," << xn << "," << yn << "," << zn << "," << E << endl;
			//cout <<m<<"  " <<n << " " << xn << " " << yn << " " << zn << " " << E << "\n";

			alpha = compute_alpha(Z, E);           //screening factor
			crosssection = 5.21*pow(10, -21)*Z*Z / E / E * 4 * 3.1415926 / (alpha*(1 + alpha))*pow((E + 511) / (E + 1024), 2);        //cm2/atom
			meanfreepath = A / (6.02*pow(10, 23)*density*crosssection);

			RND1 = distribution(generator);
			s = -meanfreepath*log(RND1);               //cm

			RND3 = distribution(generator);
			angle = 2 * 3.1415926*RND3;

			alpha = compute_alpha(Z, E);           //screening factor
			RND2 = distribution(generator);
            cosphi = 1 - 2 * alpha*RND2 / (1 + alpha - RND2);

			AN = -(cx / cz);
			AM = 1 / sqrt(1 + AN*AN);

			V1 = AM*sqrt(1 - cosphi*cosphi);
			V2 = AM*AN*sqrt(1 - cosphi*cosphi);
			V3 = cos(angle);
			V4 = sin(angle);
			ca = cx*cosphi + V1*V3 + cy*V2*V4;
			cb = cy*cosphi + (V4*(cz*V1 - cx*V2));
			cc = cz*cosphi + V2*V3 - cy*V1*V4;

			xn = xn + s*ca;
			yn = yn + s*cb;
			zn = zn + s*cc;          // cm

			cx = ca;
			cy = cb;
			cz = cc;

			J = (9.76*Z + 58.5 / pow(Z, 0.19)) / 1000;                // keV
			dES = 78500 * Z / A / E*log(1.166*(E + 0.85*J) / J);     // keV/cm
			E = E - s*density*dES;

			if (zn >= 0)
			{
				backscatter = backscatter + 1;
				break;
			}

		}
		distance = sqrt((xn - x0)*(xn - x0) + (yn - y0)*(yn - y0) + (zn - z0)*(zn - z0));
		//cout<<m<<"  " << "distance=" << distance << "\n";
		//ofile<<m<<"  " << "distance=" << distance << endl;

	}

	double backscatteredcoefficient = backscatter / repeat;
	cout <<  "backscattered coefficient=" << backscatteredcoefficient << "\n";
	ofile << "backscattered coefficient=" << backscatteredcoefficient << endl;

	ofile.close();

	return 0;
}

double compute_alpha(double atomic_number, double energy_keV)
{
	double alpha = 3.4 / 1000 * pow(atomic_number, 0.67) / energy_keV;           //screening factor
	return alpha;
}

