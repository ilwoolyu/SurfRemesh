#pragma once
#include "Geom.h"

class SphericalHarmonics
{
public:
	static void basis(int degree, float *p, float *Y)
	{
		// real spherical harmonics basis functions
		// polar coordinate
		float phi, theta;
		Coordinate::cart2sph(p, &phi, &theta);
		theta = PI / 2 - theta;  // convert to interval [0, PI]
		float *Pm = new float[degree + 1];

		for (int l = 0; l <= degree; l++)
		{
			// legendre part
			Series::legendre(l, cos(theta), Pm);
			float lconstant = sqrt((2 * l + 1) / (4 * PI));

			int center = (l + 1) * (l + 1) - l - 1;

			Y[center] = lconstant * Pm[0];

			// square root of 2
			float sqr2 = sqrt(2.0f);

			for (int m = 1; m <= l; m++)
			{
				float precoeff = lconstant * (float)sqrt(1 / Series::factorial(l + m, l - m + 1));

				if (m % 2 == 1) precoeff = -precoeff;
				Y[center + m] = sqr2 * precoeff * Pm[m] * cos(m * phi);
				Y[center - m] = sqr2 * precoeff * Pm[m] * sin(m * phi);
			}
		}

		delete [] Pm;
	}
};

