#include "../GrilleOmatic.h"
#include <iostream>

static const size_t nElem = 10;

struct Density {
	double density;
	double getDensity(size_t) const { return density; }
};

struct Law {
	Law(double ym, double pr)
	{
		lambda = ym * pr / ((1.0 + pr) * (1.0 - 2.0 * pr));
		mu = ym / (2.0 * (1.0 + pr));
		mu2 = 2.0 * mu;
	}

	std::array<double, 3> applyConstitutiveLaw(size_t, const std::array<double, 3>& e)
	{
		const double trace = e[0] + e[2];
		return { lambda * trace + mu2 * e[0], mu * e[1], lambda * trace + mu2 * e[2] };
	}

	double lambda, mu, mu2;
};

int main()
{
	GrilleOmatic::Model2D<Density, Law> grillo(nElem, Density{1.0}, Law{1.0, 0.25});
	grillo.init();
	std::cout << grillo.computeStableTimeStep() << std::endl;
}
