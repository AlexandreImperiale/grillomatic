#include "../GrilleOmatic.h"
#include <iostream>
#include <string>

static const size_t nElem = 400;
static const size_t nDoFX = nElem + 1;
static const double step = 1.0 / nElem;

struct Density {
	double density0, density1;
	double getDensity(size_t g) const
	{
		const size_t gy = g / nDoFX;
		if( gy * step < 0.3) return density0;
		else return density1;
	}
};

struct Law {
	Law(double ym, double pr)
	{
		lambda = ym * pr / ((1.0 + pr) * (1.0 - 2.0 * pr));
		mu = ym / (2.0 * (1.0 + pr));
		mu2 = 2.0 * mu;
	}

	std::array<double, 3> applyConstitutiveLaw(size_t, const std::array<double, 3>& e) const
	{
		const double trace = e[0] + e[2];
		return { lambda * trace + mu2 * e[0], mu * e[1], lambda * trace + mu2 * e[2] };
	}

	double lambda, mu, mu2;
};

struct Source {
	std::array<double, 2> eval(size_t g, double t) const
	{
		const size_t gx = g % nDoFX;
		const size_t gy = g / nDoFX;
		const std::array<double, 2> xyz = { step * gx, step * gy };

		const double norm = (xyz[0] - 0.15) * (xyz[0] - 0.15) + (xyz[1] - 0.15) * (xyz[1] - 0.15);
		const double sr = 10000.0;

		const double t0 = 0.1;
		const double tr = 1000.0;
		const double omega = 50.0;

		const double val = std::sin(2.0 * M_PI * omega * (t - t0) ) * std::exp(-(t - t0) * (t - t0) * tr) * std::exp(-norm * sr);

		return {val, 0.};
	}
};

int main()
{
	GrilleOmatic::Model2D<Density, Law, Source> grillo(nElem, Density{1.0, 1.0 / 3.0}, Law{1.0, 0.25});
	grillo.init();
	grillo.setTimeStep(grillo.computeStableTimeStep());

	const size_t nStep = 1000;
	const size_t outputFreq = 10;
	size_t nOutput = 0, step = 0;

	std::cout << "Running model...";

	while(step < nStep)
	{
			grillo.forward();
			if (step % outputFreq == 0)
			{
				grillo.write("./vtk/solution_" + std::to_string(nOutput) + ".vtk");
				++nOutput;
			}
			grillo.swap();
			++step;
	}

	std::cout<< " finished." << std::endl;
}
