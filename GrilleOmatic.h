#include <array>
#include <vector>
#include <fstream>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <limits>

namespace GrilleOmatic {

	/// GRILLOMATIC MANIFESTO
	/// GRILLOMATIC is an example of implementation of a finite element method to
	/// illustrate WHAT PEOPLE SHOULDN'T DO when trying to solve linar elastodynamics
	/// equations for ultrasound propagation. The goal of the GRILLOMATIC module
	/// is to be able to take into account any homogeneities in the material defined.
	/// The equations are solved using first order finite element method on a regular
	/// grid. Hence, eventhough a lot of standard FEM sub routine can be optimized out,
	/// this type of method WILL LEAD TO POOR ORDER OF CONVERGENCE & POOR MODELLING
	/// PRECISION OF ARBITRARY SHAPED INHOMOGENEITIES AND GEOMETRIES.
	///
	/// In one sentence, the mood around this code is:
	///
	///							GRILLOMATIC, THE MODULE YOU WILL HATE TO LOVE....
	template <
		typename Density, typename Law, typename Source
	> class Model2D {

	public:

		//----------------------------------------------------------------------------//
		/*! \brief Constructor.
		*/
		Model2D(size_t nElemX, Density&& density, Law&& law, Source&& source = Source{});

		/*! \brief Initializing.
		*/
		void init();

		/*! \brief Computing necessary condition for time-scheme stability (CFL condition).
		*/
		double computeStableTimeStep(double CLFTolerance = 0.95);

		/*! \brief Setting time step.
		*/
		void setTimeStep(double ts);

		/*! \brief Moving scheme forward.
		*/
		void forward();

		/*! \brief Swaping solutions.
		*/
		void swap();

		/*! \brief Writing solution.
		*/
		void write(const std::string& filePath) const;
		//----------------------------------------------------------------------------//

	private:

		//----------------------------------------------------------------------------//
		/*! \brief Accessor of the density. The static interface is restricted to returning the
		value of the density at a specific DoF, namely:

		double getDensity(size_t glob) const;
		*/
		Density density_;

		/*! \brief Accessor for constitutive law. The static interface is restricted to computing the
		(symmetric stress tensor) from the symmetric green-lagrange tensor of the deformation at DoF,
		namely:

		std::array<double, 3> applyConstitutiveLaw(size_t g, const std::array<double, 3>& green_lagrange) const;
		*/
		Law law_;

		/*! \brief Accessor for source. If not defined as GrilleOmatic::NoSource, the source should
		satisfies the following static interface:

		std::array<double, 2> eval(size_t g, double time) const;
		*/
		Source source_;
		//----------------------------------------------------------------------------//

	private:

		//----------------------------------------------------------------------------//
		/*! \brief Compute numbering elements.
		*/
		void initNumbering(size_t nx);

		/*! \brief Number of DoF.
		*/
		size_t nDoF_, nDoFX_;

		/*! \brief Number of elements in the grid.
		*/
		size_t nElem_, nElemX_;

		/*! \brief Associated jacobian of the grid.
		*/
		double h_, hx_;
		//----------------------------------------------------------------------------//

	private:

		//----------------------------------------------------------------------------//
		/*! \brief Applying stiffness.
		*/
		void applyStiffness();

		/*! \brief Applying stiffness operator at specific element.
		*/
		void applyLocalStifness(size_t iex, size_t iey);

		/*! \brief Accessing value of mass operator integrating density values.
		*/
		double getMass(size_t glob) const;

		/*! \brief Associted time step and squared time step.
		*/
		double ts_, ts2_;
		size_t nStep_;

		/*! \brief Associated solution at scheme step.
		*/
		std::vector<double> y0_, y1_, y2_;
		//----------------------------------------------------------------------------//
	};

	/// IMPLEMENTATION OF HELPERS.
	///
	///
	///
	///
	///
	///
	///
	///
	///
	///
	template<typename T> struct BinaryWrapper
	{
		BinaryWrapper(T value) : value_(value) {}
		T value_;
	};

	template<typename T> std::ostream& operator<<(std::ostream& out, const BinaryWrapper<T>& wrappedValue)
	{
		const auto chars = reinterpret_cast<const char*>(&wrappedValue.value_);
		for (size_t i = 0, n = sizeof(wrappedValue.value_); i < n; ++i) out.write(chars + (n - i - 1), 1);
		return out;
	}

	template<typename T> BinaryWrapper<T> toBinary(const T& value) { return BinaryWrapper<T>(value); }

	/// GRILLOMATIC IMPLEMENTATION OF TIME SCHEME OPERATIONS.
	///
	///
	///
	///
	///
	///
	///
	///
	///
	///
	template<typename D, typename L, typename S> Model2D<D, L, S>::Model2D(size_t nx, D&& density, L&& law, S&& source)
		: density_(std::move(density)), law_(std::move(law)), source_(std::move(source))
	{
		if (nx <= 1) throw std::exception();
		else
		{
			nElemX_ = nx; nElem_ = nx * nx;
			nDoFX_ = nElemX_ + 1;
			nDoF_ = nDoFX_ * nDoFX_;
		}

		hx_ = 1.0 / nx;
		h_ = hx_ * hx_;
	}

	template<typename D, typename L, typename S> void Model2D<D, L, S>::init()
	{
		y0_.resize(2 * nDoF_, 0.);
		y1_.resize(2 * nDoF_, 0.);
		y2_.resize(2 * nDoF_, 0.);
		nStep_ = 0;
	}

	template<typename D, typename L, typename S> double Model2D<D, L, S>::computeStableTimeStep(double CLFTolerance)
	{
		// Iteration counters.
		const size_t iterMax = 100;
		size_t iter = 0;

		// Algorithm variables.
		const double epsTol = 1e-6;
		double lambda0 = 0., lambda1 = 0., tol;

		// Initializing a random like vector.
		for (size_t i = 0, ni = 2 * nDoF_; i < ni; ++i)
		{
			y1_[i] = i % 7 + (i % 2) * 1e-6 + i % 50;
			y0_[i] = 0.;
		}

		// Iterated power algorithm.
		do {
			lambda0 = lambda1;

			// Applying stiffness.
			applyStiffness();

			// Applying mass inverse.
			for (size_t i = 0; i < nDoF_; ++i)
			{
				const size_t i2 = 2 * i;
				const double mass_rho_inv = 1.0 / (getMass(i) * density_.getDensity(i));
				y1_[i2/**/] = y0_[i2/**/] * mass_rho_inv;
				y1_[i2 + 1] = y0_[i2 + 1] * mass_rho_inv;
			}

			// Computing infinite norm.
			lambda1 = fabs(*std::max_element(y1_.begin(), y1_.end(), [](double a, double b) -> bool { return fabs(a) < fabs(b); }));

			// Scaling and re-init.
			for (size_t i = 0, ni = 2 * nDoF_; i < ni; ++i)
			{
				y1_[i] /= lambda1;
				y0_[i] = 0.0;
			}

			// Computing new tolerance.
			if (lambda0 != 0.0) tol = fabs(lambda1 - lambda0) / fabs(lambda0);
			else tol = std::numeric_limits<double>::max();
			++iter;

		} while (tol > epsTol && iter <= iterMax);

		// Re-init.
		for (size_t i = 0, ni = 2 * nDoF_; i < ni; ++i) y1_[i] = y0_[i] = 0.0;

		return CLFTolerance * 2.0 / sqrt(lambda1);
	}

	template<typename D, typename L, typename S> void Model2D<D, L, S>::setTimeStep(double ts)
	{
		ts_ = ts;
		ts2_ = ts_ * ts_;
	}

	template<typename D, typename L, typename S> void Model2D<D, L, S>::forward()
	{
		const double time = ts_ * nStep_;

		// Applying stiffness.
		applyStiffness();

		// Inverting mass and adding previous solutions.
		for (size_t i = 0; i < nDoF_; ++i)
		{
			const size_t i2 = 2 * i;
			const double mass = getMass(i);
			const double mass_rho_inv = 1.0 / (mass * density_.getDensity(i));
			const auto src = source_.eval(i, time);

			y0_[i2/**/] = ts2_ * (mass * src[0] - y0_[i2/**/]) * mass_rho_inv + 2.0 * y1_[i2/**/] - y2_[i2/**/];
			y0_[i2 + 1] = ts2_ * (mass * src[1] - y0_[i2 + 1]) * mass_rho_inv + 2.0 * y1_[i2 + 1] - y2_[i2 + 1];
		}
	}

	template<typename D, typename L, typename S> void Model2D<D, L, S>::swap()
	{
		y2_.swap(y1_);
		y1_.swap(y0_);
		std::fill(y0_.begin(), y0_.end(), 0.);
		++nStep_;
	}

	template<typename D, typename L, typename S> void Model2D<D, L, S>::write(const std::string& filePath) const
	{
		std::ofstream ofs(filePath, std::ios_base::binary);

		// Writing header.
		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Ondo output" << std::endl;
		ofs << "BINARY" << std::endl;
		ofs << "DATASET STRUCTURED_POINTS" << std::endl;
		ofs << "DIMENSIONS " << nDoFX_ << " " << nDoFX_ << " " << 1 << std::endl;
		ofs << "ORIGIN " << 0. << " " << 0. << " " << 0. << std::endl;
		ofs << "SPACING " << hx_ << " " << hx_ << " " << 0. << std::endl;

		// Writing solution.
		ofs << "POINT_DATA " << nDoF_ << std::endl;
		ofs << "VECTORS SOLUTION double" << std::endl;
		for (size_t i = 0; i < nDoF_; ++i)
				ofs << toBinary<double>(y0_[2 * i]) << toBinary<double>(y0_[2 * i + 1]) << toBinary<double>(0.0);
		ofs << std::endl;
	}

	/// GRILLOMATIC IMPLEMENTATION OF FINITE ELEMENT OPERATOR OPERATIONS
	///
	///
	///
	///
	///
	///
	///
	///
	///
	///
	template<typename D, typename L, typename S> void Model2D<D, L, S>::applyStiffness()
	{
		for (size_t iey = 0; iey < nElemX_; ++iey)
			for (size_t iex = 0; iex < nElemX_; ++iex)
				applyLocalStifness(iex, iey);
	}

	template<typename D, typename L, typename S> void Model2D<D, L, S>::applyLocalStifness(size_t iex, size_t iey)
	{
		using Arr2 = std::array<double, 2>;
		using Arr3 = std::array<double, 3>;
		using Arr8 = std::array<double, 8>;

		// Extracting first global index of element.
		const auto g = iey * nDoFX_ + iex;

		// Extracting local solution.
		const size_t g2 = 2 * g;
		const size_t gg2 = 2 * (g + nDoFX_);

		const Arr8 sol_loc = {
			y1_[g2], y1_[g2 + 1], y1_[g2 + 2], y1_[g2 + 3],
			y1_[gg2], y1_[gg2 + 1], y1_[gg2 + 2], y1_[gg2 + 3] };

		// Applying gradient operator.
		const Arr2 g10 = { sol_loc[2] - sol_loc[0], sol_loc[3] - sol_loc[1] };
		const Arr2 g20 = { sol_loc[4] - sol_loc[0], sol_loc[5] - sol_loc[1] };
		const Arr2 g31 = { sol_loc[6] - sol_loc[2], sol_loc[7] - sol_loc[3] };
		const Arr2 g32 = { sol_loc[6] - sol_loc[4], sol_loc[7] - sol_loc[5] };

		// Computing green-lagrange tensor and applying constitutive law.
		const Arr3 stress0 = law_.applyConstitutiveLaw(g, Arr3{ g10[0], g20[0] + g10[1], g20[1] });
		const Arr3 stress1 = law_.applyConstitutiveLaw(g, Arr3{ g10[0], g31[0] + g10[1], g31[1] });
		const Arr3 stress2 = law_.applyConstitutiveLaw(g, Arr3{ g32[0], g20[0] + g32[1], g20[1] });
		const Arr3 stress3 = law_.applyConstitutiveLaw(g, Arr3{ g32[0], g31[0] + g32[1], g31[1] });

		// Applying transposed gradient operator.
		y0_[g2/**/] -= 0.25 * (stress0[0] + stress0[1] + stress1[0] + stress2[1]);
		y0_[g2 + 1] -= 0.25 * (stress0[1] + stress0[2] + stress1[1] + stress2[2]);
		y0_[g2 + 2] += 0.25 * (stress0[0] + stress1[0] - stress1[1] - stress3[1]);
		y0_[g2 + 3] += 0.25 * (stress0[1] + stress1[1] - stress1[2] - stress3[2]);
		y0_[gg2/**/] += 0.25 * (stress0[1] - stress2[0] + stress2[1] - stress3[0]);
		y0_[gg2 + 1] += 0.25 * (stress0[2] - stress2[1] + stress2[2] - stress3[1]);
		y0_[gg2 + 2] += 0.25 * (stress1[1] + stress2[0] + stress3[0] + stress3[1]);
		y0_[gg2 + 3] += 0.25 * (stress1[2] + stress2[1] + stress3[1] + stress3[2]);
	}

	template<typename D, typename L, typename S> double Model2D<D, L, S>::getMass(size_t g) const
	{
		const size_t gx = g % nDoFX_;
		const size_t gy = g / nDoFX_;

		if (gy == 0 || gy == (nDoFX_ - 1))
		{
			if (gx == 0 || gx == (nDoFX_ - 1)) return 0.25 * h_;
			return 0.5 * h_;
		}
		else
		{
			if (gx == 0 || gx == (nDoFX_ - 1)) return 0.5 * h_;
			return h_;
		}
	}
	///
	///
	///
	///
	///
	///
	///
	///
	///
	///
	/// END GRILLOMATIC.
}
