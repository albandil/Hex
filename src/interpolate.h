#ifndef HEX_INTERPOLATE
#define HEX_INTERPOLATE

#include <algorithm>

#define INTERPOLATE_AUTO	-1
#define INTERPOLATE_LINEAR	0
#define INTERPOLATE_CSPLINE	1

/**
 * Return interpolated value.
 * \param energies STL-vector-like real array containing ascending (!) energies.
 * \param values STL-vector-like real or complex array containing values.
 * \param energy Where to evaluate interpolated function.
 * \param transform T-to-T' function used when we want to interpolate
 *                  a derived variable rather than the "values" themselves.
 */
template <typename T, typename Transform> auto interpolate(
	rArray energies,
	Array<T> values,
	double energy,
	Transform transform,
	int flag = -1
) -> decltype(transform(T(0)))
{
	// guardians; E_lbound points to first E' ≥ E, whereas E_ubound points
	// to first E" > E; the energies are unique, so E" ≠ E' is equivalent to
	// E' = E ≠ E"
	rArray::const_iterator E_lbound, E_ubound;
	E_lbound = std::lower_bound(energies.begin(), energies.end(), energy);
	E_ubound = std::upper_bound(energies.begin(), energies.end(), energy);
	
	// if THE specified energy IS in "energies", use corresponding value
	if (E_lbound != E_ubound)
		return transform(values[E_lbound - energies.begin()]);
	
	// determine which interpolation rule to use
	//		1) use linear for E < E_ion
	//		2) use single cubic spline for E > E_ion
	if (flag == INTERPOLATE_AUTO)
	{
		if (energy < 1)
			flag = INTERPOLATE_LINEAR;
		else
			flag = INTERPOLATE_CSPLINE;
	}
	
	// we need 2*N points for interpolation
	switch (flag)
	{
		case INTERPOLATE_CSPLINE:	// single cubic spline
			if (E_ubound - energies.begin() < 2 or energies.end() - E_ubound < 2)
			{
				fprintf(stderr, "not enough samples for cubic interpolation (energy: %g)\n", energy);
				flag = INTERPOLATE_LINEAR;	// decrease order
				/* nobreak */
			}
			else
				break;
		default:
		case INTERPOLATE_LINEAR:
			if (E_ubound - energies.begin() < 1 or energies.end() - E_ubound < 1)
			{
				fprintf(stderr, "not enough samples for linear interpolation (energy: %g)\n", energy);
				return transform(0);	// can't interpolate
			}
			break;
	}
	
	// otherwise compute the interpolation
	switch (flag)
	{
		default:
		case INTERPOLATE_LINEAR:
		{
			double energy_left = *(E_ubound - 1);
			double energy_right = *E_ubound;
			auto val_left = transform(values[E_ubound - 1 - energies.begin()]);
			auto val_right = transform(values[E_ubound - energies.begin()]);
			return ((energy_right - energy) * val_left + (energy - energy_left) * val_right) / (energy_right - energy_left);
		}
		case INTERPOLATE_CSPLINE:
		{
			// get 4 interpolated points
			double x0 = *(E_ubound - 2);	auto y0 = transform(values[E_ubound - 2 - energies.begin()]);
			double x1 = *(E_ubound - 1);	auto y1 = transform(values[E_ubound - 1 - energies.begin()]);
			double x2 = *(E_ubound    );	auto y2 = transform(values[E_ubound     - energies.begin()]);
			double x3 = *(E_ubound + 1);	auto y3 = transform(values[E_ubound + 1 - energies.begin()]);
			
			// compute elements of the matrix
			double a11 = (x2-x0)/3.;		double a12 = (x2-x1)/6.;
			double a21 = (x2-x1)/6.;		double a22 = (x3-x1)/3.;
			
			// compute det A
			double det = a11*a22 - a12*a21;
			
			// compute right hand side
			auto b1 = (y2-y1)/(x2-x1) - (y1-y0)/(x1-x0);
			auto b2 = (y3-y2)/(x3-x2) - (y2-y1)/(x2-x1);
			
			// compute solution
			auto d2y_dx2_1 = (b1*a22-b2*a12) / det;
			auto d2y_dx2_2 = (b2*a11-b1*a21) / det;
			
			// evaluate interpolated value
			double x = energy;
			double A = (x2-x)/(x2-x1);
			double B = (x-x1)/(x2-x1);
			double C = (A*A-1)*A*(x2-x1)*(x2-x1)/6.;
			double D = (B*B-1)*B*(x2-x1)*(x2-x1)/6.;
			return A * y1 + B * y2 + C * d2y_dx2_1 + D * d2y_dx2_2;
		}
	}
}

/**
 * Return interpolated value.
 * \param energies STL-vector-like real array containing ascending (!) energies.
 * \param values STL-vector-like real or complex array containing values.
 * \param energy Where to evaluate interpolated function.
 */
template <typename T> T interpolate(rArray energies, Array<T> values, double energy, int flag = -1)
{
	return interpolate(energies, values, energy, [](T x) -> T { return x; }, flag);
}

#endif
