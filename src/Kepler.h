#pragma once
#include "Eigen/Dense"

#define MU 3.9860044188e14
#define MU_MOON 4.90486959e12
#define MU_SUN 1.327124400189e20
#define AU_TO_M 149597870700.0

struct EmptyType
{

};

struct KeplerElements
{
	// Semi-major axis
	double a;
	// Eccentricity
	double e;
	// Right ascension of ascending node
	double raan;
	// Argument of perigee
	double arg_per;
	// Inclination
	double inc;
	// True anomaly
	double true_anom;
};

template<bool has_vel, bool has_time = false>
struct EulerElements
{
	Eigen::Vector3d pos;
	std::conditional_t<has_time, double, EmptyType> time;
	std::conditional_t<has_vel, Eigen::Vector3d, EmptyType> vel;
};

static KeplerElements euler_to_kepler(const EulerElements<true>& euler)
{
	KeplerElements out;
	Eigen::Vector3d hv = euler.pos.cross(euler.vel);
	double h2 = hv.squaredNorm();
	double h = std::sqrt(h2);

	double r = euler.pos.norm();
	double v = euler.vel.norm();

	double E = v * v / 2.0 - MU / r;
	out.a = -MU / (2.0 * E);
	out.e = std::sqrt(1.0 - h2 / (out.a * MU));

	// Correct quadrant acos
	out.inc = std::acos(hv(3) / h);
	out.raan = std::atan2(hv(1), -hv(2));

	double p = out.a * (1.0 - out.e * out.e);
	double dp = euler.vel.dot(euler.pos);
	out.true_anom = std::atan2(std::sqrt(p / MU) * dp, p - r);

	out.arg_per = std::atan2(euler.pos(2) / std::sin(out.inc),
							 euler.pos(0) * std::cos(out.raan) + euler.pos(1) * std::sin(out.raan));

	return out;

}

// We only support the elliptical case
template<bool has_vel>
static EulerElements<has_vel> kepler_to_euler(const KeplerElements& kepler)
{
	double p = kepler.a * (1.0 - kepler.e * kepler.e);
	double r = p / (1.0 + kepler.e * std::cos(kepler.true_anom));

	EulerElements<has_vel> out;

	out.pos(0) = r * (std::cos(kepler.raan) * std::cos(kepler.arg_per + kepler.true_anom)
			- std::sin(kepler.raan) * std::sin(kepler.arg_per + kepler.true_anom) * std::cos(kepler.inc));
	out.pos(1) = r * (std::sin(kepler.raan) * std::cos(kepler.arg_per + kepler.true_anom)
					  + std::cos(kepler.raan) * std::sin(kepler.arg_per + kepler.true_anom) * std::cos(kepler.inc));
	out.pos(2) = r * (std::sin(kepler.inc) * std::sin(kepler.arg_per + kepler.true_anom));

	if constexpr (has_vel)
	{
		double h = std::sqrt(MU * p);

		double t1 = h * kepler.e / (r * p) * std::sin(kepler.true_anom);
		out.vel(0) = out.pos(0) * t1 - h / r * (std::cos(kepler.raan) * std::sin(kepler.arg_per + kepler.true_anom) +
												std::sin(kepler.raan) * std::cos(kepler.arg_per + kepler.true_anom) *
												std::cos(kepler.inc));
		out.vel(1) = out.pos(1) * t1 - h / r * (std::sin(kepler.raan) * std::sin(kepler.arg_per + kepler.true_anom) -
												std::cos(kepler.raan) * std::cos(kepler.arg_per + kepler.true_anom) *
												std::cos(kepler.inc));
		out.vel(2) = out.pos(2) * t1 + h / r * std::sin(kepler.inc) * std::cos(kepler.arg_per + kepler.true_anom);
	}

	return out;
}