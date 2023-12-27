#pragma once
#include "Eigen/Dense"

#define MU 3.9860044188e14
#define MU_MOON 4.90486959e12
#define MU_SUN 1.327124400189e20
#define AU_TO_M 149597870700.0
// #define J2 1.75553e25 WIKIPEDIA
// #define J2 1.75162e25 APUNTES
#define J2 1.75553e25
#define DEG_TO_RAD 0.01745329
#define RAD_TO_DEG 57.29578


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

template<bool has_time>
static KeplerElements euler_to_kepler(const EulerElements<true,has_time>& euler)
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
	out.inc = std::acos(hv(2) / h);
	out.raan = std::atan2(hv(0), -hv(1));

	double p = out.a * (1.0 - out.e * out.e);
	double dp = euler.vel.dot(euler.pos);
	out.true_anom = std::atan2(std::sqrt(p / MU) * dp, p - r);

	out.arg_per = std::atan2(euler.pos(2) / std::sin(out.inc),
							 euler.pos(0) * std::cos(out.raan) + euler.pos(1) * std::sin(out.raan)) - out.true_anom;

	return out;

}

template<bool has_time>
static KeplerElements euler_to_kepler_circular(const EulerElements<false, has_time>& euler, double a, double i)
{
	KeplerElements out;
	out.e = 0.0;
	out.a = a;
	out.inc = i;

	double r2 = euler.pos.squaredNorm();

	// We can predict the angular momentum from given data, note that
	// h^2 / mu = p = a in circular orbits
	double h = std::sqrt(a * MU);

	// Note that cos(i) = hz / h
	Eigen::Vector3d hv;
	hv(2) = h * std::cos(i);

	// h = sqrt(hv(0)^2 + hv(1)^2 + hv(2)^2)
	// We know that h is perpendicular to the orbit plane, so
	// h dot p = 0
	// which gives
	// hv(0) * p(0) + hv(1) * p(1) + hv(2) * p(2) = 0
	// So we can solve hv(0) and hv(1) from these two equations, giving
	// a quadratic, we take only one of the solutions:
	// (These were obtained using Mathematica)

	double t1 = euler.pos(0) * euler.pos(0) + euler.pos(1) * euler.pos(1);
	double root = std::sqrt(euler.pos(1) * euler.pos(1) * (h * h * t1 - hv(2) * hv(2) * r2));
	hv(0) = (-hv(2) * euler.pos(0) * euler.pos(2) + root) / t1;
	hv(1) = (-hv(2) * euler.pos(1) * euler.pos(1) * euler.pos(2) + euler.pos(0) * root) / (euler.pos(1) * t1);

	// Correct quadrant acos
	out.raan = std::atan2(hv(0), -hv(1));

	out.arg_per = std::atan2(euler.pos(2) / std::sin(out.inc),
							 euler.pos(0) * std::cos(out.raan) + euler.pos(1) * std::sin(out.raan)) - out.true_anom;

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