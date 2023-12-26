#include "Propagator.h"


Propagator::Propagator()
{
	use_geopotential = true;
	use_ephemerides = true;

}

void Propagator::init(double start_time, const EulerElements<true>& initial)
{
	t = start_time;
	st = 0.0;
	orbiter_elems = initial;
}

template<bool eval_time>
void Propagator::f(EulerElements<true> &prime, const EulerElements<true> &eval, double t)
{
	prime.pos = eval.vel;

	// Standard newtonian gravity
	double pnorm = eval.pos.norm();
	double pnorm3 = pnorm * pnorm * pnorm;
	prime.vel = -MU * eval.pos / pnorm3;

	if constexpr (eval_time)
	{
		if(use_ephemerides)
		{
			// x, y, z in AU, J2000 sun centered
			double out_earth[3];
			double out_emb[3];
			double out_moon[3];
			// time is expected in julian centuries
			double ephT = t / 365250.0;
			vsop::getEarth(ephT, out_earth);
			vsop::getEmb(ephT, out_emb);
			vsop::getMoon(out_earth, out_emb, out_moon);

			// (Naming not correct just yet)
			Eigen::Vector3d sun_pos(out_earth[0], out_earth[1], out_earth[2]);
			Eigen::Vector3d moon_pos(out_moon[0], out_moon[1], out_moon[2]);

			// From this it's trivial to obtain positions relative to earth in meters
			moon_pos -= sun_pos;
			moon_pos *= AU_TO_M;
			sun_pos *= -AU_TO_M;

			double lmoon_pos = moon_pos.norm();
			double lsun_pos = sun_pos.norm();

			Eigen::Vector3d sat_to_moon = moon_pos - eval.pos;
			Eigen::Vector3d sat_to_sun = sun_pos - eval.pos;

			double lsat_to_moon = sat_to_moon.norm();
			double lsat_to_sun = sat_to_sun.norm();
			// Newton law on these two bodies
			ephemeris_acc = MU_MOON * sat_to_moon / (lsat_to_moon * lsat_to_moon * lsat_to_moon);
			ephemeris_acc += MU_SUN * sat_to_sun / (lsat_to_sun * lsat_to_sun * lsat_to_sun);

			// The bodies also attract the Earth, include secondary tidal acceleration
			ephemeris_acc -= MU_MOON * moon_pos / (lmoon_pos * lmoon_pos * lmoon_pos);
			ephemeris_acc -= MU_SUN * sun_pos / (lsun_pos * lsun_pos * lsun_pos);
		}
	}

	if(use_ephemerides)
	{
		prime.vel += ephemeris_acc;
	}

	if(use_geopotential)
	{
		// Note that nutation and precession causes earth's rotation axis to change slightly
		// we do implement the simple model of these effects:
		Eigen::Vector3d earth_rot_axis(0, 0, 1);

		//  Note that u dot v = |u| |v| cos(angle)
		// and phi is simply the latitude, so that we have
		double phi = std::acos(earth_rot_axis.dot(eval.pos) / pnorm);
		double sphi = std::sin(phi);
		// Compute J2 effect
		//prime.vel -= MU * J2 * RT * RT / (2.0 * pnorm3) * (1.0 - 3.0 * sphi * sphi) * eval.pos;

	}

}

void
Propagator::set_b(EulerElements<true> &b, const EulerElements<true> &prime, const EulerElements<true> &b0, double h)
{
	b.pos = b0.pos + prime.pos * h;
	b.vel = b0.vel + prime.vel * h;
}

template<bool use_vel, bool use_time>
std::vector<EulerElements<use_vel, use_time>> Propagator::propagate(double tfor, double tstep, double sstep)
{
	std::vector<EulerElements<use_vel, use_time>> out;
	out.reserve((size_t)std::ceil(tfor / sstep));

	double propagated = 0.0;
	st = 0.0;

	double htstep = tstep * 0.5;
	double h = tstep * (1.0 / 6.0);

	EulerElements<true> C1, C2, C3, C4;
	EulerElements<true> b;

	while(propagated < tfor)
	{
		// RK4 propagate
		f<true>(C1, orbiter_elems, t);
		set_b(b, C1, orbiter_elems, htstep);
		f<true>(C2, b, t + htstep);
		set_b(b, C2, orbiter_elems, htstep);
		f<false>(C3, b, t + htstep);
		set_b(b, C3, orbiter_elems, tstep);
		f<true>(C4, b, t + tstep);

		orbiter_elems.pos += h * (C1.pos + 2.0 * C2.pos + 2.0 * C3.pos + C4.pos);
		orbiter_elems.vel += h * (C1.vel + 2.0 * C2.vel + 2.0 * C3.vel + C4.vel);

		st -= tstep;
		if(st <= 0.0)
		{
			EulerElements<use_vel, use_time> sample;
			sample.pos = orbiter_elems.pos;
			if constexpr (use_vel)
			{
				sample.vel = orbiter_elems.vel;
			}
			if constexpr (use_time)
			{
				sample.time = t;
			}
			out.push_back(sample);
			st = sstep;
		}
		t += tstep;
		propagated += tstep;
	}



	return out;
}

// Instantiations
template std::vector<EulerElements<true, true>> Propagator::propagate(double tfor, double tstep, double sstep);
template std::vector<EulerElements<false, true>> Propagator::propagate(double tfor, double tstep, double sstep);
template std::vector<EulerElements<true, false>> Propagator::propagate(double tfor, double tstep, double sstep);
template std::vector<EulerElements<false, false>> Propagator::propagate(double tfor, double tstep, double sstep);
