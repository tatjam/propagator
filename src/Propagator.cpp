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

void Propagator::f(EulerElements<true> &prime, const EulerElements<true> &eval)
{
	prime.pos = eval.vel;
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
	st = sstep;

	double htstep = tstep * 0.5;
	double h = tstep * (1.0 / 6.0);

	EulerElements<true> C1, C2, C3, C4;
	EulerElements<true> b;

	while(propagated < tfor)
	{
		// RK4 propagate
		f(C1, orbiter_elems);
		set_b(b, C1, orbiter_elems, htstep);
		f(C2, b);
		set_b(b, C2, orbiter_elems, htstep);
		f(C3, b);
		set_b(b, C3, orbiter_elems, htstep);
		f(C4, b);

		orbiter_elems.pos += h * (C1.pos + 2.0 * C2.pos + 2.0 * C3.pos + C4.pos);
		orbiter_elems.vel += h * (C1.vel + 2.0 * C2.vel + 2.0 * C3.vel + C4.vel);

		st -= tstep;
		if(st < 0.0)
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
