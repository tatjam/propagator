#pragma once
#include "Kepler.h"
#include "Eigen/Dense"
#include "vsop87a_large.h"

class Propagator
{
private:

	EulerElements<true> orbiter_elems;
	std::vector<EulerElements<true>> history;

	double t;
	double st;

	Eigen::Vector3d ephemeris_acc;

	using vsop = vsop87a_large;

	// Note, prime is derivatives! pos -> vel  and   vel -> acc
	template<bool eval_time>
	void f(EulerElements<true>& prime, const EulerElements<true>& eval, double t);
	void set_b(EulerElements<true>& b, const EulerElements<true>& prime, const EulerElements<true>& b0, double h);


public:

	bool use_geopotential;
	bool use_ephemerides;

	// tfor: How long to propagate for
	// tstep: Timestep to use during propagation
	// sstep: Saving interval for output vector
	// Returns the saved positions, velocities (if use_vel = true) and time
	// for each sampling position (if use_time = true)
	template<bool use_vel, bool use_time>
	std::vector<EulerElements<use_vel, use_time>> propagate(double tfor, double tstep, double sstep);


	// Start time is seconds since J2000
	void init(double start_time, const EulerElements<true>& initial);

	Propagator();

};

