#pragma once
#include "Kepler.h"
#include <fstream>


// Format is
// has_vel = true, has_time = true 		-> TIME POS_X POS_Y POS_Z VEL_X VEL_Y VEL_Z
// has_vel = false, has_time = true 	-> TIME POS_X POS_Y POS_Z
// has_vel = true, has_time = false 	-> TIME POS_X POS_Y POS_Z
// has_vel = false, has_time = false 	-> POS_X POS_Y POS_Z
// Time in seconds since J2000, position in meters, velocity in meters per second
// (Space separated numbers, records separated by new line)
template<bool has_vel, bool has_time>
void append_table(const std::vector<EulerElements<has_vel, has_time>>& elems, const std::string& file)
{
	std::ofstream f;
	f.open(file, std::ofstream::out | std::ofstream::app);

	for(const EulerElements<has_vel, has_time>& elem : elems)
	{
		if constexpr (has_time)
		{
			f << elem.time << " ";
		}
		f << elem.pos(0) << " " << elem.pos(1) << " " << elem.pos(2);
		if constexpr (has_vel)
		{
			f << " " << elem.vel(0) << " " << elem.vel(1) << " " << elem.vel(2);
		}
		f << std::endl;
	}

	f.close();
}

// Format is
// TIME? a e inc raan arg_per true_anom
// All angle in radians, a in meters
// (Space separated numbers, records separated by new line)
template<bool has_time>
void append_osculating(const std::vector<EulerElements<true, has_time>>& elems, const std::string& file)
{
	std::ofstream f;
	f.open(file, std::ofstream::out | std::ofstream::app);

	for(const EulerElements<true, has_time>& elem : elems)
	{
		if constexpr (has_time)
		{
			f << elem.time << " ";
		}
		KeplerElements osc = euler_to_kepler(elem);
		f << osc.a << " " << osc.e << " " << osc.inc << " " << osc.raan << " " << osc.arg_per << " " << osc.true_anom << std::endl;
	}

	f.close();
}

void clear_file(const std::string& file)
{
	std::ofstream f;
	f.open(file, std::ofstream::out | std::ofstream::trunc);
	f.close();
}