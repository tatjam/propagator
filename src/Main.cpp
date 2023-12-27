#include "Propagator.h"
#include "Output.h"
#include <iostream>

#define STEP 1000.0

int main(void)
{
	Propagator prop;
	KeplerElements kepler;
	kepler.a = 7633e3;
	kepler.raan = 0.0;
	kepler.arg_per = 0.0;
	kepler.inc = 100.75 * DEG_TO_RAD;
	kepler.true_anom = 0.0;
	EulerElements<true> start = kepler_to_euler<true>(kepler);
	prop.init(0.0, start);

	prop.use_ephemerides = true;
	prop.use_geopotential = true;

	clear_file("out.txt");
	clear_file("out_osc.txt");

	double final_length = 10.0 * 60.0 * 60.0 * 24.0;
	double t = 0.0;
	while(t < final_length)
	{
		auto results = prop.propagate<true, true>(STEP, 100.0, 100.0);
		append_table(results, "out.txt");
		append_osculating(results, "out_osc.txt");

		t += STEP;
		int percent = (int)std::floor(t / final_length * 100.0);
		std::cout << t << " / " << final_length << " (" << percent << "%)" << std::endl;
	}


}