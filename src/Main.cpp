#include "Propagator.h"
#include "Output.h"
#include <iostream>

#define STEP 1000.0

int main(void)
{
	Propagator prop;
	KeplerElements kepler;
	kepler.a = 10600e3;
	kepler.e = 0.01;
	kepler.raan = 0.4;
	kepler.arg_per = 0.2;
	kepler.inc = 0.4;
	kepler.true_anom = 0.3;
	EulerElements<true> start = kepler_to_euler<true>(kepler);
	prop.init(0.0, start);

	prop.use_ephemerides = true;
	prop.use_geopotential = true;

	clear_file("out.txt");
	clear_file("out_osc.txt");

	double final_length = 1.0 * 60.0 * 60.0 * 24.0;
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