#include "Propagator.h"
#include "Output.h"
#include <iostream>

#define STEP 100000.0

int main(void)
{
	std::time_t start_t = std::time(nullptr);

	Propagator prop;
	KeplerElements kepler;
	kepler.a = 7258.69e3;
	kepler.e = 0.000001;
	kepler.arg_per = 0.0;
	kepler.inc = 98.9283 * DEG_TO_RAD;

	EulerElements<false> ci;
	ci.pos(0) =  1.56275e6;
	ci.pos(1) =  5.55841e6;
	ci.pos(2) =  4.38619e6;
	KeplerElements b = euler_to_kepler_circular(ci, kepler.a, kepler.inc);

	kepler.raan = b.raan;
	kepler.arg_per = b.arg_per;
	kepler.true_anom = b.true_anom;

	EulerElements<true> start = kepler_to_euler<true>(kepler);
	prop.init(820578120.0, start);


	prop.use_ephemerides = true;
	prop.use_geopotential = true;

	clear_file("out.txt");
	clear_file("out_osc.txt");

	double final_length = 365.0 * 60.0 * 60.0 * 24.0;
	double t = 0.0;
	while(t < final_length)
	{
		auto results = prop.propagate<true, true>(STEP, 10.0, 120.0);
		append_table(results, "out.txt");
		append_osculating(results, "out_osc.txt");

		t += STEP;
		int percent = (int)std::floor(t / final_length * 100.0);
		std::cout << t << " / " << final_length << " (" << percent << "%)" << std::endl;
	}

	std::time_t end_t = std::time(nullptr);
	std::cout << "Simulation took: " << end_t - start_t << "s " << std::endl;

}